/**
 * @file minuit2_fitter.cpp
 * @brief Minuit2 MIGRAD/MINOS fits for post-sample HBT shape histograms.
 */

#include "hbt/fits/minuit2_fitter.h"

#include "hbt/fits/binned_models.h"
#include "hbt/fits/statistical_analysis.h"

#include <Minuit2/FCNBase.h>
#include <Minuit2/FunctionMinimum.h>
#include <Minuit2/MinosError.h>
#include <Minuit2/MnMigrad.h>
#include <Minuit2/MnMinos.h>
#include <Minuit2/MnUserParameters.h>

#include <array>
#include <cmath>
#include <limits>
#include <optional>
#include <stdexcept>
#include <utility>
#include <vector>

namespace hbt {
namespace {

using ROOT::Minuit2::FunctionMinimum;
using ROOT::Minuit2::MinosError;
using ROOT::Minuit2::MnMigrad;
using ROOT::Minuit2::MnMinos;
using ROOT::Minuit2::MnUserParameters;

/**
 * @brief Deterministic data-derived radius seeds for both model components.
 */
struct RadiusSeeds {
    double gaussian;    ///< Gaussian radius seed.
    double exponential; ///< Exponential radius seed.
};

/**
 * @brief Return an unattempted MIGRAD diagnostic.
 * @return Zeroed stable diagnostic state.
 */
MigradDiagnostic unattempted_migrad() {
    return {false, false, false, false, false, false, 0, std::nullopt};
}

/**
 * @brief Return an unattempted MINOS diagnostic.
 * @return Zeroed stable diagnostic state.
 */
MinosDiagnostic unattempted_minos() {
    return {false, false, false, false, false, false, false, false, false};
}

/**
 * @brief Build an invalid Gaussian result before any minimization.
 * @param reason Exact invalidity reason.
 * @param estimator Statistical objective associated with the empty result.
 * @return Empty numerical result carrying @p reason and @p estimator.
 */
GaussianFitResult invalid_gaussian(
    FitFailureReason reason,
    FitEstimator estimator
) {
    const MigradDiagnostic empty = unattempted_migrad();
    std::array<MigradDiagnostic, GaussianFitResult::kStartCount> starts{};
    starts.fill(empty);
    return {
        false,
        reason,
        estimator,
        starts,
        0U,
        0U,
        std::nullopt,
        empty,
        unattempted_minos(),
        unattempted_minos(),
        std::nullopt,
        std::nullopt,
        std::nullopt,
        {}
    };
}

/**
 * @brief Build an invalid mixed result before any minimization.
 * @param reason Exact invalidity reason.
 * @param estimator Statistical objective associated with the empty result.
 * @return Empty numerical result carrying @p reason and @p estimator.
 */
MixedFitResult invalid_mixed(
    FitFailureReason reason,
    FitEstimator estimator
) {
    const MigradDiagnostic empty = unattempted_migrad();
    std::array<MigradDiagnostic, MixedFitResult::kCoreStartCount> starts{};
    starts.fill(empty);
    std::array<
        MixedStartEndpointDiagnostic,
        MixedFitResult::kCoreStartCount
    > start_endpoints{};
    return {
        false,
        reason,
        estimator,
        starts,
        start_endpoints,
        0U,
        0U,
        0U,
        std::nullopt,
        empty,
        empty,
        unattempted_minos(),
        unattempted_minos(),
        unattempted_minos(),
        unattempted_minos(),
        std::nullopt,
        std::nullopt,
        std::nullopt,
        std::nullopt,
        std::nullopt,
        {}
    };
}

/**
 * @brief Derive deterministic radius seeds from selected raw-bin moments.
 * @param family OSL or radial physical model family.
 * @param bins Slot-major raw uint64_t histogram storage.
 * @param offset First raw counter for the logical histogram.
 * @param binning Validated uniform histogram binning.
 * @param region Selected contiguous statistical region.
 * @return Gaussian and exponential radius seeds, or std::nullopt if invalid.
 *
 * OSL uses E[x^2]=2R^2 for the Gaussian and E[x]=R for the
 * exponential. Radial uses E[r^2]=6R^2 and E[r]=3R respectively.
 */
std::optional<RadiusSeeds> derive_radius_seeds(
    FitObservableFamily family,
    const std::vector<std::uint64_t>& bins,
    std::size_t offset,
    const HistogramBinningConfig& binning,
    const StatisticalRegion& region
) {
    if (offset > bins.size() ||
        binning.nbins > bins.size() - offset ||
        region.last_bin >= binning.nbins ||
        region.selected_count == 0U) {
        return std::nullopt;
    }

    double weighted = 0.0;
    double weighted_square = 0.0;
    for (std::size_t bin = region.first_bin;
         bin <= region.last_bin;
         ++bin) {
        const double center = histogram_bin_center(binning, bin);
        const double count = static_cast<double>(bins[offset + bin]);
        weighted += count * center;
        weighted_square += count * center * center;
    }

    const double n = static_cast<double>(region.selected_count);
    const double mean = weighted / n;
    const double mean_square = weighted_square / n;
    double gaussian = 0.0;
    double exponential = 0.0;

    switch (family) {
        case FitObservableFamily::OSL:
            gaussian = std::sqrt(mean_square / 2.0);
            exponential = mean;
            break;
        case FitObservableFamily::Radial:
            gaussian = std::sqrt(mean_square / 6.0);
            exponential = mean / 3.0;
            break;
    }

    if (!std::isfinite(gaussian) || !std::isfinite(exponential) ||
        gaussian <= 0.0 || exponential <= 0.0) {
        return std::nullopt;
    }
    return RadiusSeeds{gaussian, exponential};
}

/**
 * @brief Exact integral of the unnormalized Gaussian over one fit region.
 */
double gaussian_region_integral(
    FitObservableFamily family,
    const HistogramBinningConfig& binning,
    const StatisticalRegion& region,
    double radius
) {
    const double integral = gaussian_component_integral(
        family,
        histogram_bin_lower_edge(binning, region.first_bin),
        histogram_bin_upper_edge(binning, region.last_bin),
        radius
    );
    if (!std::isfinite(integral) || integral <= 0.0) {
        throw std::invalid_argument(
            "HBT Gaussian fit: invalid integrated Gaussian normalization"
        );
    }
    return integral;
}

/**
 * @brief Minuit2 objective for a Gaussian with free positive amplitude.
 *
 * All references are borrowed only for the synchronous lifetime of one fit.
 * No reference or pointer escapes the fit call.
 */
class GaussianObjective final : public ROOT::Minuit2::FCNBase {
public:
    /**
     * @brief Construct a synchronous Gaussian estimator objective.
     * @param family OSL or radial physical family.
     * @param bins Borrowed raw uint64_t count storage.
     * @param offset First logical histogram counter.
     * @param binning Borrowed validated uniform binning.
     * @param region Selected contiguous statistical region, copied by value.
     * @param estimator Statistical objective evaluated for this fit only.
     */
    GaussianObjective(
        FitObservableFamily family,
        const std::vector<std::uint64_t>& bins,
        std::size_t offset,
        const HistogramBinningConfig& binning,
        StatisticalRegion region,
        FitEstimator estimator
    )
        : family_(family),
          bins_(bins),
          offset_(offset),
          binning_(binning),
          region_(region),
          estimator_(estimator) {}

    /**
     * @brief Evaluate the configured objective from log-radius and log-amplitude.
     * @param parameters Minuit external parameters: log(R), log(A_G), where
     *        A_G multiplies the unnormalized exact-bin Gaussian integral.
     * @return Finite deviance, or the largest finite double for an invalid
     *         model evaluation. Returned minima are validated explicitly.
     */
    double operator()(const std::vector<double>& parameters) const override {
        if (parameters.size() != 2U) {
            return std::numeric_limits<double>::max();
        }
        const double radius = std::exp(parameters[0]);
        const double amplitude = std::exp(parameters[1]);
        try {
            const std::vector<double> log_expected =
                gaussian_bin_log_expected_counts(
                    family_,
                    binning_,
                    region_,
                    radius,
                    amplitude
                );
            switch (estimator_) {
                case FitEstimator::Poisson:
                    return binned_poisson_deviance_from_log_expected(
                        bins_, offset_, region_, log_expected
                    );
                case FitEstimator::Neyman:
                    return binned_neyman_chi_square_from_log_expected(
                        bins_, offset_, region_, log_expected
                    );
                case FitEstimator::Pearson:
                    return binned_pearson_chi_square_from_log_expected(
                        bins_, offset_, region_, log_expected
                    );
            }
        } catch (const std::exception&) {
            return std::numeric_limits<double>::max();
        }
        return std::numeric_limits<double>::max();
    }

    /**
     * @brief Return Minuit error definition for the configured objective.
     * @return Exactly 1.0 for deviance or chi-square MINOS intervals.
     */
    double Up() const override {
        return 1.0;
    }

private:
    FitObservableFamily family_; ///< Physical model family.
    const std::vector<std::uint64_t>& bins_; ///< Borrowed raw counts.
    std::size_t offset_; ///< First logical raw counter.
    const HistogramBinningConfig& binning_; ///< Borrowed binning metadata.
    StatisticalRegion region_; ///< Selected region copied by value.
    FitEstimator estimator_; ///< Objective used only by this independent fit.
};

/**
 * @brief Three-parameter Neyman objective for the mixed shape.
 *
 * Parameter order is log(R_core), log(R_tail), f_core. For every tuple the
 * exact unnormalized mixed bin integrals p_i are evaluated and the positive
 * amplitude A is recalculated analytically at the Neyman minimum before the
 * objective is returned. Borrowed references are retained only during one
 * synchronous fit invocation.
 */
class MixedObjective final : public ROOT::Minuit2::FCNBase {
public:
    /**
     * @brief Construct a synchronous mixed-model objective.
     * @param family OSL or radial physical family.
     * @param bins Borrowed raw uint64_t count storage.
     * @param offset First logical histogram counter.
     * @param binning Borrowed validated uniform binning.
     * @param region Selected contiguous statistical region, copied by value.
     */
    MixedObjective(
        FitObservableFamily family,
        const std::vector<std::uint64_t>& bins,
        std::size_t offset,
        const HistogramBinningConfig& binning,
        StatisticalRegion region
    )
        : family_(family),
          bins_(bins),
          offset_(offset),
          binning_(binning),
          region_(region) {}

    /**
     * @brief Evaluate the configured objective from mixed-model parameters.
     * @param parameters Minuit external parameters in fixed model order.
     * @return Finite objective value, or the largest finite double for an
     *         invalid model/objective evaluation.
     */
    double operator()(const std::vector<double>& parameters) const override {
        if (parameters.size() != 3U) {
            return std::numeric_limits<double>::max();
        }
        const double core_radius = std::exp(parameters[0]);
        const double tail_radius = std::exp(parameters[1]);
        try {
            const std::vector<double> bin_integrals = mixed_bin_integrals(
                family_,
                binning_,
                region_,
                core_radius,
                tail_radius,
                parameters[2]
            );
            const double amplitude = neyman_optimal_mixed_amplitude(
                bins_, offset_, region_, bin_integrals
            );
            return binned_neyman_chi_square_from_integrals(
                bins_, offset_, region_, bin_integrals, amplitude
            );
        } catch (const std::exception&) {
            return std::numeric_limits<double>::max();
        }
    }

    /**
     * @brief Return the MINOS error definition for Neyman chi-square.
     * @return Exactly 1.0, corresponding to Delta(chi-square)=1.
     */
    double Up() const override {
        return 1.0;
    }

private:
    FitObservableFamily family_; ///< Physical model family.
    const std::vector<std::uint64_t>& bins_; ///< Borrowed raw counts.
    std::size_t offset_; ///< First logical raw counter.
    const HistogramBinningConfig& binning_; ///< Borrowed binning metadata.
    StatisticalRegion region_; ///< Selected region copied by value.
};

/**
 * @brief Four-parameter Neyman objective used only for the final A profile.
 *
 * The production mixed search never fits A numerically. This explicit
 * log-amplitude coordinate exists solely so MINOS can profile A while varying
 * R_core, R_tail, and f_core after the analytically profiled 3D minimum has
 * been selected.
 */
class MixedAmplitudeProfileObjective final : public ROOT::Minuit2::FCNBase {
public:
    MixedAmplitudeProfileObjective(
        FitObservableFamily family,
        const std::vector<std::uint64_t>& bins,
        std::size_t offset,
        const HistogramBinningConfig& binning,
        StatisticalRegion region
    )
        : family_(family),
          bins_(bins),
          offset_(offset),
          binning_(binning),
          region_(region) {}

    double operator()(const std::vector<double>& parameters) const override {
        if (parameters.size() != 4U) {
            return std::numeric_limits<double>::max();
        }
        const double core_radius = std::exp(parameters[0]);
        const double tail_radius = std::exp(parameters[1]);
        const double amplitude = std::exp(parameters[3]);
        try {
            const std::vector<double> bin_integrals = mixed_bin_integrals(
                family_,
                binning_,
                region_,
                core_radius,
                tail_radius,
                parameters[2]
            );
            return binned_neyman_chi_square_from_integrals(
                bins_, offset_, region_, bin_integrals, amplitude
            );
        } catch (const std::exception&) {
            return std::numeric_limits<double>::max();
        }
    }

    double Up() const override {
        return 1.0;
    }

private:
    FitObservableFamily family_;
    const std::vector<std::uint64_t>& bins_;
    std::size_t offset_;
    const HistogramBinningConfig& binning_;
    StatisticalRegion region_;
};

/**
 * @brief Convert one FunctionMinimum into stable MIGRAD diagnostics.
 * @param minimum Completed MIGRAD result.
 * @param objective_failure Whether the returned minimum is unevaluable.
 * @return Stable diagnostic without retaining Minuit state.
 */
MigradDiagnostic migrad_diagnostic(
    const FunctionMinimum& minimum,
    bool objective_failure
) {
    const double q_min = minimum.Fval();
    return {
        true,
        minimum.IsValid(),
        minimum.UserState().HasCovariance(),
        minimum.HasReachedCallLimit(),
        minimum.IsAboveMaxEdm(),
        objective_failure,
        minimum.NFcn(),
        std::isfinite(q_min) && !objective_failure
            ? std::optional<double>{q_min}
            : std::nullopt
    };
}

/**
 * @brief Convert one MinosError into stable side diagnostics.
 * @param error Completed MINOS result for one parameter.
 * @return Stable diagnostic without retaining Minuit state.
 */
MinosDiagnostic minos_diagnostic(const MinosError& error) {
    return {
        true,
        error.LowerValid(),
        error.UpperValid(),
        error.AtLowerLimit(),
        error.AtUpperLimit(),
        error.AtLowerMaxFcn(),
        error.AtUpperMaxFcn(),
        error.LowerNewMin(),
        error.UpperNewMin()
    };
}

/**
 * @brief Test whether one returned Gaussian minimum has a valid objective.
 * @param objective Borrowed synchronous objective.
 * @param minimum Completed MIGRAD result.
 * @return true when the returned parameter state evaluates to a finite value.
 */
bool gaussian_minimum_is_evaluable(
    const GaussianObjective& objective,
    const FunctionMinimum& minimum
) {
    try {
        const std::vector<double> parameters{
            minimum.UserState().Value(0U),
            minimum.UserState().Value(1U)
        };
        const double value = objective(parameters);
        return std::isfinite(value) &&
            value != std::numeric_limits<double>::max();
    } catch (const std::exception&) {
        return false;
    }
}

/**
 * @brief Test whether one returned mixed minimum has a valid objective.
 * @param objective Borrowed synchronous objective.
 * @param minimum Completed MIGRAD result.
 * @return true when the returned parameter state evaluates to a finite value.
 */
bool mixed_minimum_is_evaluable(
    const MixedObjective& objective,
    const FunctionMinimum& minimum
) {
    try {
        const std::vector<double> parameters{
            minimum.UserState().Value(0U),
            minimum.UserState().Value(1U),
            minimum.UserState().Value(2U)
        };
        const double value = objective(parameters);
        return std::isfinite(value) &&
            value != std::numeric_limits<double>::max();
    } catch (const std::exception&) {
        return false;
    }
}

/**
 * @brief Test whether one explicit-amplitude profile minimum is evaluable.
 */
bool mixed_amplitude_profile_minimum_is_evaluable(
    const MixedAmplitudeProfileObjective& objective,
    const FunctionMinimum& minimum
) {
    try {
        const std::vector<double> parameters{
            minimum.UserState().Value(0U),
            minimum.UserState().Value(1U),
            minimum.UserState().Value(2U),
            minimum.UserState().Value(3U)
        };
        const double value = objective(parameters);
        return std::isfinite(value) &&
            value != std::numeric_limits<double>::max();
    } catch (const std::exception&) {
        return false;
    }
}

/**
 * @brief Transform MINOS log-radius errors into physical R-space distances.
 * @param log_radius Central fitted log radius.
 * @param error Fully valid MINOS result for log radius.
 * @return Physical radius and non-negative asymmetric physical errors.
 * @throws std::invalid_argument If the transformed endpoints are invalid.
 */
FitParameterEstimate physical_radius_estimate(
    double log_radius,
    const MinosError& error
) {
    const double radius = std::exp(log_radius);
    const double lower = std::exp(log_radius + error.Lower());
    const double upper = std::exp(log_radius + error.Upper());
    if (!std::isfinite(radius) || !std::isfinite(lower) ||
        !std::isfinite(upper) || radius <= 0.0 || lower <= 0.0 ||
        upper <= 0.0 || lower > radius || upper < radius) {
        throw std::invalid_argument(
            "HBT analysis MINOS: invalid physical radius interval"
        );
    }
    return {radius, radius - lower, upper - radius};
}

/**
 * @brief Convert a direct physical-parameter MINOS interval to distances.
 * @param value Central fitted physical value.
 * @param error Fully valid MINOS result for the same parameter.
 * @return Physical value and non-negative asymmetric distances.
 * @throws std::invalid_argument If the transformed interval is invalid.
 */
FitParameterEstimate direct_parameter_estimate(
    double value,
    const MinosError& error
) {
    const double lower = value + error.Lower();
    const double upper = value + error.Upper();
    if (!std::isfinite(value) || !std::isfinite(lower) ||
        !std::isfinite(upper) || lower > value || upper < value) {
        throw std::invalid_argument(
            "HBT analysis MINOS: invalid direct parameter interval"
        );
    }
    return {value, value - lower, upper - value};
}

/**
 * @brief Execute one independent pure-Gaussian MIGRAD start.
 * @param objective Borrowed synchronous Minuit objective.
 * @param radius_seed Strictly positive physical radius seed.
 * @param amplitude_seed Strictly positive Gaussian amplitude seed.
 * @return FunctionMinimum plus its stable diagnostic.
 */
std::pair<FunctionMinimum, MigradDiagnostic> run_gaussian_migrad(
    GaussianObjective& objective,
    double radius_seed,
    double amplitude_seed
) {
    MnUserParameters parameters;
    parameters.Add("log_r", std::log(radius_seed), 0.1);
    parameters.Add("log_A_G", std::log(amplitude_seed), 0.1);
    MnMigrad migrad(objective, parameters);
    FunctionMinimum minimum = migrad();
    const bool objective_failure =
        !gaussian_minimum_is_evaluable(objective, minimum);
    return {
        minimum,
        migrad_diagnostic(minimum, objective_failure)
    };
}

/**
 * @brief One completed deterministic Gaussian start retained until selection.
 */
struct GaussianStartOutcome {
    FunctionMinimum minimum;      ///< Complete Minuit minimum for later MINOS.
    MigradDiagnostic diagnostic;  ///< Stable serialized start diagnostic.
};

/**
 * @brief Execute one independent mixed-model MIGRAD start.
 * @param objective Borrowed synchronous Minuit objective.
 * @param core_seed Strictly positive Gaussian core radius seed.
 * @param tail_seed Strictly positive exponential tail radius seed.
 * @param core_fraction_seed Initial Gaussian mixing coefficient in (0,1).
 * @return FunctionMinimum plus its stable diagnostic.
 *
 * Numerical parameter steps are Minuit controls, not scientific configuration.
 * Scientific start values are supplied explicitly by the Cartesian-product
 * caller.
 */
std::pair<FunctionMinimum, MigradDiagnostic> run_mixed_migrad(
    MixedObjective& objective,
    double core_seed,
    double tail_seed,
    double core_fraction_seed
) {
    MnUserParameters parameters;
    parameters.Add("log_r_core", std::log(core_seed), 0.1);
    parameters.Add("log_r_tail", std::log(tail_seed), 0.1);
    parameters.Add("f_core", core_fraction_seed, 0.05, 0.0, 1.0);
    MnMigrad migrad(objective, parameters);
    FunctionMinimum minimum = migrad();
    const bool objective_failure =
        !mixed_minimum_is_evaluable(objective, minimum);
    return {
        minimum,
        migrad_diagnostic(minimum, objective_failure)
    };
}

/**
 * @brief Build the explicit four-parameter state used only for the A profile.
 * @param objective Explicit-amplitude Neyman objective.
 * @param core_radius Final positive R_core from the 3D profiled search.
 * @param tail_radius Final positive R_tail from the 3D profiled search.
 * @param core_fraction Final mixing coefficient in (0,1).
 * @param amplitude Analytic positive A at the same 3D minimum.
 * @return Fresh 4D MIGRAD minimum plus stable diagnostics.
 *
 * The starting point is already the exact conditional minimum in A. This
 * auxiliary 4D state is never used for basin selection or central shape
 * parameters; it exists only so MINOS can profile the final amplitude.
 */
std::pair<FunctionMinimum, MigradDiagnostic> run_mixed_amplitude_profile_migrad(
    MixedAmplitudeProfileObjective& objective,
    double core_radius,
    double tail_radius,
    double core_fraction,
    double amplitude
) {
    MnUserParameters parameters;
    parameters.Add("log_r_core", std::log(core_radius), 0.1);
    parameters.Add("log_r_tail", std::log(tail_radius), 0.1);
    parameters.Add("f_core", core_fraction, 0.05, 0.0, 1.0);
    parameters.Add("log_A", std::log(amplitude), 0.1);
    MnMigrad migrad(objective, parameters);
    FunctionMinimum minimum = migrad();
    const bool objective_failure =
        !mixed_amplitude_profile_minimum_is_evaluable(objective, minimum);
    return {
        minimum,
        migrad_diagnostic(minimum, objective_failure)
    };
}

/**
 * @brief One completed deterministic mixed start retained until selection.
 */
struct MixedStartOutcome {
    FunctionMinimum minimum;      ///< Complete Minuit minimum for later MINOS.
    MigradDiagnostic diagnostic;  ///< Stable serialized start diagnostic.
    MixedBasinPoint endpoint;     ///< External endpoint for basin grouping.
};

/**
 * @brief Convert one mixed terminal point to finite physical diagnostics.
 * @param endpoint Final Minuit external coordinates for one deterministic start.
 * @return Optional physical R_core, R_tail and f_core values.
 *
 * The conversion is diagnostic only. It is performed for every attempted start,
 * including starts that later fail the MIGRAD acceptance contract. Non-finite or
 * non-positive radius transforms and non-finite fractions are left empty.
 */
MixedStartEndpointDiagnostic mixed_start_endpoint_diagnostic(
    const MixedBasinPoint& endpoint
) {
    MixedStartEndpointDiagnostic diagnostic{};
    const double core_radius = std::exp(endpoint.log_core_radius);
    const double tail_radius = std::exp(endpoint.log_tail_radius);
    if (std::isfinite(core_radius) && core_radius > 0.0) {
        diagnostic.core_radius = core_radius;
    }
    if (std::isfinite(tail_radius) && tail_radius > 0.0) {
        diagnostic.tail_radius = tail_radius;
    }
    if (std::isfinite(endpoint.core_fraction)) {
        diagnostic.core_fraction = endpoint.core_fraction;
    }
    return diagnostic;
}


}  // namespace

GaussianFitResult fit_gaussian_model(
    FitObservableFamily family,
    const std::vector<std::uint64_t>& bins,
    std::size_t offset,
    const HistogramBinningConfig& binning,
    const StatisticalRegion& region,
    FitEstimator estimator,
    double half_maximum_seed
) {
    const std::size_t selected_bins =
        region.last_bin - region.first_bin + 1U;
    if (selected_bins < 3U || region.selected_count == 0U) {
        return invalid_gaussian(FitFailureReason::InsufficientBins, estimator);
    }

    const std::optional<RadiusSeeds> moment_seeds = derive_radius_seeds(
        family,
        bins,
        offset,
        binning,
        region
    );
    if (!moment_seeds.has_value()) {
        return invalid_gaussian(FitFailureReason::InvalidMomentSeed, estimator);
    }
    if (!std::isfinite(half_maximum_seed) || half_maximum_seed <= 0.0) {
        return invalid_gaussian(
            FitFailureReason::InvalidHalfMaximumSeed,
            estimator
        );
    }

    const std::array<double, GaussianFitResult::kStartCount> radius_seeds{
        moment_seeds->gaussian,
        half_maximum_seed
    };
    GaussianFitResult result = invalid_gaussian(
        FitFailureReason::MigradInvalid,
        estimator
    );
    std::vector<GaussianStartOutcome> outcomes;
    outcomes.reserve(GaussianFitResult::kStartCount);
    std::vector<std::size_t> valid_indices;
    valid_indices.reserve(GaussianFitResult::kStartCount);

    for (std::size_t index = 0U;
         index < GaussianFitResult::kStartCount;
         ++index) {
        GaussianObjective objective(
            family, bins, offset, binning, region, estimator
        );
        double amplitude_seed = 0.0;
        try {
            amplitude_seed = 1.0 / gaussian_region_integral(
                family, binning, region, radius_seeds[index]
            );
        } catch (const std::exception&) {
            result.failure_reason = FitFailureReason::ObjectiveEvaluation;
            return result;
        }
        auto start = run_gaussian_migrad(
            objective, radius_seeds[index], amplitude_seed
        );
        result.starts[index] = start.second;
        ++result.starts_attempted;
        outcomes.push_back({std::move(start.first), start.second});
        if (fit_failure_from_migrad(start.second) == FitFailureReason::None) {
            valid_indices.push_back(index);
        }
    }
    result.valid_starts = valid_indices.size();
    if (valid_indices.empty()) {
        for (const MigradDiagnostic& diagnostic : result.starts) {
            const FitFailureReason reason = fit_failure_from_migrad(diagnostic);
            if (reason == FitFailureReason::ObjectiveEvaluation) {
                result.failure_reason = reason;
                return result;
            }
        }
        result.failure_reason = fit_failure_from_migrad(result.starts[0U]);
        return result;
    }

    std::size_t selected_index = valid_indices.front();
    for (const std::size_t index : valid_indices) {
        if (outcomes[index].diagnostic.q_min.value() <
            outcomes[selected_index].diagnostic.q_min.value()) {
            selected_index = index;
        }
    }
    result.selected_start = selected_index;
    result.migrad = outcomes[selected_index].diagnostic;
    result.q_min = result.migrad.q_min;

    FunctionMinimum& selected = outcomes[selected_index].minimum;
    const double log_radius = selected.UserState().Value(0U);
    const double log_amplitude = selected.UserState().Value(1U);
    const double radius = std::exp(log_radius);
    const double amplitude = std::exp(log_amplitude);
    if (!std::isfinite(log_radius) || !std::isfinite(log_amplitude) ||
        !std::isfinite(radius) || !std::isfinite(amplitude) ||
        radius <= 0.0 || amplitude <= 0.0 || !result.q_min.has_value()) {
        result.failure_reason = FitFailureReason::NonFiniteMinimum;
        return result;
    }

    GaussianObjective selected_objective(
        family, bins, offset, binning, region, estimator
    );
    MnMinos minos(selected_objective, selected);
    const MinosError radius_error = minos.Minos(0U);
    result.minos_radius = minos_diagnostic(radius_error);
    const FitFailureReason radius_minos_reason = fit_failure_from_minos(
        result.minos_radius,
        false
    );
    if (radius_minos_reason != FitFailureReason::None) {
        result.failure_reason = radius_minos_reason;
        return result;
    }
    const MinosError amplitude_error = minos.Minos(1U);
    result.minos_amplitude = minos_diagnostic(amplitude_error);
    const FitFailureReason amplitude_minos_reason = fit_failure_from_minos(
        result.minos_amplitude,
        false
    );
    if (amplitude_minos_reason != FitFailureReason::None) {
        result.failure_reason = amplitude_minos_reason;
        return result;
    }

    try {
        result.radius = physical_radius_estimate(log_radius, radius_error);
        result.amplitude = physical_radius_estimate(
            log_amplitude, amplitude_error
        );
        const double bin_width = 1.0 / binning.inverse_bin_width;
        const double log_density_scale = std::log(amplitude) -
            std::log(bin_width);
        result.fitted_pdf.clear();
        result.fitted_pdf.reserve(
            region.last_bin - region.first_bin + 1U
        );
        const double log_min_double = std::log(
            std::numeric_limits<double>::denorm_min()
        );
        const double log_max_double = std::log(
            std::numeric_limits<double>::max()
        );
        for (std::size_t bin = region.first_bin; bin <= region.last_bin; ++bin) {
            const double log_density = log_density_scale +
                gaussian_component_log_integral(
                    family,
                    histogram_bin_lower_edge(binning, bin),
                    histogram_bin_upper_edge(binning, bin),
                    radius
                );
            if (!std::isfinite(log_density) || log_density > log_max_double) {
                throw std::invalid_argument(
                    "HBT Gaussian fit: invalid fitted density"
                );
            }
            result.fitted_pdf.push_back(
                log_density < log_min_double ? 0.0 : std::exp(log_density)
            );
        }
    } catch (const std::exception&) {
        result.radius = std::nullopt;
        result.amplitude = std::nullopt;
        result.fitted_pdf.clear();
        result.failure_reason = FitFailureReason::NonFiniteMinimum;
        return result;
    }

    result.fully_valid = true;
    result.failure_reason = FitFailureReason::None;
    return result;
}

MixedFitResult fit_mixed_model(
    FitObservableFamily family,
    const std::vector<std::uint64_t>& bins,
    std::size_t offset,
    const HistogramBinningConfig& binning,
    const StatisticalRegion& region,
    FitEstimator estimator,
    const GaussianFitResult& gaussian_result,
    double half_maximum_seed,
    MixedCoreFractionPolicy core_fraction_policy
) {
    if (estimator != FitEstimator::Neyman) {
        return invalid_mixed(FitFailureReason::NotApplicable, estimator);
    }
    const std::size_t selected_bins =
        region.last_bin - region.first_bin + 1U;
    if (selected_bins < 5U || region.selected_count == 0U) {
        return invalid_mixed(FitFailureReason::InsufficientBins, estimator);
    }
    if (gaussian_result.estimator != estimator ||
        !gaussian_result.fully_valid || !gaussian_result.radius.has_value() ||
        !std::isfinite(gaussian_result.radius->value) ||
        gaussian_result.radius->value <= 0.0) {
        return invalid_mixed(
            FitFailureReason::InvalidGaussianCoreAnchor,
            estimator
        );
    }
    if (!std::isfinite(half_maximum_seed) || half_maximum_seed <= 0.0) {
        return invalid_mixed(
            FitFailureReason::InvalidHalfMaximumSeed,
            estimator
        );
    }

    const std::optional<RadiusSeeds> moment_seeds = derive_radius_seeds(
        family,
        bins,
        offset,
        binning,
        region
    );
    if (!moment_seeds.has_value()) {
        return invalid_mixed(FitFailureReason::InvalidMomentSeed, estimator);
    }

    const double gaussian_seed = gaussian_result.radius->value;
    const double tail_moment_seed = moment_seeds->exponential;
    const std::array<double, 4U> core_seeds{
        gaussian_seed,
        0.5 * half_maximum_seed,
        half_maximum_seed,
        2.0 * half_maximum_seed
    };
    const std::array<double, 3U> tail_scales{0.5, 1.0, 2.0};
    const std::array<double, 3U> fractions{0.25, 0.50, 0.75};

    MixedFitResult result = invalid_mixed(
        FitFailureReason::MigradInvalid,
        estimator
    );
    std::vector<MixedStartOutcome> outcomes;
    outcomes.reserve(MixedFitResult::kCoreStartCount);
    std::vector<std::size_t> valid_indices;
    valid_indices.reserve(MixedFitResult::kCoreStartCount);

    std::size_t index = 0U;
    for (const double core_seed : core_seeds) {
        for (const double tail_scale : tail_scales) {
            for (const double fraction : fractions) {
                MixedObjective objective(
                    family, bins, offset, binning, region
                );
                auto start = run_mixed_migrad(
                    objective,
                    core_seed,
                    tail_moment_seed * tail_scale,
                    fraction
                );
                result.starts[index] = start.second;
                ++result.starts_attempted;

                MixedBasinPoint endpoint{
                    start.first.UserState().Value(0U),
                    start.first.UserState().Value(1U),
                    start.first.UserState().Value(2U)
                };
                result.start_endpoints[index] =
                    mixed_start_endpoint_diagnostic(endpoint);
                outcomes.push_back({
                    std::move(start.first), start.second, endpoint
                });
                if (fit_failure_from_migrad(start.second) ==
                    FitFailureReason::None) {
                    valid_indices.push_back(index);
                }
                ++index;
            }
        }
    }
    result.valid_starts = valid_indices.size();
    if (valid_indices.empty()) {
        for (const MigradDiagnostic& diagnostic : result.starts) {
            const FitFailureReason reason = fit_failure_from_migrad(diagnostic);
            if (reason == FitFailureReason::ObjectiveEvaluation) {
                result.failure_reason = reason;
                return result;
            }
        }
        result.failure_reason = fit_failure_from_migrad(result.starts[0U]);
        return result;
    }

    std::vector<MixedBasinPoint> endpoints;
    std::vector<double> q_values;
    endpoints.reserve(outcomes.size());
    q_values.reserve(outcomes.size());
    for (const MixedStartOutcome& outcome : outcomes) {
        endpoints.push_back(outcome.endpoint);
        q_values.push_back(
            outcome.diagnostic.q_min.value_or(
                std::numeric_limits<double>::infinity()
            )
        );
    }
    const std::vector<std::size_t> candidate_order =
        rank_mixed_starts_in_selected_basin(
            endpoints,
            q_values,
            valid_indices,
            core_fraction_policy
        );
    if (candidate_order.empty()) {
        result.failure_reason = FitFailureReason::DegenerateCoreFraction;
        return result;
    }

    result.consensus_size = candidate_order.size();
    result.selected_core_start = candidate_order.front();

    // The physical basin is fixed once, from the original deterministic start
    // grid. Its members are then tried in increasing terminal q. Each member
    // supplies only fresh coordinates for a clean post-selection MIGRAD and
    // MINOS pass. A numerical failure therefore falls back to another endpoint
    // of the same already-selected basin; it never changes the physical basin.
    const auto evaluate_candidate = [&] (
        std::size_t candidate_start
    ) -> MixedFitResult {
        MixedFitResult candidate = result;
        candidate.selected_core_start = candidate_start;
        candidate.selected_migrad = unattempted_migrad();
        candidate.amplitude_profile_migrad = unattempted_migrad();
        candidate.minos_core_radius = unattempted_minos();
        candidate.minos_tail_radius = unattempted_minos();
        candidate.minos_core_fraction = unattempted_minos();
        candidate.minos_amplitude = unattempted_minos();
        candidate.q_min = std::nullopt;
        candidate.core_radius = std::nullopt;
        candidate.tail_radius = std::nullopt;
        candidate.core_fraction = std::nullopt;
        candidate.amplitude = std::nullopt;
        candidate.fitted_pdf.clear();
        candidate.fully_valid = false;
        candidate.failure_reason = FitFailureReason::MigradInvalid;

        const FunctionMinimum& basin_minimum =
            outcomes[candidate_start].minimum;
        const double basin_log_core = basin_minimum.UserState().Value(0U);
        const double basin_log_tail = basin_minimum.UserState().Value(1U);
        const double basin_core_fraction = basin_minimum.UserState().Value(2U);
        const double basin_core_radius = std::exp(basin_log_core);
        const double basin_tail_radius = std::exp(basin_log_tail);
        if (!std::isfinite(basin_log_core) ||
            !std::isfinite(basin_log_tail) ||
            !std::isfinite(basin_core_fraction) ||
            !std::isfinite(basin_core_radius) ||
            !std::isfinite(basin_tail_radius) ||
            basin_core_radius <= 0.0 || basin_tail_radius <= 0.0) {
            candidate.failure_reason = FitFailureReason::NonFiniteMinimum;
            return candidate;
        }
        const FitFailureReason basin_fraction_reason =
            mixed_core_fraction_failure(basin_core_fraction);
        if (basin_fraction_reason != FitFailureReason::None) {
            candidate.failure_reason = basin_fraction_reason;
            return candidate;
        }

        MixedObjective selected_objective(
            family, bins, offset, binning, region
        );
        auto polished = run_mixed_migrad(
            selected_objective,
            basin_core_radius,
            basin_tail_radius,
            basin_core_fraction
        );
        candidate.selected_migrad = polished.second;
        candidate.q_min = candidate.selected_migrad.q_min;
        const FitFailureReason polishing_reason =
            fit_failure_from_migrad(candidate.selected_migrad);
        if (polishing_reason != FitFailureReason::None) {
            candidate.failure_reason = polishing_reason;
            return candidate;
        }

        FunctionMinimum& selected = polished.first;
        const double log_core = selected.UserState().Value(0U);
        const double log_tail = selected.UserState().Value(1U);
        const double core_fraction = selected.UserState().Value(2U);
        const double core_radius = std::exp(log_core);
        const double tail_radius = std::exp(log_tail);
        if (!std::isfinite(log_core) || !std::isfinite(log_tail) ||
            !std::isfinite(core_radius) || !std::isfinite(tail_radius) ||
            !std::isfinite(core_fraction) || core_radius <= 0.0 ||
            tail_radius <= 0.0 || !candidate.q_min.has_value()) {
            candidate.failure_reason = FitFailureReason::NonFiniteMinimum;
            return candidate;
        }
        const FitFailureReason fraction_reason =
            mixed_core_fraction_failure(core_fraction);
        if (fraction_reason != FitFailureReason::None) {
            candidate.failure_reason = fraction_reason;
            return candidate;
        }

        std::vector<double> bin_integrals;
        double amplitude = 0.0;
        try {
            bin_integrals = mixed_bin_integrals(
                family,
                binning,
                region,
                core_radius,
                tail_radius,
                core_fraction
            );
            amplitude = neyman_optimal_mixed_amplitude(
                bins, offset, region, bin_integrals
            );
        } catch (const std::exception&) {
            candidate.failure_reason = FitFailureReason::ObjectiveEvaluation;
            return candidate;
        }

        MnMinos minos(selected_objective, selected);
        const MinosError core_error = minos.Minos(0U);
        candidate.minos_core_radius = minos_diagnostic(core_error);
        FitFailureReason minos_reason = fit_failure_from_minos(
            candidate.minos_core_radius,
            false
        );
        if (minos_reason != FitFailureReason::None) {
            candidate.failure_reason = minos_reason;
            return candidate;
        }

        const MinosError tail_error = minos.Minos(1U);
        candidate.minos_tail_radius = minos_diagnostic(tail_error);
        minos_reason = fit_failure_from_minos(
            candidate.minos_tail_radius,
            false
        );
        if (minos_reason != FitFailureReason::None) {
            candidate.failure_reason = minos_reason;
            return candidate;
        }

        const MinosError fraction_error = minos.Minos(2U);
        candidate.minos_core_fraction = minos_diagnostic(fraction_error);
        minos_reason = fit_failure_from_minos(
            candidate.minos_core_fraction,
            true
        );
        if (minos_reason != FitFailureReason::None) {
            candidate.failure_reason = minos_reason;
            return candidate;
        }
        if (core_fraction + fraction_error.Lower() <= 0.0) {
            candidate.failure_reason = FitFailureReason::MinosLowerLimit;
            return candidate;
        }
        if (core_fraction + fraction_error.Upper() >= 1.0) {
            candidate.failure_reason = FitFailureReason::MinosUpperLimit;
            return candidate;
        }

        MixedAmplitudeProfileObjective amplitude_profile_objective(
            family, bins, offset, binning, region
        );
        auto amplitude_profile = run_mixed_amplitude_profile_migrad(
            amplitude_profile_objective,
            core_radius,
            tail_radius,
            core_fraction,
            amplitude
        );
        candidate.amplitude_profile_migrad = amplitude_profile.second;
        const FitFailureReason amplitude_profile_migrad_reason =
            fit_failure_from_migrad(candidate.amplitude_profile_migrad);
        if (amplitude_profile_migrad_reason != FitFailureReason::None) {
            candidate.failure_reason = amplitude_profile_migrad_reason;
            return candidate;
        }

        FunctionMinimum& amplitude_profile_minimum = amplitude_profile.first;
        const MixedBasinPoint profile_endpoint{
            amplitude_profile_minimum.UserState().Value(0U),
            amplitude_profile_minimum.UserState().Value(1U),
            amplitude_profile_minimum.UserState().Value(2U)
        };
        const MixedBasinPoint selected_endpoint{
            log_core,
            log_tail,
            core_fraction
        };
        if (!same_mixed_basin(profile_endpoint, selected_endpoint)) {
            candidate.failure_reason = FitFailureReason::MigradInvalid;
            return candidate;
        }
        MnMinos amplitude_minos(
            amplitude_profile_objective,
            amplitude_profile_minimum
        );
        const MinosError amplitude_error = amplitude_minos.Minos(3U);
        candidate.minos_amplitude = minos_diagnostic(amplitude_error);
        minos_reason = fit_failure_from_minos(
            candidate.minos_amplitude,
            false
        );
        if (minos_reason != FitFailureReason::None) {
            candidate.failure_reason = minos_reason;
            return candidate;
        }

        try {
            candidate.core_radius = physical_radius_estimate(
                log_core,
                core_error
            );
            candidate.tail_radius = physical_radius_estimate(
                log_tail,
                tail_error
            );
            candidate.core_fraction = direct_parameter_estimate(
                core_fraction,
                fraction_error
            );
            const double profile_log_amplitude =
                amplitude_profile_minimum.UserState().Value(3U);
            const double lower_amplitude = std::exp(
                profile_log_amplitude + amplitude_error.Lower()
            );
            const double upper_amplitude = std::exp(
                profile_log_amplitude + amplitude_error.Upper()
            );
            if (!std::isfinite(amplitude) || amplitude <= 0.0 ||
                !std::isfinite(lower_amplitude) || lower_amplitude <= 0.0 ||
                !std::isfinite(upper_amplitude) || upper_amplitude <= 0.0 ||
                lower_amplitude > amplitude || upper_amplitude < amplitude) {
                throw std::invalid_argument(
                    "HBT mixed fit: invalid amplitude profile interval"
                );
            }
            candidate.amplitude = FitParameterEstimate{
                amplitude,
                amplitude - lower_amplitude,
                upper_amplitude - amplitude
            };
            candidate.fitted_pdf = mixed_integrals_to_pdf(
                bin_integrals,
                amplitude,
                binning
            );
        } catch (const std::exception&) {
            candidate.core_radius = std::nullopt;
            candidate.tail_radius = std::nullopt;
            candidate.core_fraction = std::nullopt;
            candidate.amplitude = std::nullopt;
            candidate.fitted_pdf.clear();
            candidate.failure_reason = FitFailureReason::NonFiniteMinimum;
            return candidate;
        }

        candidate.fully_valid = true;
        candidate.failure_reason = FitFailureReason::None;
        return candidate;
    };

    // Preserve the historical lowest-q endpoint failure as the primary
    // diagnostic if every endpoint in the winning basin fails. If a later
    // endpoint succeeds, its start index and polished diagnostics are the ones
    // published in the result.
    std::optional<MixedFitResult> primary_failure;
    for (const std::size_t candidate_start : candidate_order) {
        MixedFitResult candidate = evaluate_candidate(candidate_start);
        if (!primary_failure.has_value()) {
            primary_failure = candidate;
        }
        if (candidate.fully_valid) {
            return candidate;
        }
    }
    return primary_failure.value();
}

}  // namespace hbt
