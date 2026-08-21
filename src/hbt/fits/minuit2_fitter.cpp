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
    return {
        false,
        reason,
        estimator,
        starts,
        0U,
        0U,
        0U,
        std::nullopt,
        empty,
        unattempted_minos(),
        unattempted_minos(),
        unattempted_minos(),
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
 * @brief Minuit2 objective for a one-radius normalized Gaussian model.
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
     * @brief Evaluate the configured objective from one log-radius parameter.
     * @param parameters Minuit external parameters; index zero is log(R).
     * @return Finite deviance, or the largest finite double for an invalid
     *         model evaluation. Returned minima are validated explicitly.
     */
    double operator()(const std::vector<double>& parameters) const override {
        if (parameters.size() != 1U) {
            return std::numeric_limits<double>::max();
        }
        const double radius = std::exp(parameters[0]);
        try {
            const std::vector<double> probabilities =
                gaussian_bin_probabilities(
                    family_,
                    binning_,
                    region_,
                    radius
                );
            switch (estimator_) {
                case FitEstimator::Poisson:
                    return binned_poisson_deviance(
                        bins_, offset_, region_, probabilities
                    );
                case FitEstimator::Neyman:
                    return binned_neyman_chi_square(
                        bins_, offset_, region_, probabilities
                    );
                case FitEstimator::Pearson:
                    return binned_pearson_chi_square(
                        bins_, offset_, region_, probabilities
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
 * @brief Minuit2 objective for one normalized mixed-model estimator.
 *
 * Parameter order is log(R_core), log(R_tail), f_core. The model probability
 * calculation is identical for every estimator; only the scalar objective
 * applied to the same raw counts and expected counts changes. Borrowed
 * references are retained only during one synchronous fit invocation.
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
     * @param estimator Statistical objective evaluated for this fit only.
     */
    MixedObjective(
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
            const std::vector<double> probabilities = mixed_bin_probabilities(
                family_,
                binning_,
                region_,
                core_radius,
                tail_radius,
                parameters[2]
            );
            switch (estimator_) {
                case FitEstimator::Poisson:
                    return binned_poisson_deviance(
                        bins_, offset_, region_, probabilities
                    );
                case FitEstimator::Neyman:
                    return binned_neyman_chi_square(
                        bins_, offset_, region_, probabilities
                    );
                case FitEstimator::Pearson:
                    return binned_pearson_chi_square(
                        bins_, offset_, region_, probabilities
                    );
            }
        } catch (const std::exception&) {
            return std::numeric_limits<double>::max();
        }
        return std::numeric_limits<double>::max();
    }

    /**
     * @brief Return the MINOS error definition for every supported objective.
     * @return Exactly 1.0, corresponding to Delta(-2 log L)=1 for Poisson and
     *         Delta(chi-square)=1 for the two chi-square estimators.
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
            minimum.UserState().Value(0U)
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
 * @return FunctionMinimum plus its stable diagnostic.
 */
std::pair<FunctionMinimum, MigradDiagnostic> run_gaussian_migrad(
    GaussianObjective& objective,
    double radius_seed
) {
    MnUserParameters parameters;
    parameters.Add("log_r", std::log(radius_seed), 0.1);
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
 * @param core_fraction_seed Initial Gaussian mixture fraction in (0,1).
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
 * @brief One completed deterministic mixed start retained until selection.
 */
struct MixedStartOutcome {
    FunctionMinimum minimum;      ///< Complete Minuit minimum for later MINOS.
    MigradDiagnostic diagnostic;  ///< Stable serialized start diagnostic.
    MixedBasinPoint endpoint;     ///< External endpoint for basin grouping.
};

/**
 * @brief Count valid starts in the connected basin containing one start.
 * @param outcomes Endpoints for every attempted start.
 * @param valid_indices Indices with valid evaluable MIGRAD minima.
 * @param selected_index Valid selected-minimum index.
 * @return Connected same-basin component size containing @p selected_index.
 */
std::size_t selected_mixed_basin_size(
    const std::vector<MixedStartOutcome>& outcomes,
    const std::vector<std::size_t>& valid_indices,
    std::size_t selected_index
) {
    std::vector<bool> visited(outcomes.size(), false);
    std::vector<std::size_t> stack{selected_index};
    visited[selected_index] = true;
    std::size_t count = 0U;
    while (!stack.empty()) {
        const std::size_t current = stack.back();
        stack.pop_back();
        ++count;
        for (const std::size_t candidate : valid_indices) {
            if (!visited[candidate] && same_mixed_basin(
                    outcomes[current].endpoint,
                    outcomes[candidate].endpoint
                )) {
                visited[candidate] = true;
                stack.push_back(candidate);
            }
        }
    }
    return count;
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
    if (selected_bins < 2U || region.selected_count == 0U) {
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
        auto start = run_gaussian_migrad(objective, radius_seeds[index]);
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
    const double radius = std::exp(log_radius);
    if (!std::isfinite(log_radius) || !std::isfinite(radius) ||
        radius <= 0.0 || !result.q_min.has_value()) {
        result.failure_reason = FitFailureReason::NonFiniteMinimum;
        return result;
    }

    std::vector<double> probabilities;
    try {
        probabilities = gaussian_bin_probabilities(
            family,
            binning,
            region,
            radius
        );
    } catch (const std::exception&) {
        result.failure_reason = FitFailureReason::ObjectiveEvaluation;
        return result;
    }

    GaussianObjective selected_objective(
        family, bins, offset, binning, region, estimator
    );
    MnMinos minos(selected_objective, selected);
    const MinosError radius_error = minos.Minos(0U);
    result.minos_radius = minos_diagnostic(radius_error);
    const FitFailureReason minos_reason = fit_failure_from_minos(
        result.minos_radius,
        false
    );
    if (minos_reason != FitFailureReason::None) {
        result.failure_reason = minos_reason;
        return result;
    }

    try {
        result.radius = physical_radius_estimate(log_radius, radius_error);
        result.fitted_pdf = probabilities_to_pdf(probabilities, binning);
    } catch (const std::exception&) {
        result.radius = std::nullopt;
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
    double half_maximum_seed
) {
    const std::size_t selected_bins =
        region.last_bin - region.first_bin + 1U;
    if (selected_bins < 4U || region.selected_count == 0U) {
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
                    family, bins, offset, binning, region, estimator
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

    std::size_t selected_index = valid_indices.front();
    for (const std::size_t valid_index : valid_indices) {
        if (outcomes[valid_index].diagnostic.q_min.value() <
            outcomes[selected_index].diagnostic.q_min.value()) {
            selected_index = valid_index;
        }
    }
    result.selected_core_start = selected_index;
    result.selected_migrad = outcomes[selected_index].diagnostic;
    result.q_min = result.selected_migrad.q_min;
    result.consensus_size = selected_mixed_basin_size(
        outcomes,
        valid_indices,
        selected_index
    );

    FunctionMinimum& selected = outcomes[selected_index].minimum;
    const double log_core = selected.UserState().Value(0U);
    const double log_tail = selected.UserState().Value(1U);
    const double core_fraction = selected.UserState().Value(2U);
    const double core_radius = std::exp(log_core);
    const double tail_radius = std::exp(log_tail);
    if (!std::isfinite(core_radius) || !std::isfinite(tail_radius) ||
        !std::isfinite(core_fraction) || core_radius <= 0.0 ||
        tail_radius <= 0.0 || !result.q_min.has_value()) {
        result.failure_reason = FitFailureReason::NonFiniteMinimum;
        return result;
    }
    const FitFailureReason fraction_reason =
        mixed_core_fraction_failure(core_fraction);
    if (fraction_reason != FitFailureReason::None) {
        result.failure_reason = fraction_reason;
        return result;
    }

    std::vector<double> probabilities;
    try {
        probabilities = mixed_bin_probabilities(
            family,
            binning,
            region,
            core_radius,
            tail_radius,
            core_fraction
        );
    } catch (const std::exception&) {
        result.failure_reason = FitFailureReason::ObjectiveEvaluation;
        return result;
    }

    MixedObjective selected_objective(
        family, bins, offset, binning, region, estimator
    );
    MnMinos minos(selected_objective, selected);
    const MinosError core_error = minos.Minos(0U);
    result.minos_core_radius = minos_diagnostic(core_error);
    FitFailureReason minos_reason = fit_failure_from_minos(
        result.minos_core_radius,
        false
    );
    if (minos_reason != FitFailureReason::None) {
        result.failure_reason = minos_reason;
        return result;
    }

    const MinosError tail_error = minos.Minos(1U);
    result.minos_tail_radius = minos_diagnostic(tail_error);
    minos_reason = fit_failure_from_minos(result.minos_tail_radius, false);
    if (minos_reason != FitFailureReason::None) {
        result.failure_reason = minos_reason;
        return result;
    }

    const MinosError fraction_error = minos.Minos(2U);
    result.minos_core_fraction = minos_diagnostic(fraction_error);
    minos_reason = fit_failure_from_minos(result.minos_core_fraction, true);
    if (minos_reason != FitFailureReason::None) {
        result.failure_reason = minos_reason;
        return result;
    }
    if (core_fraction + fraction_error.Lower() <= 0.0) {
        result.failure_reason = FitFailureReason::MinosLowerLimit;
        return result;
    }
    if (core_fraction + fraction_error.Upper() >= 1.0) {
        result.failure_reason = FitFailureReason::MinosUpperLimit;
        return result;
    }

    try {
        result.core_radius = physical_radius_estimate(log_core, core_error);
        result.tail_radius = physical_radius_estimate(log_tail, tail_error);
        result.core_fraction = direct_parameter_estimate(
            core_fraction,
            fraction_error
        );
        result.fitted_pdf = probabilities_to_pdf(probabilities, binning);
    } catch (const std::exception&) {
        result.core_radius = std::nullopt;
        result.tail_radius = std::nullopt;
        result.core_fraction = std::nullopt;
        result.fitted_pdf.clear();
        result.failure_reason = FitFailureReason::NonFiniteMinimum;
        return result;
    }

    result.fully_valid = true;
    result.failure_reason = FitFailureReason::None;
    return result;
}

}  // namespace hbt
