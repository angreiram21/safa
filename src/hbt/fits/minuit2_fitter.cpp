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
    return {false, false, false, false, false, 0, std::nullopt};
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
 * @return Empty numerical result carrying @p reason.
 */
GaussianFitResult invalid_gaussian(FitFailureReason reason) {
    return {
        false,
        reason,
        unattempted_migrad(),
        unattempted_minos(),
        std::nullopt,
        std::nullopt,
        {}
    };
}

/**
 * @brief Build an invalid mixed result before any minimization.
 * @param reason Exact invalidity reason.
 * @return Empty numerical result carrying @p reason.
 */
MixedFitResult invalid_mixed(FitFailureReason reason) {
    const MigradDiagnostic empty = unattempted_migrad();
    return {
        false,
        reason,
        {empty, empty},
        0U,
        std::nullopt,
        false,
        false,
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
     * @brief Construct a synchronous Gaussian likelihood objective.
     * @param family OSL or radial physical family.
     * @param bins Borrowed raw uint64_t count storage.
     * @param offset First logical histogram counter.
     * @param binning Borrowed validated uniform binning.
     * @param region Selected contiguous statistical region, copied by value.
     */
    GaussianObjective(
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
     * @brief Evaluate -2 log L from one log-radius parameter.
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
            return binned_poisson_deviance(
                bins_,
                offset_,
                region_,
                gaussian_bin_probabilities(
                    family_,
                    binning_,
                    region_,
                    radius
                )
            );
        } catch (const std::exception&) {
            return std::numeric_limits<double>::max();
        }
    }

    /**
     * @brief Return Minuit error definition for a -2 log likelihood.
     * @return Exactly 1.0 as required by the post-sample MINOS contract.
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
 * @brief Minuit2 objective for the normalized mixed shape model.
 *
 * Parameter order is log(R_core), log(R_tail), f_core. Borrowed references are
 * retained only during one synchronous fit invocation.
 */
class MixedObjective final : public ROOT::Minuit2::FCNBase {
public:
    /**
     * @brief Construct a synchronous mixed-model likelihood objective.
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
     * @brief Evaluate -2 log L from log radii and core fraction.
     * @param parameters Minuit external parameters in fixed model order.
     * @return Finite deviance, or the largest finite double for an invalid
     *         model evaluation. Returned minima are validated explicitly.
     */
    double operator()(const std::vector<double>& parameters) const override {
        if (parameters.size() != 3U) {
            return std::numeric_limits<double>::max();
        }
        const double core_radius = std::exp(parameters[0]);
        const double tail_radius = std::exp(parameters[1]);
        try {
            return binned_poisson_deviance(
                bins_,
                offset_,
                region_,
                mixed_bin_probabilities(
                    family_,
                    binning_,
                    region_,
                    core_radius,
                    tail_radius,
                    parameters[2]
                )
            );
        } catch (const std::exception&) {
            return std::numeric_limits<double>::max();
        }
    }

    /**
     * @brief Return Minuit error definition for a -2 log likelihood.
     * @return Exactly 1.0 as required by the post-sample MINOS contract.
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
 * @brief Execute one independent mixed-model MIGRAD start.
 * @param objective Borrowed synchronous Minuit objective.
 * @param core_seed Strictly positive Gaussian core radius seed.
 * @param tail_seed Strictly positive exponential tail radius seed.
 * @return FunctionMinimum plus its stable diagnostic.
 *
 * The numerical parameter steps are Minuit controls rather than scientific
 * configuration. Standard strategy, call budget, and tolerance are retained.
 */
std::pair<FunctionMinimum, MigradDiagnostic> run_mixed_migrad(
    MixedObjective& objective,
    double core_seed,
    double tail_seed
) {
    MnUserParameters parameters;
    parameters.Add("log_r_core", std::log(core_seed), 0.1);
    parameters.Add("log_r_tail", std::log(tail_seed), 0.1);
    parameters.Add("f_core", 0.5, 0.05, 0.0, 1.0);
    MnMigrad migrad(objective, parameters);
    FunctionMinimum minimum = migrad();
    const bool objective_failure =
        !mixed_minimum_is_evaluable(objective, minimum);
    return {
        minimum,
        migrad_diagnostic(minimum, objective_failure)
    };
}

}  // namespace

GaussianFitResult fit_gaussian_model(
    FitObservableFamily family,
    const std::vector<std::uint64_t>& bins,
    std::size_t offset,
    const HistogramBinningConfig& binning,
    const StatisticalRegion& region
) {
    const std::size_t selected_bins =
        region.last_bin - region.first_bin + 1U;
    if (selected_bins < 2U || region.selected_count == 0U) {
        return invalid_gaussian(FitFailureReason::InsufficientBins);
    }

    const std::optional<RadiusSeeds> seeds = derive_radius_seeds(
        family,
        bins,
        offset,
        binning,
        region
    );
    if (!seeds.has_value()) {
        return invalid_gaussian(FitFailureReason::InvalidMomentSeed);
    }

    GaussianObjective objective(family, bins, offset, binning, region);
    MnUserParameters parameters;
    parameters.Add("log_r", std::log(seeds->gaussian), 0.1);
    MnMigrad migrad(objective, parameters);
    FunctionMinimum minimum = migrad();

    GaussianFitResult result = invalid_gaussian(
        FitFailureReason::MigradInvalid
    );
    const bool objective_failure =
        !gaussian_minimum_is_evaluable(objective, minimum);
    result.migrad = migrad_diagnostic(minimum, objective_failure);
    result.q_min = result.migrad.q_min;
    const FitFailureReason migrad_reason =
        fit_failure_from_migrad(result.migrad);
    if (migrad_reason != FitFailureReason::None) {
        result.failure_reason = migrad_reason;
        return result;
    }

    const double log_radius = minimum.UserState().Value(0U);
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

    MnMinos minos(objective, minimum);
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
    const GaussianFitResult& gaussian_result
) {
    const std::size_t selected_bins =
        region.last_bin - region.first_bin + 1U;
    if (selected_bins < 4U || region.selected_count == 0U) {
        return invalid_mixed(FitFailureReason::InsufficientBins);
    }

    const std::optional<RadiusSeeds> seeds = derive_radius_seeds(
        family,
        bins,
        offset,
        binning,
        region
    );
    if (!seeds.has_value()) {
        return invalid_mixed(FitFailureReason::InvalidMomentSeed);
    }

    MixedFitResult result = invalid_mixed(
        FitFailureReason::MigradInvalid
    );
    MixedObjective objective_a(family, bins, offset, binning, region);
    auto start_a = run_mixed_migrad(
        objective_a,
        seeds->gaussian,
        seeds->exponential
    );
    result.starts[0U] = start_a.second;
    result.starts_attempted = 1U;

    std::optional<std::pair<FunctionMinimum, MigradDiagnostic>> start_b;
    if (gaussian_result.fully_valid && gaussian_result.radius.has_value()) {
        MixedObjective objective_b(family, bins, offset, binning, region);
        start_b = run_mixed_migrad(
            objective_b,
            gaussian_result.radius->value,
            seeds->exponential
        );
        result.starts[1U] = start_b->second;
        result.starts_attempted = 2U;
    }

    const FitFailureReason start_a_reason =
        fit_failure_from_migrad(start_a.second);
    const FitFailureReason start_b_reason = start_b.has_value()
        ? fit_failure_from_migrad(start_b->second)
        : FitFailureReason::MigradInvalid;
    const bool valid_a = start_a_reason == FitFailureReason::None;
    const bool valid_b = start_b.has_value() &&
        start_b_reason == FitFailureReason::None;
    if (!valid_a && !valid_b) {
        if (start_a_reason == FitFailureReason::ObjectiveEvaluation ||
            start_b_reason == FitFailureReason::ObjectiveEvaluation) {
            result.failure_reason = FitFailureReason::ObjectiveEvaluation;
        } else if (start_b.has_value()) {
            result.failure_reason = start_b_reason;
        } else {
            result.failure_reason = start_a_reason;
        }
        return result;
    }

    const FunctionMinimum* selected = nullptr;
    std::size_t selected_index = 0U;
    if (valid_a && valid_b) {
        const double q_a = start_a.second.q_min.value();
        const double q_b = start_b->second.q_min.value();
        result.exact_q_tie = q_a == q_b;
        if (q_b < q_a) {
            selected = &start_b->first;
            selected_index = 1U;
        } else {
            selected = &start_a.first;
        }
    } else if (valid_a) {
        selected = &start_a.first;
    } else {
        selected = &start_b->first;
        selected_index = 1U;
    }

    result.selected_start = selected_index;
    result.selected_migrad = result.starts[selected_index];
    result.q_min = result.selected_migrad.q_min;

    const double log_core = selected->UserState().Value(0U);
    const double log_tail = selected->UserState().Value(1U);
    const double core_fraction = selected->UserState().Value(2U);
    const double core_radius = std::exp(log_core);
    const double tail_radius = std::exp(log_tail);
    if (!std::isfinite(core_radius) || !std::isfinite(tail_radius) ||
        !std::isfinite(core_fraction) || core_radius <= 0.0 ||
        tail_radius <= 0.0 || !result.q_min.has_value()) {
        result.failure_reason = FitFailureReason::NonFiniteMinimum;
        return result;
    }
    result.tail_below_core = tail_radius < core_radius;
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
        family,
        bins,
        offset,
        binning,
        region
    );
    MnMinos minos(selected_objective, *selected);
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
    minos_reason = fit_failure_from_minos(
        result.minos_core_fraction,
        true
    );
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
        result.core_radius = physical_radius_estimate(
            log_core,
            core_error
        );
        result.tail_radius = physical_radius_estimate(
            log_tail,
            tail_error
        );
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
