/**
 * @file fit_results.cpp
 * @brief Stable diagnostic token mappings for post-sample result types.
 */

#include "hbt/fits/fit_results.h"

#include <cmath>
#include <stdexcept>

namespace hbt {

bool same_mixed_basin(
    const MixedBasinPoint& lhs,
    const MixedBasinPoint& rhs
) {
    return std::fabs(lhs.log_core_radius - rhs.log_core_radius) <=
            kMixedBasinLogRadiusTolerance &&
        std::fabs(lhs.log_tail_radius - rhs.log_tail_radius) <=
            kMixedBasinLogRadiusTolerance &&
        std::fabs(lhs.core_fraction - rhs.core_fraction) <=
            kMixedBasinCoreFractionTolerance;
}

std::vector<std::size_t> largest_mixed_basin_group(
    const std::vector<MixedBasinPoint>& endpoints,
    const std::vector<std::size_t>& valid_indices
) {
    std::vector<bool> visited(endpoints.size(), false);
    std::vector<std::size_t> best;
    for (const std::size_t seed : valid_indices) {
        if (seed >= endpoints.size()) {
            throw std::out_of_range(
                "mixed basin grouping: valid start index is out of range"
            );
        }
        if (visited[seed]) {
            continue;
        }
        std::vector<std::size_t> component;
        std::vector<std::size_t> stack{seed};
        visited[seed] = true;
        while (!stack.empty()) {
            const std::size_t current = stack.back();
            stack.pop_back();
            component.push_back(current);
            for (const std::size_t candidate : valid_indices) {
                if (candidate >= endpoints.size()) {
                    throw std::out_of_range(
                        "mixed basin grouping: valid start index is out of range"
                    );
                }
                if (!visited[candidate] && same_mixed_basin(
                        endpoints[current],
                        endpoints[candidate]
                    )) {
                    visited[candidate] = true;
                    stack.push_back(candidate);
                }
            }
        }
        if (component.size() > best.size()) {
            best = std::move(component);
        }
    }
    return best;
}

FitFailureReason fit_failure_from_migrad(
    const MigradDiagnostic& diagnostic
) {
    if (diagnostic.objective_failure) {
        return FitFailureReason::ObjectiveEvaluation;
    }
    if (diagnostic.reached_call_limit) {
        return FitFailureReason::MigradCallLimit;
    }
    if (diagnostic.above_max_edm) {
        return FitFailureReason::MigradAboveMaxEdm;
    }
    if (!diagnostic.q_min.has_value()) {
        return FitFailureReason::NonFiniteMinimum;
    }
    if (!diagnostic.valid_covariance) {
        return FitFailureReason::MigradInvalid;
    }
    if (!diagnostic.valid) {
        return FitFailureReason::MigradInvalid;
    }
    return FitFailureReason::None;
}

FitFailureReason fit_failure_from_minos(
    const MinosDiagnostic& diagnostic,
    bool reject_limits
) {
    if (diagnostic.lower_new_minimum) {
        return FitFailureReason::MinosLowerNewMinimum;
    }
    if (diagnostic.upper_new_minimum) {
        return FitFailureReason::MinosUpperNewMinimum;
    }
    if (diagnostic.lower_call_limit) {
        return FitFailureReason::MinosLowerCallLimit;
    }
    if (diagnostic.upper_call_limit) {
        return FitFailureReason::MinosUpperCallLimit;
    }
    if (reject_limits && diagnostic.at_lower_limit) {
        return FitFailureReason::MinosLowerLimit;
    }
    if (reject_limits && diagnostic.at_upper_limit) {
        return FitFailureReason::MinosUpperLimit;
    }
    if (!diagnostic.lower_valid) {
        return FitFailureReason::MinosLowerInvalid;
    }
    if (!diagnostic.upper_valid) {
        return FitFailureReason::MinosUpperInvalid;
    }
    return FitFailureReason::None;
}

FitFailureReason mixed_core_fraction_failure(double core_fraction) {
    if (!std::isfinite(core_fraction) || core_fraction < 0.0 ||
        core_fraction > 1.0) {
        return FitFailureReason::NonFiniteMinimum;
    }
    if (core_fraction == 0.0 || core_fraction == 1.0) {
        return FitFailureReason::DegenerateCoreFraction;
    }
    return FitFailureReason::None;
}

const char* fit_failure_reason_token(FitFailureReason reason) {
    switch (reason) {
        case FitFailureReason::None:
            return "none";
        case FitFailureReason::NotApplicable:
            return "not_applicable";
        case FitFailureReason::EmptyHistogram:
            return "empty_histogram";
        case FitFailureReason::InsufficientBins:
            return "insufficient_bins";
        case FitFailureReason::InsufficientStatistics:
            return "insufficient_statistics";
        case FitFailureReason::InvalidMomentSeed:
            return "invalid_moment_seed";
        case FitFailureReason::InvalidGaussianCoreAnchor:
            return "invalid_gaussian_core_anchor";
        case FitFailureReason::NoBasinConsensus:
            return "no_basin_consensus";
        case FitFailureReason::ObjectiveEvaluation:
            return "objective_evaluation_failure";
        case FitFailureReason::MigradInvalid:
            return "migrad_invalid";
        case FitFailureReason::MigradCallLimit:
            return "migrad_call_limit";
        case FitFailureReason::MigradAboveMaxEdm:
            return "migrad_above_max_edm";
        case FitFailureReason::NonFiniteMinimum:
            return "non_finite_minimum";
        case FitFailureReason::DegenerateCoreFraction:
            return "degenerate_core_fraction";
        case FitFailureReason::MinosLowerInvalid:
            return "minos_lower_invalid";
        case FitFailureReason::MinosUpperInvalid:
            return "minos_upper_invalid";
        case FitFailureReason::MinosLowerCallLimit:
            return "minos_lower_call_limit";
        case FitFailureReason::MinosUpperCallLimit:
            return "minos_upper_call_limit";
        case FitFailureReason::MinosLowerLimit:
            return "minos_lower_limit";
        case FitFailureReason::MinosUpperLimit:
            return "minos_upper_limit";
        case FitFailureReason::MinosLowerNewMinimum:
            return "minos_lower_new_minimum";
        case FitFailureReason::MinosUpperNewMinimum:
            return "minos_upper_new_minimum";
    }

    throw std::invalid_argument("invalid FitFailureReason");
}

const char* delta_t_status_token(DeltaTStatisticsStatus status) {
    switch (status) {
        case DeltaTStatisticsStatus::Valid:
            return "valid";
        case DeltaTStatisticsStatus::EmptyHistogram:
            return "empty_histogram";
        case DeltaTStatisticsStatus::InsufficientCount:
            return "insufficient_count";
        case DeltaTStatisticsStatus::InvalidVariance:
            return "invalid_variance";
    }

    throw std::invalid_argument("invalid DeltaTStatisticsStatus");
}

}  // namespace hbt
