/**
 * @file fit_results.cpp
 * @brief Stable diagnostic token mappings for post-sample result types.
 */

#include "hbt/fits/fit_results.h"

#include <algorithm>
#include <cmath>
#include <limits>
#include <stdexcept>

namespace hbt {
namespace {

/**
 * @brief Partition valid mixed endpoints into connected numerical basins.
 * @param endpoints Final endpoint for every attempted start.
 * @param valid_indices Indices included in the numerical-basin graph.
 * @return Connected components in deterministic valid-index traversal order.
 * @throws std::out_of_range If a valid start index has no matching endpoint.
 */
std::vector<std::vector<std::size_t>> mixed_basin_components(
    const std::vector<MixedBasinPoint>& endpoints,
    const std::vector<std::size_t>& valid_indices
) {
    std::vector<bool> visited(endpoints.size(), false);
    std::vector<std::vector<std::size_t>> components;
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
        components.push_back(std::move(component));
    }
    return components;
}

}  // namespace

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
    const auto components = mixed_basin_components(endpoints, valid_indices);
    std::vector<std::size_t> best;
    for (const auto& component : components) {
        if (component.size() > best.size()) {
            best = component;
        }
    }
    return best;
}

std::size_t select_mixed_start_by_half_maximum_basin(
    const std::vector<MixedBasinPoint>& endpoints,
    const std::vector<double>& q_values,
    const std::vector<std::size_t>& valid_indices,
    double half_maximum_radius
) {
    if (endpoints.size() != q_values.size()) {
        throw std::invalid_argument(
            "mixed basin selection: endpoint and q arrays differ in size"
        );
    }
    if (!std::isfinite(half_maximum_radius) || half_maximum_radius <= 0.0) {
        throw std::invalid_argument(
            "mixed basin selection: R_HM must be finite and positive"
        );
    }
    if (valid_indices.empty()) {
        throw std::invalid_argument(
            "mixed basin selection: at least one valid start is required"
        );
    }

    for (const std::size_t index : valid_indices) {
        if (index >= endpoints.size()) {
            throw std::out_of_range(
                "mixed basin selection: valid start index is out of range"
            );
        }
        if (!std::isfinite(q_values[index])) {
            throw std::invalid_argument(
                "mixed basin selection: valid start q must be finite"
            );
        }
    }

    const auto components = mixed_basin_components(endpoints, valid_indices);
    const double log_half_maximum = std::log(half_maximum_radius);
    std::size_t selected_component = 0U;
    double selected_distance = std::numeric_limits<double>::infinity();
    double selected_component_q = std::numeric_limits<double>::infinity();
    std::size_t selected_component_index =
        std::numeric_limits<std::size_t>::max();

    for (std::size_t component_index = 0U;
         component_index < components.size();
         ++component_index) {
        const auto& component = components[component_index];
        double log_core_sum = 0.0;
        double component_q = std::numeric_limits<double>::infinity();
        std::size_t component_lowest_index =
            std::numeric_limits<std::size_t>::max();
        for (const std::size_t index : component) {
            log_core_sum += endpoints[index].log_core_radius;
            component_q = std::min(component_q, q_values[index]);
            component_lowest_index = std::min(component_lowest_index, index);
        }
        const double mean_log_core =
            log_core_sum / static_cast<double>(component.size());
        const double distance = std::fabs(mean_log_core - log_half_maximum);
        if (distance < selected_distance ||
            (distance == selected_distance &&
             (component_q < selected_component_q ||
              (component_q == selected_component_q &&
               component_lowest_index < selected_component_index)))) {
            selected_component = component_index;
            selected_distance = distance;
            selected_component_q = component_q;
            selected_component_index = component_lowest_index;
        }
    }

    const auto& basin = components[selected_component];
    std::size_t selected_start = basin.front();
    for (const std::size_t index : basin) {
        if (q_values[index] < q_values[selected_start] ||
            (q_values[index] == q_values[selected_start] &&
             index < selected_start)) {
            selected_start = index;
        }
    }
    return selected_start;
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
        case FitFailureReason::InvalidHalfMaximumSeed:
            return "invalid_half_maximum_seed";
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

const char* fit_estimator_token(FitEstimator estimator) {
    switch (estimator) {
        case FitEstimator::Poisson:
            return "poisson";
        case FitEstimator::Neyman:
            return "neyman";
        case FitEstimator::Pearson:
            return "pearson";
    }

    throw std::invalid_argument("invalid FitEstimator");
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
