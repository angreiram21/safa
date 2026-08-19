/**
 * @file pair_slice_routing.cpp
 * @brief Pair kT/mT slice-routing implementation.
 */

#include "hbt/pair/pair_slice_routing.h"

#include <algorithm>
#include <cmath>
#include <cstddef>
#include <optional>
#include <stdexcept>
#include <string>
#include <vector>

namespace hbt {
namespace {

/**
 * @brief Require the minimum shape needed to route one enabled axis safely.
 * @param axis Already validated slicing-axis configuration.
 * @param axis_name Name used in structural-error diagnostics.
 * @throws std::invalid_argument If an enabled axis has fewer than two edges.
 *
 * Full scientific edge validation belongs to configuration loading and is not
 * repeated in the per-pair routing hot path.
 */
void require_routable_axis_shape(
    const PairSlicingAxisConfig& axis,
    const char* axis_name
) {
    if (axis.enabled && axis.bin_edges_gev.size() < 2U) {
        throw std::invalid_argument(
            std::string("route_pair_to_slices(): enabled ") +
            axis_name + " axis has fewer than two edges"
        );
    }
}

/**
 * @brief Find the half-open slice containing one finite scalar.
 * @param value Finite kinematic scalar.
 * @param edges Strictly increasing validated bin edges.
 * @return Zero-based slice index, or std::nullopt when outside all slices.
 *
 * The lookup uses binary search and therefore does not scan all configured
 * slices for each physical pair.
 */
std::optional<std::size_t> locate_slice(
    double value,
    const std::vector<double>& edges
) noexcept {
    if (value < edges.front() || value >= edges.back()) {
        return std::nullopt;
    }

    const auto upper = std::upper_bound(
        edges.begin(),
        edges.end(),
        value
    );
    return static_cast<std::size_t>(upper - edges.begin() - 1);
}

}  // namespace

std::optional<PairSliceRoute> route_pair_to_slices(
    const PairKinematics& kinematics,
    const PairSlicingConfig& slicing
) {
    if (!slicing.kt.enabled && !slicing.mt.enabled) {
        return std::nullopt;
    }

    require_routable_axis_shape(slicing.kt, "kT");
    require_routable_axis_shape(slicing.mt, "mT");

    PairSliceRoute route{std::nullopt, std::nullopt, 0U};

    if (slicing.kt.enabled) {
        if (!std::isfinite(kinematics.kt_gev)) {
            throw std::invalid_argument(
                "route_pair_to_slices(): non-finite kT"
            );
        }

        route.kt_slice_index = locate_slice(
            kinematics.kt_gev,
            slicing.kt.bin_edges_gev
        );
        if (!route.kt_slice_index.has_value()) {
            return std::nullopt;
        }
    }

    if (slicing.mt.enabled) {
        if (!std::isfinite(kinematics.mt_gev)) {
            throw std::invalid_argument(
                "route_pair_to_slices(): non-finite mT"
            );
        }

        route.mt_slice_index = locate_slice(
            kinematics.mt_gev,
            slicing.mt.bin_edges_gev
        );
        if (!route.mt_slice_index.has_value()) {
            return std::nullopt;
        }
    }

    if (slicing.kt.enabled && slicing.mt.enabled) {
        const std::size_t mt_slice_count =
            slicing.mt.bin_edges_gev.size() - 1U;
        route.flat_slice_index =
            route.kt_slice_index.value() * mt_slice_count +
            route.mt_slice_index.value();
    } else if (slicing.kt.enabled) {
        route.flat_slice_index = route.kt_slice_index.value();
    } else {
        route.flat_slice_index = route.mt_slice_index.value();
    }

    return route;
}

}  // namespace hbt
