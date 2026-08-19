/**
 * @file pair_slice_routing.h
 * @brief Routing of validated pair kinematics to configured kT/mT slices.
 */

#ifndef SMASH_AFTERBURNER_ANALYSIS_HBT_PAIR_PAIR_SLICE_ROUTING_H
#define SMASH_AFTERBURNER_ANALYSIS_HBT_PAIR_PAIR_SLICE_ROUTING_H

#include "hbt/config/hbt_config.h"
#include "hbt/pair/pair_kinematics.h"

#include <cstddef>
#include <optional>

namespace hbt {

/**
 * @brief Unique configured slice reached by one validated physical pair.
 *
 * An index is present only for an enabled axis. When both axes are enabled,
 * the two indices identify one Cartesian kT x mT cell. This structure stores
 * routing only; it does not own origin membership or histogram state.
 */
struct PairSliceRoute {
    std::optional<std::size_t> kt_slice_index;  ///< Routed kT slice index.
    std::optional<std::size_t> mt_slice_index;  ///< Routed mT slice index.
    /// Flat slice index in the configured deterministic slice layout.
    std::size_t flat_slice_index;
};

/**
 * @brief Route one validated pair to the configured kinetic slice.
 * @param kinematics Already calculated and numerically valid pair kinematics.
 * @param slicing Already validated kT/mT slicing configuration.
 * @return The unique slice route, or std::nullopt when slicing is disabled or
 *         the pair lies outside at least one enabled slicing axis.
 * @throws std::invalid_argument If a required kinematic value is non-finite
 *         or an enabled axis has fewer than two edges.
 *
 * Enabled bins use half-open intervals [edge_i, edge_(i+1)). Disabled axes are
 * ignored completely, including any retained validated bin edges. When both
 * axes are enabled, a pair must lie inside both axes to obtain a route. The
 * flat slice index is derived once here and reused by downstream consumers.
 *
 * Full edge validation is an upstream configuration responsibility and is not
 * repeated for each pair. Routing uses binary search over validated edges.
 *
 * This function performs no pair formation, origin routing, pair-kinematics
 * calculation, histogramming, or output.
 */
[[nodiscard]] std::optional<PairSliceRoute> route_pair_to_slices(
    const PairKinematics& kinematics,
    const PairSlicingConfig& slicing
);

}  // namespace hbt

#endif  // SMASH_AFTERBURNER_ANALYSIS_HBT_PAIR_PAIR_SLICE_ROUTING_H
