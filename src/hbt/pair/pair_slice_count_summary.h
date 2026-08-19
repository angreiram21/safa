/**
 * @file pair_slice_count_summary.h
 * @brief Pair counts grouped by kinetic slice, origin route, and channel.
 */

#ifndef SMASH_AFTERBURNER_ANALYSIS_HBT_PAIR_PAIR_SLICE_COUNT_SUMMARY_H
#define SMASH_AFTERBURNER_ANALYSIS_HBT_PAIR_PAIR_SLICE_COUNT_SUMMARY_H

#include "hbt/config/hbt_config.h"
#include "hbt/pair/pair_count_summary.h"

#include <cstddef>
#include <optional>
#include <vector>

namespace hbt {

/**
 * @brief Origin-routed pair counts stored inside one kinetic slice.
 *
 * The three summaries are nested routes, not disjoint categories. Under
 * OriginMode::All, one physical pair may increment more than one route.
 */
struct PairSliceOriginCounts {
    PairCountSummary routed_P;    ///< Primordial-route counts by channel.
    PairCountSummary routed_PR;   ///< Primordial+rescattering counts.
    PairCountSummary routed_PRD;  ///< Widest-origin counts by channel.
};

/**
 * @brief Counts associated with one configured kinetic slice or cell.
 *
 * A missing index denotes a disabled axis. At least one index is present for
 * every stored entry; no fictitious inclusive slice is created when both axes
 * are disabled.
 */
struct PairSliceCountEntry {
    std::optional<std::size_t> kt_slice_index;  ///< kT index when enabled.
    std::optional<std::size_t> mt_slice_index;  ///< mT index when enabled.
    PairSliceOriginCounts origin_counts;  ///< Counts by origin and channel.
};

/**
 * @brief Complete zero-or-more-slice pair-count layout for one routing policy.
 *
 * pair_slicing preserves the validated axis definitions needed to interpret
 * the stored indices. origin_mode records which nested origin routes were
 * requested. entries are deterministic: kT-only uses increasing kT index,
 * mT-only uses increasing mT index, and kT x mT uses kT-major order.
 */
struct PairSliceCountSummary {
    PairSlicingConfig pair_slicing;  ///< Validated slicing configuration.
    OriginMode origin_mode;          ///< Requested origin-routing policy.
    std::vector<PairSliceCountEntry> entries;  ///< Ordered slice counts.
};

}  // namespace hbt

#endif  // SMASH_AFTERBURNER_ANALYSIS_HBT_PAIR_PAIR_SLICE_COUNT_SUMMARY_H
