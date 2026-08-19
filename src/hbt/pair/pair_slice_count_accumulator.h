/**
 * @file pair_slice_count_accumulator.h
 * @brief Construction and accumulation of pair counts by slice and origin.
 */

#ifndef SMASH_AFTERBURNER_ANALYSIS_HBT_PAIR_PAIR_SLICE_COUNT_ACCUMULATOR_H
#define SMASH_AFTERBURNER_ANALYSIS_HBT_PAIR_PAIR_SLICE_COUNT_ACCUMULATOR_H

#include "hbt/channels/primitive_channel.h"
#include "hbt/config/hbt_config.h"
#include "hbt/pair/pair_origin_routing.h"
#include "hbt/pair/pair_slice_count_summary.h"
#include "hbt/pair/pair_slice_routing.h"

#include <cstddef>
#include <vector>

namespace hbt {

/**
 * @brief Build deterministic zeroed counts for all configured kinetic slices.
 * @param slicing Validated pair-slicing configuration.
 * @param origin_mode Requested nested origin-routing policy.
 * @param channels Unique valid primitive channels in output order.
 * @return Zeroed slice/origin/channel counts. Both-disabled slicing produces
 *         an empty entries vector rather than an inclusive dummy slice.
 * @throws std::invalid_argument If origin_mode or the channel list is invalid.
 */
PairSliceCountSummary make_zero_pair_slice_count_summary(
    const PairSlicingConfig& slicing,
    OriginMode origin_mode,
    const std::vector<PrimitiveChannelId>& channels
);

/**
 * @brief Count one already routed pair without recalculating its slice.
 * @param summary Slice-count summary to update.
 * @param slice_route Unique kinetic route calculated once for the pair.
 * @param origin_routes Requested compatible origin routes for the same pair.
 * @param channel_index Ordered channel index already known by the caller.
 * @param channel Primitive channel expected at @p channel_index.
 * @throws std::invalid_argument If routing state or summary structure is
 *         inconsistent with the summary configuration and origin mode.
 * @throws std::overflow_error If any selected uint64_t count would overflow.
 *
 * The function performs only counter updates. It does not inspect kT/mT,
 * search bin edges, calculate origin memberships, or form a physical pair.
 */
void increment_pair_slice_count(
    PairSliceCountSummary& summary,
    const PairSliceRoute& slice_route,
    const PairOriginRoutes& origin_routes,
    std::size_t channel_index,
    PrimitiveChannelId channel
);

/**
 * @brief Add one structurally identical slice-count summary into a total.
 * @param total Existing accumulated counts to update.
 * @param subevent Local counts to add.
 * @throws std::invalid_argument If slicing, origin mode, slice layout, or
 *         channel identity/order differs between the two summaries.
 * @throws std::overflow_error If any uint64_t addition would overflow.
 *
 * The complete operation is validated before mutation, so failure leaves
 * @p total unchanged.
 */
void accumulate_pair_slice_counts(
    PairSliceCountSummary& total,
    const PairSliceCountSummary& subevent
);

}  // namespace hbt

#endif  // SMASH_AFTERBURNER_ANALYSIS_HBT_PAIR_PAIR_SLICE_COUNT_ACCUMULATOR_H
