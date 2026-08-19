/**
 * @file pair_count_accumulator.h
 * @brief Construction and accumulation of primitive-channel pair counts.
 */

#ifndef SMASH_AFTERBURNER_ANALYSIS_HBT_PAIR_PAIR_COUNT_ACCUMULATOR_H
#define SMASH_AFTERBURNER_ANALYSIS_HBT_PAIR_PAIR_COUNT_ACCUMULATOR_H

#include "hbt/channels/primitive_channel.h"
#include "hbt/pair/pair_count_summary.h"

#include <vector>

namespace hbt {

/**
 * @brief Build zeroed pair counts for an ordered required-channel list.
 * @param channels Unique valid primitive channels in required output order.
 * @return One zero-count entry for every channel, preserving input order.
 * @throws std::invalid_argument If a channel identifier is invalid or repeated.
 */
PairCountSummary make_zero_pair_count_summary(
    const std::vector<PrimitiveChannelId>& channels
);

/**
 * @brief Add one local pair-count summary to an existing total.
 * @param total Existing accumulated counts to update.
 * @param subevent Local subevent counts to add.
 * @throws std::invalid_argument If either summary is structurally invalid or
 *         their channel lists differ in size, identifiers, or order.
 * @throws std::overflow_error If any uint64_t addition would overflow.
 *
 * The function validates the complete operation before mutating @p total. A
 * failed accumulation therefore leaves @p total unchanged.
 */
void accumulate_pair_counts(
    PairCountSummary& total,
    const PairCountSummary& subevent
);

}  // namespace hbt

#endif  // SMASH_AFTERBURNER_ANALYSIS_HBT_PAIR_PAIR_COUNT_ACCUMULATOR_H
