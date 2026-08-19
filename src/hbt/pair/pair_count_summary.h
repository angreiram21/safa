/**
 * @file pair_count_summary.h
 * @brief Pair-count data produced by primitive-channel iteration.
 */

#ifndef SMASH_AFTERBURNER_ANALYSIS_HBT_PAIR_PAIR_COUNT_SUMMARY_H
#define SMASH_AFTERBURNER_ANALYSIS_HBT_PAIR_PAIR_COUNT_SUMMARY_H

#include "hbt/channels/primitive_channel.h"

#include <cstdint>
#include <vector>

namespace hbt {

/**
 * @brief Number of physical pairs iterated for one primitive HBT channel.
 */
struct PairChannelCount {
    PrimitiveChannelId channel;  ///< Primitive channel identifier.
    std::uint64_t pair_count;    ///< Number of pairs delivered for the channel.
};

/**
 * @brief Ordered per-channel pair counts produced by one iteration.
 *
 * Entries preserve the order of the required primitive channels supplied to
 * the pair iterator. A valid requested channel has one entry even when its
 * physical pair count is zero.
 */
struct PairCountSummary {
    std::vector<PairChannelCount> channels;  ///< Ordered channel counts.
};

}  // namespace hbt

#endif  // SMASH_AFTERBURNER_ANALYSIS_HBT_PAIR_PAIR_COUNT_SUMMARY_H
