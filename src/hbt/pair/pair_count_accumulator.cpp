/**
 * @file pair_count_accumulator.cpp
 * @brief Primitive-channel pair-count accumulation implementation.
 */

#include "hbt/pair/pair_count_accumulator.h"

#include "hbt/channels/channel_catalog.h"

#include <cstddef>
#include <cstdint>
#include <limits>
#include <stdexcept>

namespace hbt {
namespace {

/**
 * @brief Validate one ordered pair-count channel list.
 * @param summary Pair-count summary whose channel structure is checked.
 * @throws std::invalid_argument If a channel identifier is invalid or repeated.
 */
void validate_pair_count_summary_channels(
    const PairCountSummary& summary
) {
    for (std::size_t index = 0U; index < summary.channels.size(); ++index) {
        const PrimitiveChannelId channel = summary.channels[index].channel;
        static_cast<void>(primitive_channel_definition(channel));

        for (std::size_t previous = 0U; previous < index; ++previous) {
            if (summary.channels[previous].channel == channel) {
                throw std::invalid_argument(
                    "pair count summary contains duplicate channel"
                );
            }
        }
    }
}

}  // namespace

PairCountSummary make_zero_pair_count_summary(
    const std::vector<PrimitiveChannelId>& channels
) {
    PairCountSummary summary;
    summary.channels.reserve(channels.size());

    for (const PrimitiveChannelId channel : channels) {
        static_cast<void>(primitive_channel_definition(channel));

        for (const PairChannelCount& existing : summary.channels) {
            if (existing.channel == channel) {
                throw std::invalid_argument(
                    "make_zero_pair_count_summary(): duplicate channel"
                );
            }
        }

        summary.channels.push_back({channel, 0U});
    }

    return summary;
}

void accumulate_pair_counts(
    PairCountSummary& total,
    const PairCountSummary& subevent
) {
    validate_pair_count_summary_channels(total);
    validate_pair_count_summary_channels(subevent);

    if (total.channels.size() != subevent.channels.size()) {
        throw std::invalid_argument(
            "accumulate_pair_counts(): channel-count mismatch"
        );
    }

    for (std::size_t index = 0U; index < total.channels.size(); ++index) {
        const PairChannelCount& total_entry = total.channels[index];
        const PairChannelCount& local_entry = subevent.channels[index];

        if (total_entry.channel != local_entry.channel) {
            throw std::invalid_argument(
                "accumulate_pair_counts(): channel-order mismatch"
            );
        }

        if (local_entry.pair_count >
            std::numeric_limits<std::uint64_t>::max() -
                total_entry.pair_count) {
            throw std::overflow_error(
                "accumulate_pair_counts(): pair-count overflow"
            );
        }
    }

    for (std::size_t index = 0U; index < total.channels.size(); ++index) {
        total.channels[index].pair_count += subevent.channels[index].pair_count;
    }
}

}  // namespace hbt
