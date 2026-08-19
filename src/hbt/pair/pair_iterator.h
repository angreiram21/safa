/**
 * @file pair_iterator.h
 * @brief Streaming iteration over physical pairs in one HBT subevent.
 */

#ifndef SMASH_AFTERBURNER_ANALYSIS_HBT_PAIR_PAIR_ITERATOR_H
#define SMASH_AFTERBURNER_ANALYSIS_HBT_PAIR_PAIR_ITERATOR_H

#include "hbt/channels/channel_catalog.h"
#include "hbt/event/event_buffers.h"
#include "hbt/pair/pair_count_summary.h"

#include <cstddef>
#include <stdexcept>
#include <vector>

namespace hbt {
namespace detail {

/**
 * @brief Validate and resolve channels before pair iteration begins.
 * @param required_channels Primitive channels requested for this subevent.
 * @return Catalogue definitions in the same order as @p required_channels.
 * @throws std::invalid_argument If a channel identifier is invalid or repeated.
 *
 * Every channel definition is resolved exactly once before any pair is
 * delivered to the consumer. This prevents structural input errors from
 * producing a partially processed subevent and avoids catalogue lookups in the
 * pair-emission hot path.
 */
inline std::vector<const PrimitiveChannel*> resolve_pair_iterator_channels(
    const std::vector<PrimitiveChannelId>& required_channels
) {
    std::vector<const PrimitiveChannel*> definitions;
    definitions.reserve(required_channels.size());

    for (std::size_t index = 0; index < required_channels.size(); ++index) {
        const PrimitiveChannelId channel = required_channels[index];

        for (std::size_t previous = 0; previous < index; ++previous) {
            if (required_channels[previous] == channel) {
                throw std::invalid_argument(
                    "for_each_pair(): duplicate PrimitiveChannelId"
                );
            }
        }

        definitions.push_back(&primitive_channel_definition(channel));
    }

    return definitions;
}

}  // namespace detail

/**
 * @brief Deliver every physical pair for the required primitive channels.
 * @tparam PairConsumer Callable accepting the ordered channel index, channel,
 *         and two Particle references.
 * @param buffers Accepted particles from exactly one current subevent.
 * @param required_channels Unique primitive channels to process, preserving
 *        output order.
 * @param consumer Callable invoked once for every physical particle pair.
 * @return Pair counts for every requested channel, preserving channel order.
 * @throws std::invalid_argument If a channel identifier is invalid or repeated.
 *
 * Identical-species channels produce each unordered pair exactly once using
 * indices i < j. Cross-species channels produce the complete Cartesian product
 * while preserving the canonical species-A/species-B roles from the channel
 * catalogue. Empty buffers and channels with too few particles are valid and
 * produce a zero count.
 *
 * Pairs are streamed directly to the consumer and are never stored by this
 * function. The function performs no pair kinematics, slicing, histogramming,
 * origin routing, or output. Errors raised by the consumer are not caught or
 * suppressed.
 */
template <typename PairConsumer>
PairCountSummary for_each_pair(
    const EventBuffers& buffers,
    const std::vector<PrimitiveChannelId>& required_channels,
    PairConsumer&& consumer
) {
    const std::vector<const PrimitiveChannel*> definitions =
        detail::resolve_pair_iterator_channels(required_channels);

    PairCountSummary summary;
    summary.channels.reserve(required_channels.size());

    for (std::size_t channel_index = 0U;
         channel_index < required_channels.size();
         ++channel_index) {
        const PrimitiveChannelId channel = required_channels[channel_index];
        const PrimitiveChannel& definition = *definitions[channel_index];
        const std::vector<Particle>& particles_a =
            buffers.get(definition.species_a);

        summary.channels.push_back({channel, 0U});
        std::uint64_t& pair_count = summary.channels.back().pair_count;

        if (definition.species_a == definition.species_b) {
            for (std::size_t i = 0; i < particles_a.size(); ++i) {
                for (std::size_t j = i + 1; j < particles_a.size(); ++j) {
                    consumer(
                        channel_index,
                        channel,
                        particles_a[i],
                        particles_a[j]
                    );
                    ++pair_count;
                }
            }
            continue;
        }

        const std::vector<Particle>& particles_b =
            buffers.get(definition.species_b);

        for (const Particle& particle_a : particles_a) {
            for (const Particle& particle_b : particles_b) {
                consumer(
                    channel_index,
                    channel,
                    particle_a,
                    particle_b
                );
                ++pair_count;
            }
        }
    }

    return summary;
}

}  // namespace hbt

#endif  // SMASH_AFTERBURNER_ANALYSIS_HBT_PAIR_PAIR_ITERATOR_H
