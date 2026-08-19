/**
 * @file required_species.h
 * @brief Derivation of particle species required by primitive HBT channels.
 *
 * This file declares the operation that derives the unique particle species
 * required by a collection of primitive HBT channels.
 *
 * The operation resolves channel-to-species dependencies only. It does not
 * select analysis products, construct particle pairs, filter event particles,
 * calculate observables, allocate histograms, or perform analysis processing.
 */

#ifndef HBT_SELECTION_REQUIRED_SPECIES_H
#define HBT_SELECTION_REQUIRED_SPECIES_H

#include "hbt/channels/primitive_channel.h"

#include <vector>

namespace hbt {

    /**
     * @brief Derive the unique particle species required by primitive HBT
     * channels.
     *
     * Primitive channels are inspected in their supplied order. For each
     * channel, the canonical species-A role is inspected before species B.
     *
     * Each SpeciesId appears at most once in the returned vector. When a
     * species occurs in several channels, its first occurrence determines its
     * position in the result.
     *
     * Identical-particle channels contribute their species only once.
     *
     * Preserving first-occurrence order provides deterministic species ordering
     * while respecting the canonical A/B ordering of each primitive channel.
     *
     * The returned species collection is intended to be resolved before event
     * processing begins, so that the HBT input path can reject particles whose
     * species cannot contribute to any requested primitive channel.
     *
     * This function performs dependency derivation only. The actual
     * HBT-specific
     * particle filtering belongs to the event-preparation layer.
     *
     * An empty channel collection produces an empty result.
     *
     * @param channels Primitive HBT channels whose particle species are needed.
     *
     * @return Unique particle species required by the supplied primitive
     *         channels, ordered by first occurrence.
     *
     * @throws std::invalid_argument If any channel identifier is not a valid
     *         PrimitiveChannelId value.
     */
    std::vector<SpeciesId> required_species(
        const std::vector<PrimitiveChannelId>& channels
    );

}  // namespace hbt

#endif  // HBT_SELECTION_REQUIRED_SPECIES_H
