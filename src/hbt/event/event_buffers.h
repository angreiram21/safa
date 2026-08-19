/**
 * @file event_buffers.h
 * @brief Per-subevent storage of accepted HBT particles by species.
 */

#ifndef SMASH_AFTERBURNER_ANALYSIS_HBT_EVENT_EVENT_BUFFERS_H
#define SMASH_AFTERBURNER_ANALYSIS_HBT_EVENT_EVENT_BUFFERS_H

#include "hbt/event/particle.h"

#include <array>
#include <cstddef>
#include <vector>

namespace hbt {

/**
 * @brief Group final accepted particles by canonical SpeciesId.
 *
 * One EventBuffers instance belongs to one independent subevent. Event
 * orchestration creates a fresh instance for each subevent and does not carry
 * its particle vectors into the next subevent.
 *
 * Required-species filtering, all particle acceptance decisions, and the
 * deterministic per-subevent input-order shuffle happen before add() is
 * called. This class only groups already accepted Particle values; it does not
 * inspect configuration, origin modes, channels, products, or pair physics.
 *
 * Storage relies on the SpeciesId contract that physical species are contiguous
 * from zero up to, but not including, SpeciesId::Count.
 */
class EventBuffers {
public:
    /**
     * @brief Add one accepted particle to its canonical species buffer.
     * @param particle Final accepted particle to store.
     * @throws std::invalid_argument If particle.species is not a physical
     *         SpeciesId value.
     */
    void add(Particle particle);

    /**
     * @brief Access the particles stored for one species.
     * @param species Canonical species whose current-subevent buffer is needed.
     * @return Immutable vector of accepted particles for the species.
     * @throws std::invalid_argument If species is not a physical SpeciesId.
     */
    [[nodiscard]] const std::vector<Particle>& get(SpeciesId species) const;

    /**
     * @brief Remove every particle from every species buffer.
     */
    void clear() noexcept;

private:
    /// Number of physical SpeciesId values stored by this class.
    static constexpr std::size_t kSpeciesBufferCount =
        static_cast<std::size_t>(SpeciesId::Count);

    /**
     * @brief Convert a physical SpeciesId to its storage-array index.
     * @param species Canonical species to index.
     * @return Zero-based array index for the species.
     * @throws std::invalid_argument If species is not a physical SpeciesId.
     */
    [[nodiscard]] static std::size_t index_for_species(SpeciesId species);

    /// Per-species accepted-particle vectors for the current subevent.
    std::array<std::vector<Particle>, kSpeciesBufferCount> buffers_{};
};

}  // namespace hbt

#endif  // SMASH_AFTERBURNER_ANALYSIS_HBT_EVENT_EVENT_BUFFERS_H
