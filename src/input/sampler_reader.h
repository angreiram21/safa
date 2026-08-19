/**
 * @file sampler_reader.h
 * @brief Indexed access to SMASH Sampler OSCAR emission positions.
 *
 * This file declares the input boundary that loads one Sampler
 * `particle_lists.oscar` file and indexes particle positions by the legacy
 * `(subevent_id, ID, PDG)` key.
 *
 * Sampler momentum is deliberately not exposed by this interface. The modular
 * HBT pipeline uses Afterburner momentum for kinematic cuts and for the final
 * accepted particle. Sampler is only a possible source of the emission
 * position resolved later by the emission-point module.
 */

#ifndef SMASH_AFTERBURNER_ANALYSIS_INPUT_SAMPLER_READER_H
#define SMASH_AFTERBURNER_ANALYSIS_INPUT_SAMPLER_READER_H

#include "common/four_vector.h"

#include <cstddef>
#include <filesystem>
#include <optional>
#include <unordered_map>

namespace input {

/**
 * @brief Loads and indexes one complete SMASH Sampler particle list.
 *
 * The constructor validates the supported OSCAR2013 particle-list schema,
 * every subevent boundary, every declared particle row, and uniqueness of the
 * legacy lookup key `(subevent_id, ID, PDG)`.
 *
 * The stored index contains only Sampler positions. Particle mass, momentum,
 * and charge are parsed to enforce the complete 12-column input schema but are
 * not retained or exposed as HBT data.
 *
 * This class does not decide whether a particle is primordial, whether a
 * Sampler lookup is required, how a missing lookup is handled by HBT, or which
 * position becomes the final HBT emission point.
 */
class SamplerReader {
public:
    /**
     * @brief Opens, validates, and indexes one Sampler OSCAR input file.
     * @param path Path to the Sampler `particle_lists.oscar` file.
     * @throws std::runtime_error If the file cannot be opened, its schema or
     *         subevent structure is malformed, a declared row cannot be read,
     *         or a duplicate `(subevent_id, ID, PDG)` key is encountered.
     */
    explicit SamplerReader(const std::filesystem::path& path);

    /**
     * @brief Finds the Sampler position associated with one legacy key.
     * @param subevent_id Subevent identifier from the matching Afterburner
     *        subevent.
     * @param id Raw SMASH particle identifier.
     * @param pdg Raw signed PDG particle code.
     * @return The indexed Sampler `(t, x, y, z)` position in fm when the key
     *         exists, otherwise `std::nullopt`.
     */
    std::optional<common::FourVector> find_position(
        int subevent_id,
        int id,
        int pdg
    ) const;

    /**
     * @brief Returns the number of parsed Sampler subevents.
     * @return Number of validated `out`/`end` subevent pairs in the file.
     */
    std::size_t subevent_count() const noexcept;

    /**
     * @brief Returns the number of parsed Sampler particle rows.
     * @return Number of validated 12-column particle rows in the file.
     */
    std::size_t particle_count() const noexcept;

private:
    /**
     * @brief Legacy identity used to index one Sampler particle position.
     */
    struct Key {
        int subevent_id;  ///< Subevent-local event identifier.
        int id;           ///< Raw SMASH particle identifier.
        int pdg;          ///< Raw signed PDG particle code.

        /**
         * @brief Tests exact equality of all three lookup-key components.
         * @param other Key to compare with this key.
         * @return `true` only when subevent, ID, and PDG all match.
         */
        bool operator==(const Key& other) const noexcept;
    };

    /**
     * @brief Hashes the three components of a Sampler lookup key.
     */
    struct KeyHash {
        /**
         * @brief Computes a hash value for one Sampler lookup key.
         * @param key Key whose components are hashed.
         * @return Combined hash value for unordered-map lookup.
         */
        std::size_t operator()(const Key& key) const noexcept;
    };

    /// Indexed Sampler positions keyed by `(subevent_id, ID, PDG)`.
    std::unordered_map<Key, common::FourVector, KeyHash> positions_;

    /// Number of validated Sampler subevents loaded from the file.
    std::size_t subevent_count_ = 0;

    /// Number of validated Sampler particle rows loaded from the file.
    std::size_t particle_count_ = 0;
};

}  // namespace input

#endif  // SMASH_AFTERBURNER_ANALYSIS_INPUT_SAMPLER_READER_H
