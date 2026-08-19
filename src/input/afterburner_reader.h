/**
 * @file afterburner_reader.h
 * @brief Streaming interface for SMASH Afterburner OSCAR input.
 *
 * This file defines the raw particle record, subevent metadata, and reader
 * interface for SMASH Afterburner `particle_lists.oscar` input.
 *
 * The input layer mirrors file data only. It does not identify HBT species,
 * classify particle origin, resolve emission points, calculate kinematics,
 * apply scientific cuts, construct pairs, or accumulate results.
 */

#ifndef SMASH_AFTERBURNER_ANALYSIS_INPUT_AFTERBURNER_READER_H
#define SMASH_AFTERBURNER_ANALYSIS_INPUT_AFTERBURNER_READER_H

#include "common/four_vector.h"

#include <cstddef>
#include <filesystem>
#include <fstream>
#include <optional>

namespace input {

/**
 * @brief Metadata declared by one Afterburner subevent opening marker.
 *
 * A SMASH OSCAR subevent begins with a marker of the form
 * `# event <id> ensemble <id> out <count>`.
 *
 * The subevent identifier is local to the current Afterburner file. The outer
 * event represented by the containing directory is intentionally not stored
 * here and remains an orchestration-level concept.
 */
struct AfterburnerSubeventHeader {
    int subevent_id;             ///< Integer following the `event` token.
    int ensemble_id;             ///< Integer following the `ensemble` token.
    std::size_t particle_count;  ///< Particle rows declared after `out`.
};

/**
 * @brief Raw particle row from SMASH OSCAR2013Extended Afterburner output.
 *
 * The record mirrors the columns and units declared by the supported
 * Afterburner `particle_lists.oscar` schema.
 *
 * The stored position is the raw Afterburner output position. It is not yet a
 * resolved HBT emission point. Likewise, the raw origin-related fields are
 * preserved without translating them into inclusive HBT origin selections.
 *
 * No HBT-specific filtering or derived kinematic quantity belongs to this
 * type.
 */
struct AfterburnerParticleRecord {
    common::FourVector position;  ///< Raw `(t, x, y, z)` position in fm.
    double mass;                  ///< Raw particle mass in GeV.
    common::FourVector momentum;  ///< Raw `(p0, px, py, pz)` in GeV.
    int pdg;                      ///< Raw dimensionless PDG particle code.
    int id;                       ///< Raw dimensionless particle identifier.
    int charge;                   ///< Raw electric charge in units of e.
    int ncoll;                    ///< Raw dimensionless collision count.
    double form_time;             ///< Raw formation time in fm.
    double xsecfac;               ///< Raw dimensionless cross-section factor.
    int proc_id_origin;           ///< Raw dimensionless origin process ID.
    int proc_type_origin;         ///< Raw dimensionless origin process type.
    double time_last_coll;        ///< Raw last-collision time in fm.
    int pdg_mother1;              ///< Raw dimensionless first-mother PDG code.
    int pdg_mother2;              ///< Raw dimensionless second-mother PDG code.
    int baryon_number;            ///< Raw dimensionless baryon number.
    int strangeness;              ///< Raw dimensionless strangeness.
};

/**
 * @brief Streams raw subevents and particle records from one Afterburner file.
 *
 * The reader owns one `particle_lists.oscar` input stream. It exposes subevent
 * boundaries explicitly so callers can process and release each subevent as
 * one independent analysis unit.
 *
 * Particle records are returned one at a time and are never accumulated by
 * this class. HBT-specific species filtering belongs to event preparation
 * after each raw record has been read.
 */
class AfterburnerReader {
public:
    /**
     * @brief Opens and validates one Afterburner OSCAR input file.
     * @param path Path to the Afterburner `particle_lists.oscar` file.
     * @throws std::runtime_error If the input file cannot be opened or its
     * OSCAR2013Extended particle-list header is missing, malformed, or
     * incompatible with the supported column and unit schema.
     */
    explicit AfterburnerReader(const std::filesystem::path& path);

    /**
     * @brief Begins the next subevent in the input stream.
     * @return Metadata from the next `out` marker, or `std::nullopt` after a
     * clean end of file.
     * @throws std::logic_error If the previous subevent is still open.
     * @throws std::runtime_error If the input is malformed or ends
     * unexpectedly while reading a subevent opening marker.
     */
    std::optional<AfterburnerSubeventHeader> begin_next_subevent();

    /**
     * @brief Reads the next raw particle row of the current subevent.
     * @return Parsed raw Afterburner particle record.
     * @throws std::logic_error If no subevent is open or all particle rows
     * declared by its `out` marker have already been read.
     * @throws std::runtime_error If the particle row is malformed or the input
     * ends unexpectedly.
     */
    AfterburnerParticleRecord read_particle();

    /**
     * @brief Reads and validates the closing marker of the current subevent.
     *
     * The closing marker must correspond to the same subevent and ensemble as
     * the opening marker. All particle rows declared by `particle_count` must
     * have been consumed before this operation is called.
     *
     * @throws std::logic_error If no subevent is open or declared particle
     * rows remain unread.
     * @throws std::runtime_error If the closing marker is missing, malformed,
     * inconsistent with the opening marker, or the input ends unexpectedly.
     */
    void finish_subevent();

private:
    std::ifstream input_;  ///< Owned Afterburner input stream.

    /// Metadata of the currently open subevent.
    AfterburnerSubeventHeader current_header_{};

    /// Number of particle rows still expected in the current subevent.
    std::size_t particles_remaining_ = 0;

    /// Whether a subevent has been opened and not yet closed.
    bool subevent_open_ = false;
};

}  // namespace input

#endif  // SMASH_AFTERBURNER_ANALYSIS_INPUT_AFTERBURNER_READER_H
