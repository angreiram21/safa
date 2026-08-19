/**
 * @file rejected_particle_report.h
 * @brief In-memory storage for numerically rejected HBT particles.
 */

#ifndef SMASH_AFTERBURNER_ANALYSIS_HBT_REPORTING_REJECTED_PARTICLE_REPORT_H
#define SMASH_AFTERBURNER_ANALYSIS_HBT_REPORTING_REJECTED_PARTICLE_REPORT_H

#include "common/four_vector.h"
#include "hbt/species/species.h"

#include <array>
#include <cstddef>
#include <optional>
#include <vector>

namespace hbt {

/**
 * @brief Stable numerical reason for rejecting one HBT particle.
 */
enum class ParticleRejectionReason {
    NonFiniteMomentum,  ///< At least one momentum component is invalid.
    NonPositiveEnergy,  ///< Raw Afterburner energy is zero or negative.
    NonFiniteTransverseMomentum, ///< Calculated pT is not finite.
    InvalidRapidityInput,        ///< Momentum is outside the rapidity domain.
    NonFiniteRapidity,           ///< Calculated rapidity is not finite.
    InvalidPseudorapidityInput,  ///< Momentum is outside the eta domain.
    NonFinitePseudorapidity,     ///< Calculated pseudorapidity is not finite.
    NonFiniteInvariantMassSquared, ///< Calculated invariant m^2 is not finite.
    NonPositiveInvariantMassSquared, ///< Invariant m^2 is zero or negative.
    NonFiniteInvariantMass,      ///< Calculated invariant mass is not finite.
    NonFiniteSamplerEmissionPosition, ///< Sampler position is not finite.
    NonFinitePropagationEmissionPosition, ///< Propagated position is invalid.
    NonFiniteAfterburnerEmissionPosition, ///< Raw selected position is invalid.
    Count                        ///< Number of physical rejection reasons.
};

/**
 * @brief Complete diagnostic record for one numerically rejected particle.
 */
struct RejectedParticleRecord {
    std::size_t outer_event_number;  ///< One-based outer-event number.
    int subevent_id;                 ///< Afterburner subevent identifier.
    int particle_id;                 ///< Raw Afterburner particle identifier.
    int raw_pdg;                     ///< Raw signed PDG code.
    int raw_charge;                  ///< Raw electric charge in units of e.
    SpeciesId species;               ///< Canonical identified HBT species.
    common::FourVector momentum;     ///< Raw Afterburner four-momentum in GeV.
    common::FourVector raw_position; ///< Raw Afterburner position in fm.
    double raw_mass_gev;             ///< Raw OSCAR mass column in GeV.
    int ncoll;                       ///< Raw Afterburner collision count.
    double time_last_coll;           ///< Raw last-collision time in fm.
    ParticleRejectionReason reason;  ///< Exact numerical rejection reason.
    /// Invalid calculated quantity when a specific scalar value exists.
    std::optional<double> diagnostic_value;
    /// Invalid selected emission position when position resolution failed.
    std::optional<common::FourVector> diagnostic_position;
};

/**
 * @brief Stores numerical particle rejections without performing any I/O.
 *
 * The report owns every complete rejected-particle record and maintains exact
 * aggregate counts by rejection reason. It does not print, serialize, write
 * files, or decide whether a particle should be rejected.
 */
class RejectedParticleReport {
public:
    /**
     * @brief Add one complete rejected-particle record.
     * @param record Rejection record to store.
     * @throws std::invalid_argument If record.reason is Count or invalid.
     */
    void add(RejectedParticleRecord record);

    /**
     * @brief Return all stored records in insertion order.
     * @return Immutable complete rejection-record collection.
     */
    [[nodiscard]] const std::vector<RejectedParticleRecord>& records()
        const noexcept;

    /**
     * @brief Return the total number of rejected particles.
     * @return Number of stored rejection records.
     */
    [[nodiscard]] std::size_t size() const noexcept;

    /**
     * @brief Return whether the report contains no rejected particles.
     * @return true when no rejection records are stored.
     */
    [[nodiscard]] bool empty() const noexcept;

    /**
     * @brief Return the number of records for one rejection reason.
     * @param reason Physical rejection reason to count.
     * @return Number of stored records with the requested reason.
     * @throws std::invalid_argument If reason is Count or invalid.
     */
    [[nodiscard]] std::size_t count(
        ParticleRejectionReason reason) const;

private:
    /// Number of physical particle-rejection reasons.
    static constexpr std::size_t reason_count_ =
        static_cast<std::size_t>(ParticleRejectionReason::Count);

    std::vector<RejectedParticleRecord> records_;  ///< Stored full records.
    /// Aggregate counts indexed by physical ParticleRejectionReason.
    std::array<std::size_t, reason_count_> counts_{};
};

}  // namespace hbt

#endif
