/**
 * @file rejected_pair_report.h
 * @brief In-memory storage for numerically rejected HBT particle pairs.
 */

#ifndef SMASH_AFTERBURNER_ANALYSIS_HBT_REPORTING_REJECTED_PAIR_REPORT_H
#define SMASH_AFTERBURNER_ANALYSIS_HBT_REPORTING_REJECTED_PAIR_REPORT_H

#include "common/four_vector.h"
#include "hbt/channels/channel_catalog.h"
#include "hbt/event/particle.h"
#include "hbt/pair/pair_kinematics.h"

#include <array>
#include <cstddef>
#include <cstdint>
#include <vector>

namespace hbt {

/**
 * @brief Stable numerical reason for rejecting one already formed HBT pair.
 */
enum class PairRejectionReason {
    NonFiniteKt,           ///< Calculated pair kT is not finite.
    NonFiniteMt,           ///< Calculated pair mT is not finite.
    NonFiniteDeltaTLab,    ///< Calculated Lab delta-t is not finite.
    NonFiniteDeltaTLcms,   ///< Calculated LCMS delta-t is not finite.
    NonFiniteDeltaTPrf,    ///< Calculated PRF delta-t is not finite.
    NonFiniteROutLcms,     ///< Calculated LCMS r-out is not finite.
    NonFiniteROutPrf,      ///< Calculated PRF r-out is not finite.
    NonFiniteRSide,        ///< Calculated r-side is not finite.
    NonFiniteRLong,        ///< Calculated r-long is not finite.
    NonFiniteRRadialLcms,  ///< Calculated LCMS radial radius is not finite.
    NonFiniteRRadialPrf,   ///< Calculated PRF radial radius is not finite.
    Count                  ///< Number of physical pair-rejection reasons.
};

/**
 * @brief Diagnostic particle data needed to investigate one rejected pair.
 */
struct RejectedPairParticleSnapshot {
    SpeciesId species;              ///< Canonical accepted HBT species.
    common::FourVector momentum;    ///< Final Afterburner momentum in GeV.
    double invariant_mass_gev;      ///< Validated stored invariant mass in GeV.
    int raw_pdg;                    ///< Raw signed PDG code.
    int raw_charge;                 ///< Raw electric charge in units of e.
};

/**
 * @brief Copy stable pair-rejection diagnostic data from one accepted particle.
 * @param particle Accepted HBT particle participating in the rejected pair.
 * @return Momentum, mass, species, and raw identity needed for diagnostics.
 */
[[nodiscard]] RejectedPairParticleSnapshot make_rejected_pair_particle_snapshot(
    const Particle& particle
) noexcept;

/**
 * @brief Complete diagnostic record for one numerically rejected physical pair.
 *
 * pair_ordinal_in_channel is one-based within one outer-event/subevent/channel
 * traversal. Together with outer_event_number, subevent_id, and channel it
 * identifies the rejected physical pair in the deterministic iterator order.
 */
struct RejectedPairRecord {
    std::size_t outer_event_number;  ///< One-based outer-event number.
    int subevent_id;                 ///< Afterburner subevent identifier.
    PrimitiveChannelId channel;      ///< Canonical primitive HBT channel.
    std::uint64_t pair_ordinal_in_channel;  ///< One-based traversal ordinal.
    RejectedPairParticleSnapshot particle_a;  ///< Canonical role-A snapshot.
    RejectedPairParticleSnapshot particle_b;  ///< Canonical role-B snapshot.
    PairKinematics kinematics;       ///< Calculated kinematics at rejection.
    PairRejectionReason reason;      ///< Exact numerical rejection reason.
};

/**
 * @brief Stores numerical pair rejections without performing any I/O.
 *
 * The report owns every complete rejected-pair record and maintains exact
 * aggregate counts by rejection reason. It does not decide whether a pair is
 * rejected and does not change formed-pair accounting.
 */
class RejectedPairReport {
public:
    /**
     * @brief Add one complete rejected-pair record.
     * @param record Rejection record to store.
     * @throws std::invalid_argument If record.reason is Count or invalid.
     */
    void add(RejectedPairRecord record);

    /**
     * @brief Return all stored records in insertion order.
     * @return Immutable complete rejected-pair collection.
     */
    [[nodiscard]] const std::vector<RejectedPairRecord>& records()
        const noexcept;

    /**
     * @brief Return the total number of numerically rejected pairs.
     * @return Number of stored rejection records.
     */
    [[nodiscard]] std::size_t size() const noexcept;

    /**
     * @brief Return whether the report contains no rejected pairs.
     * @return true when no rejection records are stored.
     */
    [[nodiscard]] bool empty() const noexcept;

    /**
     * @brief Return the number of records for one rejection reason.
     * @param reason Physical pair-rejection reason to count.
     * @return Number of stored records with the requested reason.
     * @throws std::invalid_argument If reason is Count or invalid.
     */
    [[nodiscard]] std::size_t count(PairRejectionReason reason) const;

private:
    /// Number of physical pair-rejection reasons.
    static constexpr std::size_t reason_count_ =
        static_cast<std::size_t>(PairRejectionReason::Count);

    std::vector<RejectedPairRecord> records_;  ///< Stored complete records.
    /// Aggregate counts indexed by physical PairRejectionReason.
    std::array<std::size_t, reason_count_> counts_{};
};

}  // namespace hbt

#endif  // SMASH_AFTERBURNER_ANALYSIS_HBT_REPORTING_REJECTED_PAIR_REPORT_H
