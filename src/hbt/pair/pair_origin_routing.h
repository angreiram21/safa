/**
 * @file pair_origin_routing.h
 * @brief Pair-origin membership and requested-origin routing.
 */

#ifndef SMASH_AFTERBURNER_ANALYSIS_HBT_PAIR_PAIR_ORIGIN_ROUTING_H
#define SMASH_AFTERBURNER_ANALYSIS_HBT_PAIR_PAIR_ORIGIN_ROUTING_H

#include "hbt/config/hbt_config.h"
#include "hbt/event/origin_flags.h"

namespace hbt {

/**
 * @brief Nested physical origin memberships of one already formed pair.
 *
 * Memberships are inclusive. A pair belongs to one selection only when both
 * particles belong to that same particle-level selection.
 */
struct PairOriginMemberships {
    bool primordial;                     ///< Both particles are primordial.
    bool primordial_rescattering;        ///< Both belong to the middle slice.
    bool primordial_rescattering_decay;  ///< Both belong to the widest slice.
};

/**
 * @brief Requested origin routes activated for one physical pair.
 *
 * These fields describe output routing, not additional physical categories.
 * OriginMode::All may activate more than one route for the same pair.
 */
struct PairOriginRoutes {
    bool primordial;                     ///< Route to primordial results.
    bool primordial_rescattering;        ///< Route to middle results.
    bool primordial_rescattering_decay;  ///< Route to widest results.
};

/**
 * @brief Calculate nested physical origin memberships for one pair.
 * @param particle_a Origin memberships of the first particle.
 * @param particle_b Origin memberships of the second particle.
 * @return Pair memberships obtained by intersecting equal nested selections.
 *
 * The exact contract is:
 *
 * - P-P belongs to all three selections;
 * - P-R, R-P, and R-R belong to the middle and widest selections;
 * - any pair containing D belongs only to the widest selection.
 *
 * This function performs no pair formation, kinematic calculation, slicing,
 * histogramming, or output.
 */
[[nodiscard]] PairOriginMemberships calculate_pair_origin_memberships(
    const OriginFlags& particle_a,
    const OriginFlags& particle_b
) noexcept;

/**
 * @brief Select requested origin routes from physical pair memberships.
 * @param memberships Physical nested memberships of the pair.
 * @param mode Configured origin mode for the HBT analysis.
 * @return Routes requested by mode and compatible with the pair.
 * @throws std::invalid_argument If mode is not a valid OriginMode value.
 *
 * Individual modes activate only their corresponding route. OriginMode::All
 * activates every compatible nested route simultaneously. The function does
 * not recalculate pair observables or perform downstream analysis work.
 */
[[nodiscard]] PairOriginRoutes route_pair_origin_memberships(
    const PairOriginMemberships& memberships,
    OriginMode mode
);

}  // namespace hbt

#endif  // SMASH_AFTERBURNER_ANALYSIS_HBT_PAIR_PAIR_ORIGIN_ROUTING_H
