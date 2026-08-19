/**
 * @file pair_origin_routing.cpp
 * @brief Pair-origin membership and requested-origin routing implementation.
 */

#include "hbt/pair/pair_origin_routing.h"

#include <stdexcept>

namespace hbt {

PairOriginMemberships calculate_pair_origin_memberships(
    const OriginFlags& particle_a,
    const OriginFlags& particle_b
) noexcept {
    return {
        particle_a.primordial && particle_b.primordial,
        particle_a.primordial_rescattering &&
            particle_b.primordial_rescattering,
        particle_a.primordial_rescattering_decay &&
            particle_b.primordial_rescattering_decay
    };
}

PairOriginRoutes route_pair_origin_memberships(
    const PairOriginMemberships& memberships,
    OriginMode mode
) {
    switch (mode) {
    case OriginMode::Primordial:
        return {memberships.primordial, false, false};
    case OriginMode::PrimordialRescattering:
        return {false, memberships.primordial_rescattering, false};
    case OriginMode::PrimordialRescatteringDecay:
        return {false, false, memberships.primordial_rescattering_decay};
    case OriginMode::All:
        return {
            memberships.primordial,
            memberships.primordial_rescattering,
            memberships.primordial_rescattering_decay
        };
    }

    throw std::invalid_argument("invalid HBT origin mode for pair routing");
}

}  // namespace hbt
