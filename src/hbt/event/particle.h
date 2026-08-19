/**
 * @file particle.h
 * @brief Final accepted particle data consumed by later HBT stages.
 */

#ifndef SMASH_AFTERBURNER_ANALYSIS_HBT_EVENT_PARTICLE_H
#define SMASH_AFTERBURNER_ANALYSIS_HBT_EVENT_PARTICLE_H

#include "common/four_vector.h"
#include "hbt/event/origin_flags.h"
#include "hbt/species/species.h"

namespace hbt {

/**
 * @brief Final accepted HBT particle after event preparation.
 *
 * position is the resolved emission position in fm. momentum is always the
 * Afterburner four-momentum in GeV that was used for particle acceptance.
 * invariant_mass_gev is calculated once from that same momentum after its
 * squared value has been validated as finite and strictly positive.
 *
 * Sampler data may determine position but never replaces momentum. origin
 * preserves all nested origin memberships for later slice routing. raw_pdg and
 * raw_charge preserve the identifying Afterburner values for diagnostics.
 *
 * This structure contains no pair observables, histogram state, product IDs,
 * or information about which source produced the resolved position.
 */
struct Particle {
    SpeciesId species;             ///< Canonical HBT particle species.
    common::FourVector position;   ///< Resolved `(t, x, y, z)` in fm.
    common::FourVector momentum;   ///< Afterburner `(p0, px, py, pz)` in GeV.
    double invariant_mass_gev;     ///< Validated invariant mass in GeV.
    OriginFlags origin;            ///< Nested origin memberships.
    int raw_pdg;                   ///< Raw signed Afterburner PDG code.
    int raw_charge;                ///< Raw Afterburner electric charge.
};

}  // namespace hbt

#endif  // SMASH_AFTERBURNER_ANALYSIS_HBT_EVENT_PARTICLE_H
