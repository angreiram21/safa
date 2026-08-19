/**
 * @file origin_flags.h
 * @brief Nested HBT origin memberships carried by accepted particles.
 */

#ifndef SMASH_AFTERBURNER_ANALYSIS_HBT_EVENT_ORIGIN_FLAGS_H
#define SMASH_AFTERBURNER_ANALYSIS_HBT_EVENT_ORIGIN_FLAGS_H

namespace hbt {

/**
 * @brief Nested origin memberships carried by an accepted HBT particle.
 *
 * The flags represent three inclusive selections. A primordial contribution
 * belongs to all three. A rescattering contribution belongs to the two wider
 * selections. A contribution present only in the widest selection has only
 * primordial_rescattering_decay set.
 *
 * These flags describe origin membership only. They do not determine particle
 * acceptance, emission-position resolution, pair construction, or routing to
 * analysis products.
 */
struct OriginFlags {
    bool primordial;                     ///< Most restrictive membership.
    bool primordial_rescattering;        ///< Inclusive middle membership.
    bool primordial_rescattering_decay;  ///< Widest inclusive membership.
};

}  // namespace hbt

#endif  // SMASH_AFTERBURNER_ANALYSIS_HBT_EVENT_ORIGIN_FLAGS_H
