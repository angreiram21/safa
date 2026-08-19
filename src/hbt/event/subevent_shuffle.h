/**
 * @file subevent_shuffle.h
 * @brief Deterministic per-subevent decorrelation of accepted HBT particles.
 */

#ifndef SMASH_AFTERBURNER_ANALYSIS_HBT_EVENT_SUBEVENT_SHUFFLE_H
#define SMASH_AFTERBURNER_ANALYSIS_HBT_EVENT_SUBEVENT_SHUFFLE_H

#include "hbt/event/particle.h"

#include <cstddef>
#include <vector>

namespace hbt {

/**
 * @brief Shuffle accepted HBT particles for one subevent deterministically.
 * @param particles Complete accepted and required particles of one subevent.
 * @param outer_event_number One-based outer-event number.
 * @param subevent_id Afterburner subevent identifier.
 *
 * SMASH input ordering can carry non-physical correlations. Pair construction
 * assigns roles through container order, so signed pair observables would
 * otherwise inherit that ordering. This operation removes that dependence by
 * permuting the accepted particles once, independently for every subevent,
 * before species grouping and pair formation.
 *
 * The complete Particle values are permuted in place; particle identity and
 * every stored scientific field remain attached to the same Particle. The
 * deterministic seed depends only on outer_event_number and subevent_id. Each
 * call owns its random-number generator, so the permutation is independent of
 * worker assignment and configured thread count and shares no mutable state.
 */
void shuffle_subevent_particles(
    std::vector<Particle>& particles,
    std::size_t outer_event_number,
    int subevent_id
);

}  // namespace hbt

#endif  // SMASH_AFTERBURNER_ANALYSIS_HBT_EVENT_SUBEVENT_SHUFFLE_H
