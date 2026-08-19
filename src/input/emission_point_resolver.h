/**
 * @file emission_point_resolver.h
 * @brief Resolution of the HBT emission position from prepared input data.
 *
 * This file declares the input-adapter operation that selects the emission
 * position from Sampler data, propagation of the raw Afterburner position, or
 * the raw Afterburner position itself.
 */

#ifndef SMASH_AFTERBURNER_ANALYSIS_INPUT_EMISSION_POINT_RESOLVER_H
#define SMASH_AFTERBURNER_ANALYSIS_INPUT_EMISSION_POINT_RESOLVER_H

#include "input/afterburner_reader.h"
#include "input/sampler_reader.h"

#include <optional>

namespace input {

/**
 * @brief Source used to obtain one emission-position candidate.
 */
enum class EmissionPointSource {
    Sampler,      ///< Position obtained from the matching Sampler record.
    Propagation,  ///< Afterburner position propagated to time_last_coll.
    Afterburner   ///< Raw Afterburner position used without propagation.
};

/**
 * @brief One emission-position candidate and its provenance.
 */
struct ResolvedEmissionPoint {
    common::FourVector position;  ///< Candidate `(t, x, y, z)` in fm.
    EmissionPointSource source;   ///< Branch that produced the candidate.
};

/**
 * @brief Outcome of one emission-position resolution attempt.
 */
enum class EmissionPointResolutionStatus {
    Resolved,  ///< Candidate exists and every component is finite.
    MissingMandatorySampler, ///< Required Sampler key was absent.
    NonFinitePosition  ///< Selected candidate contains a non-finite value.
};

/**
 * @brief Structured result of one emission-position resolution attempt.
 *
 * candidate is engaged for Resolved and NonFinitePosition. It is empty only
 * for MissingMandatorySampler. A non-finite candidate is returned unchanged so
 * the orchestration layer can reject and report that particle without fallback.
 */
struct EmissionPointResolutionResult {
    EmissionPointResolutionStatus status;  ///< Resolution outcome.
    /// Selected candidate and source when a branch produced a position.
    std::optional<ResolvedEmissionPoint> candidate;
};

/**
 * @brief Resolve the emission position for one already accepted particle.
 *
 * Resolution follows this precedence:
 *
 * 1. A primordial particle with `ncoll == 0` requires a Sampler lookup by
 *    `(subevent_id, ID, PDG)`. A missing required match returns
 *    MissingMandatorySampler; no fallback to Afterburner is permitted.
 * 2. Otherwise, propagation to `time_last_coll` is used when
 *    `time_last_coll` is finite and non-negative, and
 *    `time_last_coll <= t_afterburner`.
 * 3. If propagation is unavailable, the raw Afterburner position is selected.
 *
 * A selected non-finite position returns NonFinitePosition with the unchanged
 * candidate and its source. It is never repaired and never replaced by another
 * source.
 *
 * @param subevent_id Subevent identifier used in the Sampler lookup key.
 * @param particle Raw Afterburner record whose position may be resolved.
 * @param primordial Whether OriginSelector classified the particle primordial.
 * @param sampler_reader Sampler index for the same outer event.
 * @return Structured resolution outcome and selected candidate when present.
 * @pre particle was accepted by HBT particle acceptance and therefore has a
 *      finite four-momentum with strictly positive energy.
 */
[[nodiscard]] EmissionPointResolutionResult resolve_emission_point(
    int subevent_id,
    const AfterburnerParticleRecord& particle,
    bool primordial,
    const SamplerReader& sampler_reader
);

}  // namespace input

#endif  // SMASH_AFTERBURNER_ANALYSIS_INPUT_EMISSION_POINT_RESOLVER_H
