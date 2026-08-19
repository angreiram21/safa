/**
 * @file emission_point_resolver.cpp
 * @brief Implementation of HBT emission-position resolution.
 */

#include "input/emission_point_resolver.h"

#include <cmath>

namespace input {
namespace {

/**
 * @brief Test whether every selected emission-position component is finite.
 * @param position Candidate `(t, x, y, z)` position.
 * @return true only when all four components satisfy std::isfinite().
 */
bool is_finite_emission_position(
    const common::FourVector& position
) noexcept {
    return std::isfinite(position.x0) &&
           std::isfinite(position.x1) &&
           std::isfinite(position.x2) &&
           std::isfinite(position.x3);
}

/**
 * @brief Build a structured result from one selected position candidate.
 * @param position Selected position candidate.
 * @param source Branch that produced the candidate.
 * @return Resolved for finite positions, otherwise NonFinitePosition.
 */
EmissionPointResolutionResult candidate_result(
    const common::FourVector& position,
    EmissionPointSource source
) {
    return {
        is_finite_emission_position(position)
            ? EmissionPointResolutionStatus::Resolved
            : EmissionPointResolutionStatus::NonFinitePosition,
        ResolvedEmissionPoint{position, source}
    };
}

/**
 * @brief Test the Afterburner time inputs required for propagation.
 * @param particle Raw Afterburner particle record.
 * @return true when the propagation-time interval is valid.
 */
bool propagation_time_inputs_valid(
    const AfterburnerParticleRecord& particle
) noexcept {
    return std::isfinite(particle.time_last_coll) &&
           particle.time_last_coll >= 0.0 &&
           particle.time_last_coll <= particle.position.x0;
}

/**
 * @brief Propagate an Afterburner position back to time_last_coll.
 * @param particle Record with valid propagation time and positive energy.
 * @return Propagated `(t, x, y, z)` position in fm.
 */
common::FourVector propagate_to_last_interaction(
    const AfterburnerParticleRecord& particle
) noexcept {
    const double dt = particle.position.x0 - particle.time_last_coll;
    const double vx = particle.momentum.x1 / particle.momentum.x0;
    const double vy = particle.momentum.x2 / particle.momentum.x0;
    const double vz = particle.momentum.x3 / particle.momentum.x0;

    return {
        particle.time_last_coll,
        particle.position.x1 - vx * dt,
        particle.position.x2 - vy * dt,
        particle.position.x3 - vz * dt
    };
}

}  // namespace

EmissionPointResolutionResult resolve_emission_point(
    int subevent_id,
    const AfterburnerParticleRecord& particle,
    bool primordial,
    const SamplerReader& sampler_reader
) {
    if (primordial && particle.ncoll == 0) {
        const std::optional<common::FourVector> sampler_position =
            sampler_reader.find_position(
                subevent_id,
                particle.id,
                particle.pdg
            );

        if (!sampler_position.has_value()) {
            return {
                EmissionPointResolutionStatus::MissingMandatorySampler,
                std::nullopt
            };
        }

        return candidate_result(
            sampler_position.value(),
            EmissionPointSource::Sampler
        );
    }

    if (propagation_time_inputs_valid(particle)) {
        return candidate_result(
            propagate_to_last_interaction(particle),
            EmissionPointSource::Propagation
        );
    }

    return candidate_result(
        particle.position,
        EmissionPointSource::Afterburner
    );
}

}  // namespace input
