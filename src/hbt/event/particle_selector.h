/**
 * @file particle_selector.h
 * @brief Single-particle kinematic acceptance for HBT event preparation.
 *
 * This file declares HBT particle-acceptance evaluation that distinguishes
 * configured-cut rejection from recoverable numerical rejection.
 */

#ifndef SMASH_AFTERBURNER_ANALYSIS_HBT_EVENT_PARTICLE_SELECTOR_H
#define SMASH_AFTERBURNER_ANALYSIS_HBT_EVENT_PARTICLE_SELECTOR_H

#include "common/four_vector.h"
#include "hbt/config/hbt_config.h"
#include "hbt/species/species.h"

#include <optional>

namespace hbt {

/**
 * @brief High-level outcome of HBT single-particle acceptance evaluation.
 */
enum class ParticleAcceptanceStatus {
    Accepted,           ///< Particle passed all configured acceptance cuts.
    OutsideCuts,        ///< Valid particle failed a configured physical cut.
    NumericalRejection  ///< Particle cannot be evaluated numerically.
};

/**
 * @brief Numerical reason for rejecting one particle during acceptance.
 */
enum class ParticleAcceptanceNumericalReason {
    NonFiniteMomentum,  ///< At least one momentum component is invalid.
    NonPositiveEnergy,  ///< Energy is zero or negative.
    NonFiniteTransverseMomentum, ///< Calculated pT is not finite.
    InvalidRapidityInput,        ///< Momentum is outside the rapidity domain.
    NonFiniteRapidity,           ///< Calculated rapidity is not finite.
    InvalidPseudorapidityInput,  ///< Momentum is outside the eta domain.
    NonFinitePseudorapidity      ///< Calculated pseudorapidity is not finite.
};

/**
 * @brief Result of evaluating HBT single-particle kinematic acceptance.
 *
 * numerical_reason is engaged only for NumericalRejection. diagnostic_value is
 * engaged when a specific invalid scalar is useful for diagnostics.
 */
struct ParticleAcceptanceDecision {
    ParticleAcceptanceStatus status;  ///< High-level acceptance outcome.
    /// Numerical failure reason, present only for NumericalRejection.
    std::optional<ParticleAcceptanceNumericalReason> numerical_reason;
    /// Invalid diagnostic scalar when one exists.
    std::optional<double> diagnostic_value;
};

/**
 * @brief Evaluate one identified particle against configured HBT cuts.
 *
 * The supplied momentum is the raw Afterburner momentum selected by event
 * orchestration. Numerical failures, including non-positive energy, are
 * returned as NumericalRejection rather than aborting the run. A particle with
 * non-positive energy is rejected before any acceptance observable is
 * calculated. Valid values are tested against the open boundaries
 *
 *     |longitudinal variable| < longitudinal_abs_max
 *     pt_min_gev < pT < pt_max_gev.
 *
 * This operation does not identify species, inspect required_species, classify
 * origin, calculate invariant mass, access Sampler, resolve emission points,
 * construct final Particle values, or build pairs.
 *
 * @param momentum Four-momentum used for HBT single-particle acceptance.
 * @param species Already identified canonical HBT particle species.
 * @param config Validated particle-acceptance configuration.
 * @return Structured acceptance decision and any numerical-failure detail.
 * @throws std::invalid_argument If species or longitudinal_variable contains
 *         an invalid enum value.
 */
[[nodiscard]] ParticleAcceptanceDecision evaluate_particle_acceptance(
    const common::FourVector& momentum,
    SpeciesId species,
    const ParticleAcceptanceConfig& config
);

}  // namespace hbt

#endif  // SMASH_AFTERBURNER_ANALYSIS_HBT_EVENT_PARTICLE_SELECTOR_H
