/**
 * @file particle_selector.cpp
 * @brief Implementation of HBT single-particle kinematic acceptance.
 */

#include "hbt/event/particle_selector.h"

#include "common/kinematics.h"
#include "common/kinematics_validation.h"

#include <cmath>
#include <stdexcept>

namespace hbt {
namespace {

/**
 * @brief Return the configured acceptance cuts for one particle species.
 * @param species Canonical HBT particle species.
 * @param config Complete particle-acceptance configuration.
 * @return Reference to the species group's configured acceptance cuts.
 * @throws std::invalid_argument If species is not a physical SpeciesId.
 */
const ParticleAcceptanceCuts& acceptance_cuts_for_species(
    SpeciesId species,
    const ParticleAcceptanceConfig& config
) {
    switch (species) {
        case SpeciesId::PiPlus:
        case SpeciesId::PiMinus:
        case SpeciesId::PiZero:
            return config.pions;

        case SpeciesId::KPlus:
        case SpeciesId::KMinus:
        case SpeciesId::KZero:
        case SpeciesId::KZeroBar:
            return config.kaons;

        case SpeciesId::Proton:
        case SpeciesId::ProtonBar:
        case SpeciesId::Neutron:
        case SpeciesId::NeutronBar:
            return config.nucleons;

        case SpeciesId::SigmaPlus:
        case SpeciesId::SigmaBarMinus:
        case SpeciesId::SigmaZero:
            return config.sigmas;

        case SpeciesId::Lambda:
        case SpeciesId::LambdaBar:
            return config.lambdas;

        case SpeciesId::Count:
            break;
    }

    throw std::invalid_argument(
        "evaluate_particle_acceptance(): invalid SpeciesId"
    );
}

/**
 * @brief Build a numerical-rejection decision without a calculated value.
 * @param reason Numerical cause of the rejection.
 * @return Numerical-rejection decision with no diagnostic value.
 */
ParticleAcceptanceDecision numerical_rejection(
    ParticleAcceptanceNumericalReason reason
) {
    return {
        ParticleAcceptanceStatus::NumericalRejection,
        reason,
        std::nullopt
    };
}

/**
 * @brief Build a numerical-rejection decision with a calculated value.
 * @param reason Numerical cause of the rejection.
 * @param value Invalid calculated quantity associated with the rejection.
 * @return Numerical-rejection decision containing the diagnostic value.
 */
ParticleAcceptanceDecision numerical_rejection(
    ParticleAcceptanceNumericalReason reason,
    double value
) {
    return {
        ParticleAcceptanceStatus::NumericalRejection,
        reason,
        value
    };
}

}  // namespace

ParticleAcceptanceDecision evaluate_particle_acceptance(
    const common::FourVector& momentum,
    SpeciesId species,
    const ParticleAcceptanceConfig& config
) {
    const ParticleAcceptanceCuts& cuts =
        acceptance_cuts_for_species(species, config);

    if (!common::is_finite_four_momentum(momentum)) {
        return numerical_rejection(
            ParticleAcceptanceNumericalReason::NonFiniteMomentum
        );
    }

    if (momentum.x0 <= 0.0) {
        return numerical_rejection(
            ParticleAcceptanceNumericalReason::NonPositiveEnergy,
            momentum.x0
        );
    }

    const double transverse_momentum =
        common::transverse_momentum(momentum);

    if (!common::is_finite_kinematic_result(transverse_momentum)) {
        return numerical_rejection(
            ParticleAcceptanceNumericalReason::NonFiniteTransverseMomentum,
            transverse_momentum
        );
    }

    double longitudinal_value = 0.0;

    switch (config.longitudinal_variable) {
        case LongitudinalVariable::Rapidity:
            if (std::abs(momentum.x3) >= momentum.x0) {
                return numerical_rejection(
                    ParticleAcceptanceNumericalReason::InvalidRapidityInput
                );
            }

            longitudinal_value = common::rapidity(momentum);

            if (!common::is_finite_kinematic_result(longitudinal_value)) {
                return numerical_rejection(
                    ParticleAcceptanceNumericalReason::NonFiniteRapidity,
                    longitudinal_value
                );
            }
            break;

        case LongitudinalVariable::Pseudorapidity:
            if (transverse_momentum == 0.0) {
                return numerical_rejection(
                    ParticleAcceptanceNumericalReason::
                        InvalidPseudorapidityInput
                );
            }

            longitudinal_value = common::pseudorapidity(momentum);

            if (!common::is_finite_kinematic_result(longitudinal_value)) {
                return numerical_rejection(
                    ParticleAcceptanceNumericalReason::
                        NonFinitePseudorapidity,
                    longitudinal_value
                );
            }
            break;

        default:
            throw std::invalid_argument(
                "evaluate_particle_acceptance(): invalid "
                "LongitudinalVariable"
            );
    }

    const bool accepted =
        std::abs(longitudinal_value) < cuts.longitudinal_abs_max &&
        cuts.pt_min_gev < transverse_momentum &&
        transverse_momentum < cuts.pt_max_gev;

    return {
        accepted ? ParticleAcceptanceStatus::Accepted
                 : ParticleAcceptanceStatus::OutsideCuts,
        std::nullopt,
        std::nullopt
    };
}

}  // namespace hbt
