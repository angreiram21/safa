/**
 * @file particle_selector_test.cpp
 * @brief Unit tests for HBT single-particle kinematic acceptance.
 */

#include "hbt/event/particle_selector.h"

#include <array>
#include <cmath>
#include <cstdlib>
#include <iostream>
#include <limits>
#include <stdexcept>

namespace {

/**
 * @brief Construct one acceptance-cut set.
 * @param longitudinal_abs_max Positive longitudinal absolute-value limit.
 * @param pt_min_gev Open lower transverse-momentum bound in GeV.
 * @param pt_max_gev Open upper transverse-momentum bound in GeV.
 * @return Configured cut set.
 */
hbt::ParticleAcceptanceCuts cuts(
    double longitudinal_abs_max,
    double pt_min_gev,
    double pt_max_gev
) {
    return {longitudinal_abs_max, pt_min_gev, pt_max_gev};
}

/**
 * @brief Construct identical acceptance cuts for every particle group.
 * @param variable Configured longitudinal variable.
 * @param common_cuts Cuts copied to every species group.
 * @return Complete particle-acceptance configuration.
 */
hbt::ParticleAcceptanceConfig uniform_config(
    hbt::LongitudinalVariable variable,
    const hbt::ParticleAcceptanceCuts& common_cuts
) {
    return {
        variable,
        common_cuts,
        common_cuts,
        common_cuts,
        common_cuts,
        common_cuts
    };
}

/**
 * @brief Test whether a decision has the expected high-level status.
 * @param decision Decision returned by the selector.
 * @param expected Expected high-level status.
 * @return true when both statuses are equal.
 */
bool has_status(
    const hbt::ParticleAcceptanceDecision& decision,
    hbt::ParticleAcceptanceStatus expected
) noexcept {
    return decision.status == expected;
}

/**
 * @brief Test whether a numerical decision has the expected reason.
 * @param decision Decision returned by the selector.
 * @param expected Expected numerical reason.
 * @return true when the decision is numerical and the reason matches.
 */
bool has_numerical_reason(
    const hbt::ParticleAcceptanceDecision& decision,
    hbt::ParticleAcceptanceNumericalReason expected
) noexcept {
    return decision.status ==
               hbt::ParticleAcceptanceStatus::NumericalRejection &&
           decision.numerical_reason.has_value() &&
           decision.numerical_reason.value() == expected;
}

/**
 * @brief Verify open transverse-momentum boundaries.
 * @return true when both boundaries fail and the interior passes.
 */
bool verify_open_pt_boundaries() {
    const auto config = uniform_config(
        hbt::LongitudinalVariable::Rapidity,
        cuts(1.0, 1.0, 2.0)
    );

    const auto lower = hbt::evaluate_particle_acceptance(
        {10.0, 1.0, 0.0, 0.0}, hbt::SpeciesId::PiPlus, config);
    const auto interior = hbt::evaluate_particle_acceptance(
        {10.0, 1.5, 0.0, 0.0}, hbt::SpeciesId::PiPlus, config);
    const auto upper = hbt::evaluate_particle_acceptance(
        {10.0, 2.0, 0.0, 0.0}, hbt::SpeciesId::PiPlus, config);

    return has_status(lower, hbt::ParticleAcceptanceStatus::OutsideCuts) &&
           has_status(interior, hbt::ParticleAcceptanceStatus::Accepted) &&
           has_status(upper, hbt::ParticleAcceptanceStatus::OutsideCuts);
}

/**
 * @brief Verify rapidity acceptance and its open absolute-value boundary.
 * @return true when central, boundary, and outside cases behave correctly.
 */
bool verify_rapidity_selection() {
    const double boundary = 0.5 * std::log(4.0);
    const auto boundary_config = uniform_config(
        hbt::LongitudinalVariable::Rapidity,
        cuts(boundary, 0.5, 2.0)
    );
    const auto config = uniform_config(
        hbt::LongitudinalVariable::Rapidity,
        cuts(0.8, 0.5, 2.0)
    );

    const auto boundary_result = hbt::evaluate_particle_acceptance(
        {5.0, 1.0, 0.0, 3.0}, hbt::SpeciesId::PiPlus, boundary_config);
    const auto central = hbt::evaluate_particle_acceptance(
        {5.0, 1.0, 0.0, 0.0}, hbt::SpeciesId::PiPlus, config);
    const auto outside = hbt::evaluate_particle_acceptance(
        {2.0, 1.0, 0.0, 1.5}, hbt::SpeciesId::PiPlus, config);

    return
        has_status(
            boundary_result,
            hbt::ParticleAcceptanceStatus::OutsideCuts) &&
        has_status(central, hbt::ParticleAcceptanceStatus::Accepted) &&
        has_status(outside, hbt::ParticleAcceptanceStatus::OutsideCuts);
}

/**
 * @brief Verify pseudorapidity acceptance for both longitudinal signs.
 * @return true when central values pass and large magnitudes fail.
 */
bool verify_pseudorapidity_selection() {
    const auto config = uniform_config(
        hbt::LongitudinalVariable::Pseudorapidity,
        cuts(0.8, 0.5, 2.0)
    );

    const auto central = hbt::evaluate_particle_acceptance(
        {5.0, 1.0, 0.0, 0.0}, hbt::SpeciesId::PiPlus, config);
    const auto positive = hbt::evaluate_particle_acceptance(
        {5.0, 1.0, 0.0, 2.0}, hbt::SpeciesId::PiPlus, config);
    const auto negative = hbt::evaluate_particle_acceptance(
        {5.0, 1.0, 0.0, -2.0}, hbt::SpeciesId::PiPlus, config);

    return
        has_status(central, hbt::ParticleAcceptanceStatus::Accepted) &&
        has_status(positive, hbt::ParticleAcceptanceStatus::OutsideCuts) &&
        has_status(negative, hbt::ParticleAcceptanceStatus::OutsideCuts);
}

/**
 * @brief Pair one physical species with the pT accepted only by its group.
 */
struct SpeciesAcceptanceCase {
    hbt::SpeciesId species;  ///< Species whose group mapping is tested.
    double pt_gev;           ///< Group-specific accepted transverse momentum.
};

/**
 * @brief Verify the fixed SpeciesId-to-acceptance-group mapping.
 * @return true when all physical species use the intended cut group.
 */
bool verify_species_group_mapping() {
    const hbt::ParticleAcceptanceConfig config{
        hbt::LongitudinalVariable::Rapidity,
        cuts(1.0, 0.5, 1.5),
        cuts(1.0, 1.5, 2.5),
        cuts(1.0, 2.5, 3.5),
        cuts(1.0, 3.5, 4.5),
        cuts(1.0, 4.5, 5.5)
    };

    const std::array<SpeciesAcceptanceCase, 16> cases{{
        {hbt::SpeciesId::PiPlus, 1.0},
        {hbt::SpeciesId::PiMinus, 1.0},
        {hbt::SpeciesId::PiZero, 1.0},
        {hbt::SpeciesId::KPlus, 2.0},
        {hbt::SpeciesId::KMinus, 2.0},
        {hbt::SpeciesId::KZero, 2.0},
        {hbt::SpeciesId::KZeroBar, 2.0},
        {hbt::SpeciesId::Proton, 3.0},
        {hbt::SpeciesId::ProtonBar, 3.0},
        {hbt::SpeciesId::Neutron, 3.0},
        {hbt::SpeciesId::NeutronBar, 3.0},
        {hbt::SpeciesId::SigmaPlus, 4.0},
        {hbt::SpeciesId::SigmaBarMinus, 4.0},
        {hbt::SpeciesId::SigmaZero, 4.0},
        {hbt::SpeciesId::Lambda, 5.0},
        {hbt::SpeciesId::LambdaBar, 5.0}
    }};

    for (const SpeciesAcceptanceCase& test_case : cases) {
        const auto result = hbt::evaluate_particle_acceptance(
            {10.0, test_case.pt_gev, 0.0, 0.0},
            test_case.species,
            config
        );

        if (!has_status(result, hbt::ParticleAcceptanceStatus::Accepted)) {
            std::cerr
                << "particle_selector_test: species used wrong cut group.\n";
            return false;
        }
    }

    return true;
}

/**
 * @brief Verify recoverable numerical-rejection reasons.
 * @return true when representative numerical failures return exact reasons.
 */
bool verify_numerical_rejections() {
    const double nan = std::numeric_limits<double>::quiet_NaN();
    const double max = std::numeric_limits<double>::max();
    const double below_max = std::nextafter(max, 0.0);
    const auto rapidity_config = uniform_config(
        hbt::LongitudinalVariable::Rapidity,
        cuts(1000.0, 0.0, max)
    );
    const auto eta_config = uniform_config(
        hbt::LongitudinalVariable::Pseudorapidity,
        cuts(1000.0, 0.0, max)
    );

    const auto nonfinite_momentum = hbt::evaluate_particle_acceptance(
        {5.0, nan, 0.0, 0.0}, hbt::SpeciesId::PiPlus, rapidity_config);
    const auto negative_energy = hbt::evaluate_particle_acceptance(
        {-5.0, 1.0, 0.0, 0.0}, hbt::SpeciesId::PiPlus, eta_config);
    const auto zero_energy = hbt::evaluate_particle_acceptance(
        {0.0, 1.0, 0.0, 0.0}, hbt::SpeciesId::PiPlus, rapidity_config);
    const auto nonfinite_pt = hbt::evaluate_particle_acceptance(
        {max, max, max, 0.0}, hbt::SpeciesId::PiPlus, rapidity_config);
    const auto invalid_rapidity = hbt::evaluate_particle_acceptance(
        {1.0, 0.5, 0.0, 1.0}, hbt::SpeciesId::PiPlus, rapidity_config);
    const auto nonfinite_rapidity = hbt::evaluate_particle_acceptance(
        {max, 1.0, 0.0, below_max},
        hbt::SpeciesId::PiPlus,
        rapidity_config
    );
    const auto invalid_eta = hbt::evaluate_particle_acceptance(
        {5.0, 0.0, 0.0, 1.0}, hbt::SpeciesId::PiPlus, eta_config);
    const auto nonfinite_eta = hbt::evaluate_particle_acceptance(
        {max, 1.0, 0.0, max}, hbt::SpeciesId::PiPlus, eta_config);

    return
        has_numerical_reason(
            nonfinite_momentum,
            hbt::ParticleAcceptanceNumericalReason::NonFiniteMomentum) &&
        has_numerical_reason(
            negative_energy,
            hbt::ParticleAcceptanceNumericalReason::NonPositiveEnergy) &&
        negative_energy.diagnostic_value.has_value() &&
        negative_energy.diagnostic_value.value() == -5.0 &&
        has_numerical_reason(
            zero_energy,
            hbt::ParticleAcceptanceNumericalReason::NonPositiveEnergy) &&
        zero_energy.diagnostic_value.has_value() &&
        zero_energy.diagnostic_value.value() == 0.0 &&
        has_numerical_reason(
            nonfinite_pt,
            hbt::ParticleAcceptanceNumericalReason::
                NonFiniteTransverseMomentum) &&
        has_numerical_reason(
            invalid_rapidity,
            hbt::ParticleAcceptanceNumericalReason::InvalidRapidityInput) &&
        has_numerical_reason(
            nonfinite_rapidity,
            hbt::ParticleAcceptanceNumericalReason::NonFiniteRapidity) &&
        has_numerical_reason(
            invalid_eta,
            hbt::ParticleAcceptanceNumericalReason::
                InvalidPseudorapidityInput) &&
        has_numerical_reason(
            nonfinite_eta,
            hbt::ParticleAcceptanceNumericalReason::
                NonFinitePseudorapidity);
}

/**
 * @brief Verify invalid enum values remain structural API errors.
 * @return true when invalid species and longitudinal enum values throw.
 */
bool verify_invalid_enums_throw() {
    const auto config = uniform_config(
        hbt::LongitudinalVariable::Rapidity,
        cuts(1.0, 0.5, 2.0)
    );
    bool species_threw = false;
    bool variable_threw = false;

    try {
        static_cast<void>(hbt::evaluate_particle_acceptance(
            {5.0, 1.0, 0.0, 0.0},
            hbt::SpeciesId::Count,
            config
        ));
    } catch (const std::invalid_argument&) {
        species_threw = true;
    }

    hbt::ParticleAcceptanceConfig invalid = config;
    invalid.longitudinal_variable =
        static_cast<hbt::LongitudinalVariable>(999);

    try {
        static_cast<void>(hbt::evaluate_particle_acceptance(
            {5.0, 1.0, 0.0, 0.0},
            hbt::SpeciesId::PiPlus,
            invalid
        ));
    } catch (const std::invalid_argument&) {
        variable_threw = true;
    }

    return species_threw && variable_threw;
}

}  // namespace

/**
 * @brief Run the complete HBT particle-acceptance test collection.
 * @return EXIT_SUCCESS when every test passes, otherwise EXIT_FAILURE.
 */
int main() {
    bool success = true;
    success = verify_open_pt_boundaries() && success;
    success = verify_rapidity_selection() && success;
    success = verify_pseudorapidity_selection() && success;
    success = verify_species_group_mapping() && success;
    success = verify_numerical_rejections() && success;
    success = verify_invalid_enums_throw() && success;
    return success ? EXIT_SUCCESS : EXIT_FAILURE;
}
