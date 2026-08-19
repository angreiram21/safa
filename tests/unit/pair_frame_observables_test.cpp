/**
 * @file pair_frame_observables_test.cpp
 * @brief Unit tests for composed Lab, LCMS, OSL, and PRF pair observables.
 */

#include "hbt/pair/pair_frame_observables.h"

#include <cmath>
#include <cstdlib>
#include <iostream>

namespace {

constexpr double kTolerance = 1.0e-12;

/**
 * @brief Compare two floating-point values with scaled tolerance.
 * @param actual Calculated value.
 * @param expected Reference value.
 * @return true when values agree within the test tolerance.
 */
bool nearly_equal(double actual, double expected) {
    const double scale = 1.0 + std::abs(expected);
    return std::abs(actual - expected) <= kTolerance * scale;
}

/**
 * @brief Report one failed frame-observable condition.
 * @param message Description of the violated contract.
 * @return Always false.
 */
bool fail(const char* message) {
    std::cerr << "pair_frame_observables_test: " << message << ".\n";
    return false;
}

/**
 * @brief Build one accepted particle for composed-frame tests.
 * @param position Emission position (t, x, y, z) in fm.
 * @param momentum Final momentum (E, px, py, pz) in GeV.
 * @param mass_gev Validated stored invariant mass in GeV.
 * @return Complete accepted particle.
 */
hbt::Particle make_particle(
    common::FourVector position,
    common::FourVector momentum,
    double mass_gev
) {
    return {
        hbt::SpeciesId::PiPlus,
        position,
        momentum,
        mass_gev,
        {true, true, true},
        211,
        1
    };
}

/**
 * @brief Verify one fully nontrivial Lab-to-LCMS-to-PRF transformation.
 * @return true when all nine retained observables match analytic formulas.
 */
bool verify_nontrivial_composed_transformation() {
    const hbt::Particle particle_a = make_particle(
        {5.0, 3.0, 4.0, 2.0},
        {6.0, 3.0, 0.0, 2.0},
        std::sqrt(23.0)
    );
    const hbt::Particle particle_b = make_particle(
        {1.0, 1.0, 1.0, -1.0},
        {4.0, 1.0, 0.0, 2.0},
        std::sqrt(11.0)
    );
    const hbt::PairKinematics kinematics =
        hbt::calculate_pair_kinematics(particle_a, particle_b);
    const hbt::PairFrameObservables result =
        hbt::calculate_pair_frame_observables(
            particle_a,
            particle_b,
            kinematics
        );

    const double beta_lcms = 0.4;
    const double gamma_lcms = 1.0 / std::sqrt(1.0 - beta_lcms * beta_lcms);
    const double delta_t_lcms = gamma_lcms * (4.0 - beta_lcms * 3.0);
    const double r_long = gamma_lcms * (3.0 - beta_lcms * 4.0);
    const double beta_out = 4.0 * gamma_lcms / 10.0;
    const double gamma_out = 1.0 / std::sqrt(1.0 - beta_out * beta_out);
    const double delta_t_prf = gamma_out * (
        delta_t_lcms - beta_out * 2.0
    );
    const double r_out_prf = gamma_out * (
        2.0 - beta_out * delta_t_lcms
    );
    const double radial_lcms = std::sqrt(
        2.0 * 2.0 + 3.0 * 3.0 + r_long * r_long
    );
    const double radial_prf = std::sqrt(
        r_out_prf * r_out_prf + 3.0 * 3.0 + r_long * r_long
    );

    if (!nearly_equal(result.delta_t_lab_fm, 4.0) ||
        !nearly_equal(result.delta_t_lcms_fm, delta_t_lcms) ||
        !nearly_equal(result.delta_t_prf_fm, delta_t_prf) ||
        !nearly_equal(result.r_out_lcms_fm, 2.0) ||
        !nearly_equal(result.r_out_prf_fm, r_out_prf) ||
        !nearly_equal(result.r_side_fm, 3.0) ||
        !nearly_equal(result.r_long_fm, r_long) ||
        !nearly_equal(result.r_radial_lcms_fm, radial_lcms) ||
        !nearly_equal(result.r_radial_prf_fm, radial_prf)) {
        return fail("nontrivial composed transformation differs");
    }
    return true;
}

/**
 * @brief Verify the normal OSL side orientation for non-axis-aligned K_T.
 * @return true when out and side use the agreed right-handed convention.
 */
bool verify_normal_osl_orientation() {
    const hbt::Particle particle_a = make_particle(
        {0.0, 5.0, 10.0, 0.0},
        {12.0, 4.0, 5.0, 0.0},
        std::sqrt(103.0)
    );
    const hbt::Particle particle_b = make_particle(
        {0.0, 0.0, 0.0, 0.0},
        {12.0, 2.0, 3.0, 0.0},
        std::sqrt(131.0)
    );
    const hbt::PairKinematics kinematics =
        hbt::calculate_pair_kinematics(particle_a, particle_b);
    const hbt::PairFrameObservables result =
        hbt::calculate_pair_frame_observables(
            particle_a,
            particle_b,
            kinematics
        );

    if (!nearly_equal(kinematics.kx_gev, 3.0) ||
        !nearly_equal(kinematics.ky_gev, 4.0) ||
        !nearly_equal(kinematics.kt_gev, 5.0) ||
        !nearly_equal(result.r_out_lcms_fm, 11.0) ||
        !nearly_equal(result.r_side_fm, 2.0) ||
        !nearly_equal(result.r_long_fm, 0.0)) {
        return fail("normal OSL orientation differs");
    }
    return true;
}

/**
 * @brief Verify exact kT zero uses qT for OSL and makes PRF equal LCMS.
 * @return true when qT orientation and both identity boosts are correct.
 */
bool verify_zero_kt_uses_qt_fallback() {
    const hbt::Particle particle_a = make_particle(
        {4.0, 5.0, 1.0, 2.0},
        {5.0, 1.0, 2.0, 0.0},
        std::sqrt(20.0)
    );
    const hbt::Particle particle_b = make_particle(
        {1.0, 0.0, 0.0, 0.0},
        {5.0, -1.0, -2.0, 0.0},
        std::sqrt(20.0)
    );
    const hbt::PairKinematics kinematics =
        hbt::calculate_pair_kinematics(particle_a, particle_b);
    const hbt::PairFrameObservables result =
        hbt::calculate_pair_frame_observables(
            particle_a,
            particle_b,
            kinematics
        );

    const double sqrt_five = std::sqrt(5.0);
    const double r_out = 7.0 / sqrt_five;
    const double r_side = -9.0 / sqrt_five;
    const double radial = std::sqrt(
        r_out * r_out + r_side * r_side + 4.0
    );

    if (kinematics.kt_gev != 0.0 ||
        !nearly_equal(result.delta_t_lab_fm, 3.0) ||
        !nearly_equal(result.delta_t_lcms_fm, 3.0) ||
        !nearly_equal(result.delta_t_prf_fm, 3.0) ||
        !nearly_equal(result.r_out_lcms_fm, r_out) ||
        !nearly_equal(result.r_out_prf_fm, r_out) ||
        !nearly_equal(result.r_side_fm, r_side) ||
        !nearly_equal(result.r_long_fm, 2.0) ||
        !nearly_equal(result.r_radial_lcms_fm, radial) ||
        !nearly_equal(result.r_radial_prf_fm, radial)) {
        return fail("exact kT zero did not use the qT OSL convention");
    }
    return true;
}

/**
 * @brief Verify large finite qT preserves the exact kT-zero OSL basis.
 * @return true when qT normalization remains finite and unit length.
 */
bool verify_large_finite_qt_preserves_osl_basis() {
    constexpr double transverse = 9.0e153;
    const hbt::Particle particle_a = make_particle(
        {0.0, 1.0, 2.0, 0.0},
        {1.0e154, transverse, 0.0, 0.0},
        1.0
    );
    const hbt::Particle particle_b = make_particle(
        {0.0, 0.0, 0.0, 0.0},
        {1.0e154, -transverse, 0.0, 0.0},
        1.0
    );
    const hbt::PairKinematics kinematics =
        hbt::calculate_pair_kinematics(particle_a, particle_b);
    const hbt::PairFrameObservables result =
        hbt::calculate_pair_frame_observables(
            particle_a,
            particle_b,
            kinematics
        );

    if (kinematics.kt_gev != 0.0 ||
        !nearly_equal(result.r_out_lcms_fm, 1.0) ||
        !nearly_equal(result.r_side_fm, 2.0) ||
        !std::isfinite(result.r_radial_lcms_fm)) {
        return fail("large finite qT corrupted the kT-zero OSL basis");
    }
    return true;
}

/**
 * @brief Verify large finite spatial components keep finite radial radii.
 * @return true when both LCMS and PRF Euclidean radii remain representable.
 */
bool verify_large_finite_radial_radii_remain_finite() {
    constexpr double separation = 1.0e154;
    const hbt::Particle particle_a = make_particle(
        {0.0, separation, separation, 0.0},
        {2.0, 1.0, 0.0, 0.0},
        1.0
    );
    const hbt::Particle particle_b = make_particle(
        {0.0, 0.0, 0.0, 0.0},
        {2.0, 1.0, 0.0, 0.0},
        1.0
    );
    const hbt::PairKinematics kinematics =
        hbt::calculate_pair_kinematics(particle_a, particle_b);
    const hbt::PairFrameObservables result =
        hbt::calculate_pair_frame_observables(
            particle_a,
            particle_b,
            kinematics
        );
    const double expected_lcms = std::sqrt(2.0) * separation;
    const double expected_prf = std::sqrt(7.0 / 3.0) * separation;

    if (!std::isfinite(result.r_radial_lcms_fm) ||
        !std::isfinite(result.r_radial_prf_fm) ||
        !nearly_equal(result.r_radial_lcms_fm, expected_lcms) ||
        !nearly_equal(result.r_radial_prf_fm, expected_prf)) {
        return fail("large representable radial radius became non-finite");
    }
    return true;
}

}  // namespace

/**
 * @brief Run all composed pair-frame-observable unit tests.
 * @return EXIT_SUCCESS when every test passes, otherwise EXIT_FAILURE.
 */
int main() {
    bool success = true;
    success = verify_nontrivial_composed_transformation() && success;
    success = verify_normal_osl_orientation() && success;
    success = verify_zero_kt_uses_qt_fallback() && success;
    success = verify_large_finite_qt_preserves_osl_basis() && success;
    success = verify_large_finite_radial_radii_remain_finite() && success;
    return success ? EXIT_SUCCESS : EXIT_FAILURE;
}
