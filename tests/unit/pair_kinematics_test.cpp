/**
 * @file pair_kinematics_test.cpp
 * @brief Unit tests for reusable pair-level transverse kinematics.
 */

#include "hbt/pair/pair_kinematics.h"

#include <cmath>
#include <cstdlib>
#include <iostream>

namespace {

constexpr double kTolerance = 1.0e-12;

/**
 * @brief Construct one accepted particle for pair-kinematics tests.
 * @param px Final Afterburner px component in GeV.
 * @param py Final Afterburner py component in GeV.
 * @param mass_gev Stored validated invariant mass in GeV.
 * @return Complete particle containing the requested transverse data.
 */
hbt::Particle make_particle(double px, double py, double mass_gev) {
    return {
        hbt::SpeciesId::PiPlus,
        {},
        {10.0, px, py, 7.0},
        mass_gev,
        {true, true, true},
        211,
        1
    };
}

/**
 * @brief Compare two floating-point values with the test tolerance.
 * @param actual Calculated value.
 * @param expected Reference value.
 * @return true when the absolute difference is within tolerance.
 */
bool nearly_equal(double actual, double expected) {
    return std::abs(actual - expected) <= kTolerance;
}

/**
 * @brief Report one failed pair-kinematics condition.
 * @param message Description of the violated contract.
 * @return Always false.
 */
bool fail(const char* message) {
    std::cerr << "pair_kinematics_test: " << message << ".\n";
    return false;
}

/**
 * @brief Verify the total pair four-momentum and K components are retained.
 * @return true when all pair sums and derived K components are exact.
 */
bool verify_reusable_pair_state() {
    const hbt::PairKinematics result = hbt::calculate_pair_kinematics(
        make_particle(3.0, 4.0, 0.5),
        make_particle(1.0, -2.0, 0.7)
    );

    if (!nearly_equal(result.pair_four_momentum.x0, 20.0) ||
        !nearly_equal(result.pair_four_momentum.x1, 4.0) ||
        !nearly_equal(result.pair_four_momentum.x2, 2.0) ||
        !nearly_equal(result.pair_four_momentum.x3, 14.0) ||
        !nearly_equal(result.kx_gev, 2.0) ||
        !nearly_equal(result.ky_gev, 1.0)) {
        return fail("reusable pair four-momentum or K components differ");
    }
    return true;
}

/**
 * @brief Verify opposite transverse momenta produce zero pair kT.
 * @return true when kT is zero within floating-point tolerance.
 */
bool verify_zero_kt_for_opposite_transverse_momenta() {
    const hbt::PairKinematics result = hbt::calculate_pair_kinematics(
        make_particle(1.0, 2.0, 0.5),
        make_particle(-1.0, -2.0, 0.5)
    );

    if (!nearly_equal(result.kt_gev, 0.0)) {
        return fail("opposite transverse momenta did not produce zero kT");
    }
    return true;
}

/**
 * @brief Verify kT against an analytic transverse-momentum example.
 * @return true when the calculated value equals sqrt(8) GeV.
 */
bool verify_known_kt() {
    const hbt::PairKinematics result = hbt::calculate_pair_kinematics(
        make_particle(3.0, 4.0, 0.5),
        make_particle(1.0, 0.0, 0.5)
    );

    if (!nearly_equal(result.kt_gev, std::sqrt(8.0))) {
        return fail("analytic kT value differs");
    }
    return true;
}

/**
 * @brief Verify mT for two particles with the same stored mass.
 * @return true when mT follows sqrt(m^2 + kT^2).
 */
bool verify_mt_with_equal_masses() {
    const hbt::PairKinematics result = hbt::calculate_pair_kinematics(
        make_particle(0.8, 0.0, 0.6),
        make_particle(0.8, 0.0, 0.6)
    );

    if (!nearly_equal(result.kt_gev, 0.8) ||
        !nearly_equal(result.mt_gev, 1.0)) {
        return fail("equal-mass mT value differs");
    }
    return true;
}

/**
 * @brief Verify mT uses the arithmetic average of distinct particle masses.
 * @return true when masses 0.4 and 0.8 GeV produce m_avg = 0.6 GeV.
 */
bool verify_mt_with_distinct_masses() {
    const hbt::PairKinematics result = hbt::calculate_pair_kinematics(
        make_particle(0.8, 0.0, 0.4),
        make_particle(0.8, 0.0, 0.8)
    );

    if (!nearly_equal(result.kt_gev, 0.8) ||
        !nearly_equal(result.mt_gev, 1.0)) {
        return fail("distinct-mass mT value differs");
    }
    return true;
}

/**
 * @brief Verify a representable tiny nonzero K_T is not rounded to zero.
 * @return true when stable kT preserves the nonzero transverse magnitude.
 */
bool verify_tiny_nonzero_kt_remains_nonzero() {
    const hbt::PairKinematics result = hbt::calculate_pair_kinematics(
        make_particle(1.0e-200, 0.0, 0.5),
        make_particle(1.0e-200, 0.0, 0.5)
    );

    if (!(result.kt_gev > 0.0) ||
        !std::isfinite(result.kt_gev) ||
        result.kt_gev != std::hypot(result.kx_gev, result.ky_gev)) {
        return fail("tiny finite K_T was rounded to zero");
    }
    return true;
}

/**
 * @brief Verify a representable large mT survives intermediate squares.
 * @return true when finite mass and kT produce the finite Euclidean norm.
 */
bool verify_large_finite_mt_remains_finite() {
    constexpr double scale = 1.0e154;
    const hbt::PairKinematics result = hbt::calculate_pair_kinematics(
        make_particle(scale, 0.0, scale),
        make_particle(scale, 0.0, scale)
    );
    const double expected = std::hypot(scale, scale);

    if (!std::isfinite(result.kt_gev) ||
        !std::isfinite(result.mt_gev) ||
        !nearly_equal(result.mt_gev, expected)) {
        return fail("large representable mT became non-finite");
    }
    return true;
}

/**
 * @brief Verify pair kinematics are symmetric under particle exchange.
 * @return true when swapping A and B preserves both calculated values.
 */
bool verify_particle_exchange_symmetry() {
    const hbt::Particle particle_a = make_particle(1.2, -0.4, 0.3);
    const hbt::Particle particle_b = make_particle(-0.2, 1.0, 0.9);

    const hbt::PairKinematics ab =
        hbt::calculate_pair_kinematics(particle_a, particle_b);
    const hbt::PairKinematics ba =
        hbt::calculate_pair_kinematics(particle_b, particle_a);

    if (!nearly_equal(ab.pair_four_momentum.x0,
                      ba.pair_four_momentum.x0) ||
        !nearly_equal(ab.pair_four_momentum.x1,
                      ba.pair_four_momentum.x1) ||
        !nearly_equal(ab.pair_four_momentum.x2,
                      ba.pair_four_momentum.x2) ||
        !nearly_equal(ab.pair_four_momentum.x3,
                      ba.pair_four_momentum.x3) ||
        !nearly_equal(ab.kx_gev, ba.kx_gev) ||
        !nearly_equal(ab.ky_gev, ba.ky_gev) ||
        !nearly_equal(ab.kt_gev, ba.kt_gev) ||
        !nearly_equal(ab.mt_gev, ba.mt_gev)) {
        return fail("particle exchange changed pair kinematics");
    }
    return true;
}

/**
 * @brief Verify mT reuses stored masses instead of reconstructing them.
 * @return true when intentionally inconsistent energies do not affect mT.
 */
bool verify_stored_mass_is_reused() {
    hbt::Particle particle_a = make_particle(0.8, 0.0, 0.4);
    hbt::Particle particle_b = make_particle(0.8, 0.0, 0.8);

    particle_a.momentum.x0 = 100.0;
    particle_a.momentum.x3 = 90.0;
    particle_b.momentum.x0 = 2.0;
    particle_b.momentum.x3 = 1.5;

    const hbt::PairKinematics result =
        hbt::calculate_pair_kinematics(particle_a, particle_b);

    if (!nearly_equal(result.kt_gev, 0.8) ||
        !nearly_equal(result.mt_gev, 1.0)) {
        return fail("mT did not reuse stored particle masses");
    }
    return true;
}

}  // namespace

/**
 * @brief Run all pair-kinematics unit tests.
 * @return EXIT_SUCCESS when every test passes, otherwise EXIT_FAILURE.
 */
int main() {
    bool success = true;
    success = verify_reusable_pair_state() && success;
    success = verify_zero_kt_for_opposite_transverse_momenta() && success;
    success = verify_known_kt() && success;
    success = verify_mt_with_equal_masses() && success;
    success = verify_mt_with_distinct_masses() && success;
    success = verify_tiny_nonzero_kt_remains_nonzero() && success;
    success = verify_large_finite_mt_remains_finite() && success;
    success = verify_particle_exchange_symmetry() && success;
    success = verify_stored_mass_is_reused() && success;

    if (!success) {
        return EXIT_FAILURE;
    }
    return EXIT_SUCCESS;
}
