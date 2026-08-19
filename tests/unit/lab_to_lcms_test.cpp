/**
 * @file lab_to_lcms_test.cpp
 * @brief Unit tests for the HBT Lab-to-LCMS transformation.
 */

#include "hbt/boosts/lab_to_lcms.h"

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
 * @brief Report one failed Lab-to-LCMS test condition.
 * @param message Description of the violated contract.
 * @return Always false.
 */
bool fail(const char* message) {
    std::cerr << "lab_to_lcms_test: " << message << ".\n";
    return false;
}

/**
 * @brief Build reusable pair kinematics for one boost test.
 * @param pair Total Lab-frame pair four-momentum.
 * @param kt_gev Pair kT in GeV.
 * @return Complete kinematics with transverse values consistent with pair.
 */
hbt::PairKinematics make_kinematics(
    common::FourVector pair,
    double kt_gev
) {
    return {
        pair,
        0.5 * pair.x1,
        0.5 * pair.x2,
        kt_gev,
        kt_gev
    };
}

/**
 * @brief Verify the selected legacy-sign longitudinal boost formula.
 * @return true when all transformed components and beta_out are correct.
 */
bool verify_nonzero_longitudinal_boost() {
    const hbt::PairKinematics kinematics = make_kinematics(
        {10.0, 4.0, 0.0, 6.0},
        2.0
    );
    const common::FourVector separation{5.0, 2.0, -3.0, 4.0};
    const hbt::LabToLCMSResult result = hbt::boost_lab_to_lcms(
        separation,
        kinematics
    );

    if (!nearly_equal(result.relative_separation_fm.x0, 3.25) ||
        !nearly_equal(result.relative_separation_fm.x1, 2.0) ||
        !nearly_equal(result.relative_separation_fm.x2, -3.0) ||
        !nearly_equal(result.relative_separation_fm.x3, 1.25) ||
        !nearly_equal(result.beta_out, 0.5)) {
        return fail("nonzero longitudinal boost differs");
    }
    return true;
}

/**
 * @brief Verify exact beta_LCMS zero reuses the Lab separation unchanged.
 * @return true when the identity branch and beta_out are exact.
 */
bool verify_exact_zero_is_identity() {
    const hbt::PairKinematics kinematics = make_kinematics(
        {10.0, 4.0, 0.0, 0.0},
        2.0
    );
    const common::FourVector separation{5.0, 2.0, -3.0, 4.0};
    const hbt::LabToLCMSResult result = hbt::boost_lab_to_lcms(
        separation,
        kinematics
    );

    if (result.relative_separation_fm.x0 != separation.x0 ||
        result.relative_separation_fm.x1 != separation.x1 ||
        result.relative_separation_fm.x2 != separation.x2 ||
        result.relative_separation_fm.x3 != separation.x3 ||
        !nearly_equal(result.beta_out, 0.4)) {
        return fail("exact zero beta was not the identity");
    }
    return true;
}

/**
 * @brief Verify a tiny nonzero beta is not replaced by a legacy threshold.
 * @return true when the linear beta term remains active for beta = 1e-16.
 */
bool verify_tiny_nonzero_beta_is_transformed() {
    const hbt::PairKinematics kinematics = make_kinematics(
        {1.0e16, 0.0, 0.0, 1.0},
        0.0
    );
    const common::FourVector separation{2.0, 0.0, 0.0, 1.0e16};
    const hbt::LabToLCMSResult result = hbt::boost_lab_to_lcms(
        separation,
        kinematics
    );

    if (!nearly_equal(result.relative_separation_fm.x0, 1.0) ||
        result.relative_separation_fm.x0 == separation.x0) {
        return fail("tiny nonzero beta was treated as zero");
    }
    return true;
}

}  // namespace

/**
 * @brief Run all Lab-to-LCMS unit tests.
 * @return EXIT_SUCCESS when every test passes, otherwise EXIT_FAILURE.
 */
int main() {
    bool success = true;
    success = verify_nonzero_longitudinal_boost() && success;
    success = verify_exact_zero_is_identity() && success;
    success = verify_tiny_nonzero_beta_is_transformed() && success;
    return success ? EXIT_SUCCESS : EXIT_FAILURE;
}
