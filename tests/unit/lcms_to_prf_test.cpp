/**
 * @file lcms_to_prf_test.cpp
 * @brief Unit tests for the HBT LCMS-to-PRF out-direction boost.
 */

#include "hbt/boosts/lcms_to_prf.h"

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
 * @brief Report one failed LCMS-to-PRF test condition.
 * @param message Description of the violated contract.
 * @return Always false.
 */
bool fail(const char* message) {
    std::cerr << "lcms_to_prf_test: " << message << ".\n";
    return false;
}

/**
 * @brief Verify the one-dimensional temporal/out Lorentz transformation.
 * @return true when both transformed components match analytic values.
 */
bool verify_nonzero_out_boost() {
    const hbt::LCMSToPRFResult result = hbt::boost_lcms_to_prf(
        5.0,
        4.0,
        0.6
    );

    if (!nearly_equal(result.delta_t_prf_fm, 3.25) ||
        !nearly_equal(result.r_out_prf_fm, 1.25)) {
        return fail("nonzero out boost differs");
    }
    return true;
}

/**
 * @brief Verify the boost removes the total pair out momentum in the PRF.
 * @return true when a four-momentum-like input has zero transformed out.
 */
bool verify_pair_rest_frame_definition() {
    const hbt::LCMSToPRFResult result = hbt::boost_lcms_to_prf(
        10.0,
        6.0,
        0.6
    );

    if (!nearly_equal(result.r_out_prf_fm, 0.0) ||
        !nearly_equal(result.delta_t_prf_fm, 8.0)) {
        return fail("boost did not produce the pair rest frame");
    }
    return true;
}

/**
 * @brief Verify exact beta_out zero makes PRF identical to LCMS.
 * @return true when both values are returned bit-for-bit unchanged.
 */
bool verify_exact_zero_is_identity() {
    const hbt::LCMSToPRFResult result = hbt::boost_lcms_to_prf(
        5.0,
        -4.0,
        0.0
    );

    if (result.delta_t_prf_fm != 5.0 ||
        result.r_out_prf_fm != -4.0) {
        return fail("exact zero beta was not the identity");
    }
    return true;
}

/**
 * @brief Verify a tiny nonzero beta_out is never thresholded to zero.
 * @return true when its linear temporal correction remains active.
 */
bool verify_tiny_nonzero_beta_is_transformed() {
    const hbt::LCMSToPRFResult result = hbt::boost_lcms_to_prf(
        2.0,
        1.0e16,
        1.0e-16
    );

    if (!nearly_equal(result.delta_t_prf_fm, 1.0) ||
        result.delta_t_prf_fm == 2.0) {
        return fail("tiny nonzero beta was treated as zero");
    }
    return true;
}

}  // namespace

/**
 * @brief Run all LCMS-to-PRF unit tests.
 * @return EXIT_SUCCESS when every test passes, otherwise EXIT_FAILURE.
 */
int main() {
    bool success = true;
    success = verify_nonzero_out_boost() && success;
    success = verify_pair_rest_frame_definition() && success;
    success = verify_exact_zero_is_identity() && success;
    success = verify_tiny_nonzero_beta_is_transformed() && success;
    return success ? EXIT_SUCCESS : EXIT_FAILURE;
}
