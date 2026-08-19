/**
 * @file pair_kinematics_validation_test.cpp
 * @brief Unit tests for pure HBT pair-kinematics numerical validation.
 */

#include "hbt/pair/pair_kinematics_validation.h"

#include <cstdlib>
#include <iostream>
#include <limits>

namespace {

/**
 * @brief Build pair kinematics with selected kT and mT test values.
 * @param kt_gev kT value to store.
 * @param mt_gev mT value to store.
 * @return Complete PairKinematics for validation tests.
 */
hbt::PairKinematics make_kinematics(double kt_gev, double mt_gev) {
    return {{1.0, 0.0, 0.0, 0.0}, 0.0, 0.0, kt_gev, mt_gev};
}

/**
 * @brief Report one failed pair-kinematics validation condition.
 * @param message Description of the violated contract.
 * @return Always false.
 */
bool fail(const char* message) {
    std::cerr << "pair_kinematics_validation_test: " << message << ".\n";
    return false;
}

/**
 * @brief Verify finite kT and mT values are accepted independently.
 * @return true when both finite checks return true.
 */
bool verify_finite_values_are_accepted() {
    const hbt::PairKinematics kinematics = make_kinematics(0.4, 0.8);

    if (!hbt::is_finite_pair_kt(kinematics) ||
        !hbt::is_finite_pair_mt(kinematics)) {
        return fail("finite pair kinematics were rejected");
    }
    return true;
}

/**
 * @brief Verify non-finite kT is rejected without inspecting mT policy.
 * @return true when infinite and NaN kT values are rejected.
 */
bool verify_non_finite_kt_is_rejected() {
    const double infinity = std::numeric_limits<double>::infinity();
    const double nan = std::numeric_limits<double>::quiet_NaN();

    if (hbt::is_finite_pair_kt(make_kinematics(infinity, 1.0)) ||
        hbt::is_finite_pair_kt(make_kinematics(nan, 1.0))) {
        return fail("non-finite kT was accepted");
    }
    return true;
}

/**
 * @brief Verify non-finite mT is rejected independently from finite kT.
 * @return true when infinite and NaN mT values are rejected.
 */
bool verify_non_finite_mt_is_rejected() {
    const double infinity = std::numeric_limits<double>::infinity();
    const double nan = std::numeric_limits<double>::quiet_NaN();

    if (hbt::is_finite_pair_mt(make_kinematics(1.0, infinity)) ||
        hbt::is_finite_pair_mt(make_kinematics(1.0, nan))) {
        return fail("non-finite mT was accepted");
    }
    return true;
}

}  // namespace

/**
 * @brief Run all pair-kinematics validation unit tests.
 * @return EXIT_SUCCESS when every test passes, otherwise EXIT_FAILURE.
 */
int main() {
    bool success = true;
    success = verify_finite_values_are_accepted() && success;
    success = verify_non_finite_kt_is_rejected() && success;
    success = verify_non_finite_mt_is_rejected() && success;
    return success ? EXIT_SUCCESS : EXIT_FAILURE;
}
