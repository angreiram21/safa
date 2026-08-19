/**
 * @file kinematics_test.cpp
 * @brief Unit tests for common single-particle kinematic calculations.
 *
 * This test verifies the numerical results of the calculation-only kinematics
 * module. Input-domain and result validation are tested separately by
 * kinematics_validation_test.
 */

#include "common/kinematics.h"

#include <cmath>
#include <cstdlib>
#include <iostream>
#include <string_view>

namespace {

    constexpr double relative_tolerance = 1.0e-12;

    /**
     * @brief Compare two finite values with a relative numerical tolerance.
     *
     * @param actual Calculated value.
     * @param expected Expected reference value.
     * @return true when the values agree within the test tolerance.
     */
    bool nearly_equal(double actual, double expected) noexcept {
        const double scale = 1.0 + std::abs(expected);
        return std::abs(actual - expected) <= relative_tolerance * scale;
    }

    /**
     * @brief Verify one calculated value against its expected result.
     *
     * @param actual Calculated value.
     * @param expected Expected reference value.
     * @param case_name Human-readable name used in failure diagnostics.
     * @return true when the value is finite and agrees with the reference.
     */
    bool expect_near(
        double actual,
        double expected,
        std::string_view case_name) noexcept {
        if (!std::isfinite(actual) || !nearly_equal(actual, expected)) {
            std::cerr
                << "kinematics_test: "
                << case_name
                << " failed; expected "
                << expected
                << ", got "
                << actual
                << ".\n";
            return false;
        }

        return true;
    }

    /**
     * @brief Verify transverse-momentum calculations.
     *
     * @return true when all transverse-momentum cases pass.
     */
    bool verify_transverse_momentum() noexcept {
        const double large_expected = std::sqrt(2.0) * 1.0e308;

        return
            expect_near(
                common::transverse_momentum(
                    common::FourVector{7.0, 0.0, 0.0, 9.0}),
                0.0,
                "zero transverse momentum") &&
            expect_near(
                common::transverse_momentum(
                    common::FourVector{7.0, 3.0, 4.0, 9.0}),
                5.0,
                "3-4-5 transverse momentum") &&
            expect_near(
                common::transverse_momentum(
                    common::FourVector{-7.0, -3.0, -4.0, -9.0}),
                5.0,
                "transverse-momentum sign invariance") &&
            expect_near(
                common::transverse_momentum(
                    common::FourVector{0.0, 1.0e308, 1.0e308, 0.0}),
                large_expected,
                "large transverse momentum");
    }

    /**
     * @brief Verify three-momentum-magnitude calculations.
     *
     * @return true when all three-momentum-magnitude cases pass.
     */
    bool verify_momentum_magnitude() noexcept {
        const double large_expected = std::sqrt(3.0) * 1.0e308;

        return
            expect_near(
                common::momentum_magnitude(
                    common::FourVector{7.0, 0.0, 0.0, 0.0}),
                0.0,
                "zero three-momentum magnitude") &&
            expect_near(
                common::momentum_magnitude(
                    common::FourVector{7.0, 3.0, 4.0, 12.0}),
                13.0,
                "3-4-12 three-momentum magnitude") &&
            expect_near(
                common::momentum_magnitude(
                    common::FourVector{-7.0, -3.0, -4.0, -12.0}),
                13.0,
                "three-momentum-magnitude sign invariance") &&
            expect_near(
                common::momentum_magnitude(
                    common::FourVector{
                        0.0, 1.0e308, 1.0e308, 1.0e308}),
                large_expected,
                "large three-momentum magnitude");
    }

    /**
     * @brief Verify longitudinal-rapidity calculations.
     *
     * @return true when all rapidity cases pass.
     */
    bool verify_rapidity() noexcept {
        const double log_two = std::log(2.0);

        return
            expect_near(
                common::rapidity(
                    common::FourVector{5.0, 7.0, -9.0, 0.0}),
                0.0,
                "zero longitudinal rapidity") &&
            expect_near(
                common::rapidity(
                    common::FourVector{5.0, 7.0, -9.0, 3.0}),
                log_two,
                "positive longitudinal rapidity") &&
            expect_near(
                common::rapidity(
                    common::FourVector{5.0, -8.0, 11.0, -3.0}),
                -log_two,
                "negative longitudinal rapidity");
    }

    /**
     * @brief Verify pseudorapidity calculations.
     *
     * @return true when all pseudorapidity cases pass.
     */
    bool verify_pseudorapidity() noexcept {
        const double log_five = std::log(5.0);

        return
            expect_near(
                common::pseudorapidity(
                    common::FourVector{17.0, 3.0, 4.0, 0.0}),
                0.0,
                "zero pseudorapidity") &&
            expect_near(
                common::pseudorapidity(
                    common::FourVector{17.0, 3.0, 4.0, 12.0}),
                log_five,
                "positive pseudorapidity") &&
            expect_near(
                common::pseudorapidity(
                    common::FourVector{-23.0, 3.0, 4.0, -12.0}),
                -log_five,
                "negative pseudorapidity");
    }


    /**
     * @brief Verify invariant-mass calculations.
     *
     * @return true when direct mass-squared and validated mass cases pass.
     */
    bool verify_invariant_mass() noexcept {
        const common::FourVector positive{5.0, 3.0, 0.0, 0.0};
        const common::FourVector zero{5.0, 3.0, 4.0, 0.0};
        const common::FourVector negative{4.0, 3.0, 4.0, 0.0};

        return
            expect_near(
                common::invariant_mass_squared(positive),
                16.0,
                "positive invariant mass squared") &&
            expect_near(
                common::invariant_mass_squared(zero),
                0.0,
                "zero invariant mass squared") &&
            expect_near(
                common::invariant_mass_squared(negative),
                -9.0,
                "negative invariant mass squared preserved") &&
            expect_near(
                common::invariant_mass(16.0),
                4.0,
                "invariant mass from validated mass squared");
    }

}  // namespace

/**
 * @brief Run the common single-particle kinematics unit tests.
 *
 * @return EXIT_SUCCESS when every test passes, otherwise EXIT_FAILURE.
 */
int main() {
    bool success = true;

    if (!verify_transverse_momentum()) {
        success = false;
    }

    if (!verify_momentum_magnitude()) {
        success = false;
    }

    if (!verify_rapidity()) {
        success = false;
    }

    if (!verify_pseudorapidity()) {
        success = false;
    }

    if (!verify_invariant_mass()) {
        success = false;
    }

    return success ? EXIT_SUCCESS : EXIT_FAILURE;
}
