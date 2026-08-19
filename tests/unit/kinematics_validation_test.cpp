/**
 * @file kinematics_validation_test.cpp
 * @brief Unit tests for non-fatal common kinematic validation.
 */

#include "common/kinematics_validation.h"

#include <cstdlib>
#include <iostream>
#include <limits>

namespace {

    /**
     * @brief Verify finite four-momentum validation.
     * @return true when finite input passes and every non-finite case fails.
     */
    bool verify_finite_four_momentum() noexcept {
        const double nan = std::numeric_limits<double>::quiet_NaN();
        const double inf = std::numeric_limits<double>::infinity();

        if (!common::is_finite_four_momentum({5.0, 1.0, 2.0, 3.0})) {
            std::cerr
                << "kinematics_validation_test: finite momentum rejected.\n";
            return false;
        }

        const common::FourVector invalid[] = {
            {nan, 1.0, 2.0, 3.0},
            {inf, 1.0, 2.0, 3.0},
            {-inf, 1.0, 2.0, 3.0},
            {5.0, nan, 2.0, 3.0},
            {5.0, inf, 2.0, 3.0},
            {5.0, -inf, 2.0, 3.0},
            {5.0, 1.0, nan, 3.0},
            {5.0, 1.0, inf, 3.0},
            {5.0, 1.0, -inf, 3.0},
            {5.0, 1.0, 2.0, nan},
            {5.0, 1.0, 2.0, inf},
            {5.0, 1.0, 2.0, -inf}
        };

        for (const common::FourVector& momentum : invalid) {
            if (common::is_finite_four_momentum(momentum)) {
                std::cerr
                    << "kinematics_validation_test: non-finite momentum "
                    << "accepted.\n";
                return false;
            }
        }

        return true;
    }

    /**
     * @brief Verify rapidity input-domain validation.
     * @return true when only finite E > 0 and |p_z| < E inputs pass.
     */
    bool verify_rapidity_input() noexcept {
        const double nan = std::numeric_limits<double>::quiet_NaN();

        if (!common::is_valid_rapidity_input({5.0, 1.0, 2.0, 3.0}) ||
            !common::is_valid_rapidity_input({5.0, 1.0, 2.0, -3.0})) {
            std::cerr
                << "kinematics_validation_test: valid rapidity input "
                << "rejected.\n";
            return false;
        }

        const common::FourVector invalid[] = {
            {0.0, 1.0, 2.0, 0.0},
            {-5.0, 1.0, 2.0, 0.0},
            {5.0, 1.0, 2.0, 5.0},
            {5.0, 1.0, 2.0, -5.0},
            {5.0, 1.0, 2.0, 6.0},
            {5.0, 1.0, 2.0, -6.0},
            {nan, 1.0, 2.0, 0.0}
        };

        for (const common::FourVector& momentum : invalid) {
            if (common::is_valid_rapidity_input(momentum)) {
                std::cerr
                    << "kinematics_validation_test: invalid rapidity input "
                    << "accepted.\n";
                return false;
            }
        }

        return true;
    }

    /**
     * @brief Verify pseudorapidity input-domain validation.
     * @return true when finite non-zero transverse momentum is required.
     */
    bool verify_pseudorapidity_input() noexcept {
        const double inf = std::numeric_limits<double>::infinity();

        if (!common::is_valid_pseudorapidity_input(
                {5.0, 1.0, 0.0, 3.0}) ||
            !common::is_valid_pseudorapidity_input(
                {5.0, 0.0, -1.0, -3.0})) {
            std::cerr
                << "kinematics_validation_test: valid pseudorapidity input "
                << "rejected.\n";
            return false;
        }

        const common::FourVector invalid[] = {
            {5.0, 0.0, 0.0, 0.0},
            {5.0, 0.0, 0.0, 3.0},
            {5.0, 0.0, 0.0, -3.0},
            {5.0, inf, 0.0, 3.0}
        };

        for (const common::FourVector& momentum : invalid) {
            if (common::is_valid_pseudorapidity_input(momentum)) {
                std::cerr
                    << "kinematics_validation_test: invalid pseudorapidity "
                    << "input accepted.\n";
                return false;
            }
        }

        return true;
    }

    /**
     * @brief Verify generic finite-result validation.
     * @return true when finite values pass and NaN or infinities fail.
     */
    bool verify_finite_result() noexcept {
        const double nan = std::numeric_limits<double>::quiet_NaN();
        const double inf = std::numeric_limits<double>::infinity();

        if (!common::is_finite_kinematic_result(0.0) ||
            !common::is_finite_kinematic_result(
                std::numeric_limits<double>::max())) {
            std::cerr
                << "kinematics_validation_test: finite result rejected.\n";
            return false;
        }

        if (common::is_finite_kinematic_result(nan) ||
            common::is_finite_kinematic_result(inf) ||
            common::is_finite_kinematic_result(-inf)) {
            std::cerr
                << "kinematics_validation_test: non-finite result accepted.\n";
            return false;
        }

        return true;
    }

    /**
     * @brief Verify invariant-mass-squared validation.
     * @return true when only finite strictly positive values pass.
     */
    bool verify_invariant_mass_squared() noexcept {
        const double nan = std::numeric_limits<double>::quiet_NaN();
        const double inf = std::numeric_limits<double>::infinity();

        if (!common::is_valid_invariant_mass_squared(0.019) ||
            !common::is_valid_invariant_mass_squared(1.0)) {
            std::cerr
                << "kinematics_validation_test: valid mass squared "
                << "rejected.\n";
            return false;
        }

        const double invalid[] = {0.0, -0.001, nan, inf, -inf};

        for (const double value : invalid) {
            if (common::is_valid_invariant_mass_squared(value)) {
                std::cerr
                    << "kinematics_validation_test: invalid mass squared "
                    << "accepted.\n";
                return false;
            }
        }

        return true;
    }

}  // namespace

/**
 * @brief Run the complete non-fatal kinematic-validation unit tests.
 * @return EXIT_SUCCESS when every test passes, otherwise EXIT_FAILURE.
 */
int main() {
    bool success = true;
    success = verify_finite_four_momentum() && success;
    success = verify_rapidity_input() && success;
    success = verify_pseudorapidity_input() && success;
    success = verify_finite_result() && success;
    success = verify_invariant_mass_squared() && success;
    return success ? EXIT_SUCCESS : EXIT_FAILURE;
}
