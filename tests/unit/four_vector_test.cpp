/**
 * @file four_vector_test.cpp
 * @brief Unit tests for the minimal FourVector data contract.
 *
 * This test verifies the behavior guaranteed by four_vector.h:
 * default initialization to zero, aggregate initialization with explicit
 * component values, and independent storage of the four components.
 */

#include "common/four_vector.h"

#include <cstdlib>
#include <iostream>

namespace {

    /**
     * @brief Verify that a default-constructed FourVector contains four zeros.
     *
     * @return true when all four components are zero, otherwise false.
     */
    bool verify_default_initialization() noexcept {
        const common::FourVector vector{};

        return vector.x0 == 0.0 &&
               vector.x1 == 0.0 &&
               vector.x2 == 0.0 &&
               vector.x3 == 0.0;
    }

    /**
     * @brief Verify aggregate initialization with four explicit component
     *        values.
     *
     * @return true when each value is stored in the expected component,
     *         otherwise false.
     */
    bool verify_explicit_initialization() noexcept {
        const common::FourVector vector{
            1.0,
            -2.0,
            3.5,
            -4.25
        };

        return vector.x0 == 1.0 &&
               vector.x1 == -2.0 &&
               vector.x2 == 3.5 &&
               vector.x3 == -4.25;
    }

    /**
     * @brief Verify that the four components can be modified independently.
     *
     * Each component is assigned separately while the values stored in the
     * remaining components are preserved.
     *
     * @return true when all four components retain their independently assigned
     *         values, otherwise false.
     */
    bool verify_independent_components() noexcept {
        common::FourVector vector{};

        vector.x0 = 10.0;
        vector.x1 = 20.0;
        vector.x2 = 30.0;
        vector.x3 = 40.0;

        return vector.x0 == 10.0 &&
               vector.x1 == 20.0 &&
               vector.x2 == 30.0 &&
               vector.x3 == 40.0;
    }

}  // namespace

/**
 * @brief Run the FourVector unit tests.
 *
 * @return EXIT_SUCCESS when every test passes, otherwise EXIT_FAILURE.
 */
int main() {
    if (!verify_default_initialization()) {
        std::cerr << "four_vector_test: default initialization failed.\n";
        return EXIT_FAILURE;
    }

    if (!verify_explicit_initialization()) {
        std::cerr << "four_vector_test: explicit initialization failed.\n";
        return EXIT_FAILURE;
    }

    if (!verify_independent_components()) {
        std::cerr
            << "four_vector_test: independent component storage failed.\n";
        return EXIT_FAILURE;
    }

    return EXIT_SUCCESS;
}
