/**
 * @file kinematics_validation.cpp
 * @brief Implementation of non-fatal single-particle kinematic validation.
 */

#include "common/kinematics_validation.h"

#include <cmath>

namespace common {

    bool is_finite_four_momentum(const FourVector& momentum) noexcept {
        return std::isfinite(momentum.x0) &&
               std::isfinite(momentum.x1) &&
               std::isfinite(momentum.x2) &&
               std::isfinite(momentum.x3);
    }

    bool is_valid_rapidity_input(const FourVector& momentum) noexcept {
        if (!is_finite_four_momentum(momentum)) {
            return false;
        }

        return momentum.x0 > 0.0 && std::abs(momentum.x3) < momentum.x0;
    }

    bool is_valid_pseudorapidity_input(
        const FourVector& momentum) noexcept {
        if (!is_finite_four_momentum(momentum)) {
            return false;
        }

        return momentum.x1 != 0.0 || momentum.x2 != 0.0;
    }

    bool is_finite_kinematic_result(double result) noexcept {
        return std::isfinite(result);
    }

    bool is_valid_invariant_mass_squared(double mass_squared) noexcept {
        return std::isfinite(mass_squared) && mass_squared > 0.0;
    }

}  // namespace common
