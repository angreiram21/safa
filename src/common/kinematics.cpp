/**
 * @file kinematics.cpp
 * @brief Implementation of common single-particle kinematic operations.
 */

#include "common/kinematics.h"

#include <cmath>

namespace common {

    double transverse_momentum(const FourVector& momentum) noexcept {
        return std::hypot(momentum.x1, momentum.x2);
    }

    double momentum_magnitude(const FourVector& momentum) noexcept {
        return std::hypot(momentum.x1, momentum.x2, momentum.x3);
    }

    double rapidity(const FourVector& momentum) noexcept {
        return 0.5 * std::log(
            (momentum.x0 + momentum.x3) /
            (momentum.x0 - momentum.x3));
    }

    double pseudorapidity(const FourVector& momentum) noexcept {
        const double magnitude = momentum_magnitude(momentum);

        return 0.5 * std::log(
            (magnitude + momentum.x3) /
            (magnitude - momentum.x3));
    }

    double invariant_mass_squared(const FourVector& momentum) noexcept {
        return momentum.x0 * momentum.x0 -
               momentum.x1 * momentum.x1 -
               momentum.x2 * momentum.x2 -
               momentum.x3 * momentum.x3;
    }

    double invariant_mass(double mass_squared) noexcept {
        return std::sqrt(mass_squared);
    }

}  // namespace common
