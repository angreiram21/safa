/**
 * @file pair_kinematics.cpp
 * @brief Pair-level kinematics implementation.
 */

#include "hbt/pair/pair_kinematics.h"

#include <cmath>

namespace hbt {

PairKinematics calculate_pair_kinematics(
    const Particle& particle_a,
    const Particle& particle_b
) noexcept {
    const common::FourVector pair_four_momentum{
        particle_a.momentum.x0 + particle_b.momentum.x0,
        particle_a.momentum.x1 + particle_b.momentum.x1,
        particle_a.momentum.x2 + particle_b.momentum.x2,
        particle_a.momentum.x3 + particle_b.momentum.x3
    };

    const double kx_gev = 0.5 * pair_four_momentum.x1;
    const double ky_gev = 0.5 * pair_four_momentum.x2;
    const double kt_gev = std::hypot(kx_gev, ky_gev);

    const double average_mass_gev = 0.5 * (
        particle_a.invariant_mass_gev + particle_b.invariant_mass_gev
    );
    const double mt_gev = std::hypot(average_mass_gev, kt_gev);

    return {
        pair_four_momentum,
        kx_gev,
        ky_gev,
        kt_gev,
        mt_gev
    };
}

}  // namespace hbt
