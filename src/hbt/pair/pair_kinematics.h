/**
 * @file pair_kinematics.h
 * @brief Calculation of reusable pair-level kinematics.
 */

#ifndef SMASH_AFTERBURNER_ANALYSIS_HBT_PAIR_PAIR_KINEMATICS_H
#define SMASH_AFTERBURNER_ANALYSIS_HBT_PAIR_PAIR_KINEMATICS_H

#include "common/four_vector.h"
#include "hbt/event/particle.h"

namespace hbt {

/**
 * @brief Reusable kinematics calculated once for one physical pair.
 *
 * Values are calculated from final Afterburner momenta and the validated
 * invariant masses already stored in the two accepted particles. The total
 * four-momentum is retained so later HBT frame transformations can reuse the
 * same pair sums without reconstructing them.
 */
struct PairKinematics {
    /// Total Lab-frame pair four-momentum (E, px, py, pz) in GeV.
    common::FourVector pair_four_momentum;
    double kx_gev;  ///< Pair-average Kx = (px_a + px_b) / 2 in GeV.
    double ky_gev;  ///< Pair-average Ky = (py_a + py_b) / 2 in GeV.
    double kt_gev;  ///< Pair transverse momentum kT in GeV.
    double mt_gev;  ///< Pair transverse mass mT in GeV.
};

/**
 * @brief Calculate reusable kinematics for one already formed particle pair.
 * @param particle_a First accepted particle in the pair.
 * @param particle_b Second accepted particle in the pair.
 * @return Pair kinematics calculated from stored final particle data.
 *
 * The total Lab-frame pair four-momentum is formed once. Kx and Ky are then
 * derived from its transverse components, followed by
 *
 *     kT = sqrt(Kx^2 + Ky^2)
 *     m_avg = (m_a + m_b) / 2
 *     mT = sqrt(m_avg^2 + kT^2)
 *
 * Particle masses are reused from Particle::invariant_mass_gev and are not
 * reconstructed from the four-momenta. This function performs calculation
 * only: it does not validate, reject, clamp, route, slice, boost, or report.
 */
[[nodiscard]] PairKinematics calculate_pair_kinematics(
    const Particle& particle_a,
    const Particle& particle_b
) noexcept;

}  // namespace hbt

#endif  // SMASH_AFTERBURNER_ANALYSIS_HBT_PAIR_PAIR_KINEMATICS_H
