/**
 * @file pair_frame_observables.h
 * @brief HBT pair separation observables in Lab, LCMS, and PRF.
 */

#ifndef SMASH_AFTERBURNER_ANALYSIS_HBT_PAIR_PAIR_FRAME_OBSERVABLES_H
#define SMASH_AFTERBURNER_ANALYSIS_HBT_PAIR_PAIR_FRAME_OBSERVABLES_H

#include "hbt/event/particle.h"
#include "hbt/pair/pair_kinematics.h"

namespace hbt {

/**
 * @brief Reusable frame observables calculated once for one physical pair.
 *
 * The three relative emission times are retained separately for frame-to-frame
 * temporal comparisons. Side and long are stored once because the subsequent
 * LCMS-to-PRF boost is purely along out and leaves those components unchanged.
 */
struct PairFrameObservables {
    double delta_t_lab_fm;       ///< Relative emission time in Lab, in fm.
    double delta_t_lcms_fm;      ///< Relative emission time in LCMS, in fm.
    double delta_t_prf_fm;       ///< Relative emission time in PRF, in fm.
    double r_out_lcms_fm;        ///< OSL out separation in LCMS, in fm.
    double r_out_prf_fm;         ///< OSL out separation in PRF, in fm.
    double r_side_fm;            ///< Shared LCMS/PRF side separation, in fm.
    double r_long_fm;            ///< Shared LCMS/PRF long separation, in fm.
    double r_radial_lcms_fm;     ///< Spatial LCMS radial separation, in fm.
    double r_radial_prf_fm;      ///< Spatial PRF radial separation, in fm.
};

/**
 * @brief Calculate all Phase-5 pair frame observables exactly once.
 * @param particle_a First accepted particle in canonical pair order.
 * @param particle_b Second accepted particle in canonical pair order.
 * @param kinematics Validated pair kinematics for the same ordered pair.
 * @return Lab/LCMS/PRF temporal and spatial separation observables.
 *
 * @pre particle_a and particle_b are accepted HBT particles from one subevent.
 * @pre Their individual transverse momenta are strictly positive under the
 *      already applied particle-acceptance contract.
 * @pre kinematics is the already validated result for this exact ordered pair.
 *
 * The relative separation is xa - xb. Lab-to-LCMS is performed once on that
 * relative four-vector rather than on two individual positions. For kT > 0,
 * the OSL out axis is K_T / kT. For exact kT == 0, qT is calculated lazily from
 * pa_T - pb_T and used as the out direction; the accepted-particle pT contract
 * guarantees qT > 0 in that branch. No legacy epsilon or second fallback is
 * used. LCMS-to-PRF then boosts only the temporal/out pair of components.
 */
[[nodiscard]] PairFrameObservables calculate_pair_frame_observables(
    const Particle& particle_a,
    const Particle& particle_b,
    const PairKinematics& kinematics
) noexcept;

}  // namespace hbt

#endif  // SMASH_AFTERBURNER_ANALYSIS_HBT_PAIR_PAIR_FRAME_OBSERVABLES_H
