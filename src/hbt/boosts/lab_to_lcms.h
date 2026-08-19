/**
 * @file lab_to_lcms.h
 * @brief HBT Lab-to-LCMS relative-separation transformation.
 */

#ifndef SMASH_AFTERBURNER_ANALYSIS_HBT_BOOSTS_LAB_TO_LCMS_H
#define SMASH_AFTERBURNER_ANALYSIS_HBT_BOOSTS_LAB_TO_LCMS_H

#include "common/four_vector.h"
#include "hbt/pair/pair_kinematics.h"

namespace hbt {

/**
 * @brief Transient result of transforming one pair from Lab to LCMS.
 *
 * beta_out is derived while the Lab-to-LCMS boost parameters are available so
 * the later LCMS-to-PRF transformation does not reconstruct them.
 */
struct LabToLCMSResult {
    /// Relative separation (dt, dx, dy, dz) in LCMS, in fm.
    common::FourVector relative_separation_fm;
    /// Pair speed along the LCMS out direction, in units with c = 1.
    double beta_out;
};

/**
 * @brief Transform one Lab-frame relative separation to the pair LCMS.
 * @param relative_separation_lab_fm Relative separation xa - xb in Lab, fm.
 * @param kinematics Validated reusable Lab-frame pair kinematics.
 * @return LCMS relative separation and the derived PRF boost speed.
 *
 * @pre kinematics belongs to the same ordered physical pair as the supplied
 *      relative separation.
 * @pre pair_four_momentum.x0 is strictly positive and the total pair
 *      four-momentum is future timelike.
 * @pre kt_gev is finite and non-negative.
 *
 * beta_LCMS is Pz / P0. If beta_LCMS is exactly zero, Lab already is the LCMS
 * and the separation is returned unchanged. Any representable non-zero beta,
 * however small, uses the complete Lorentz transformation without an epsilon
 * threshold. No numerical validation, clamp, repair, or fallback is performed.
 */
[[nodiscard]] LabToLCMSResult boost_lab_to_lcms(
    const common::FourVector& relative_separation_lab_fm,
    const PairKinematics& kinematics
) noexcept;

}  // namespace hbt

#endif  // SMASH_AFTERBURNER_ANALYSIS_HBT_BOOSTS_LAB_TO_LCMS_H
