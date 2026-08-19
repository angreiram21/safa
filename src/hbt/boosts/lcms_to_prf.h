/**
 * @file lcms_to_prf.h
 * @brief HBT LCMS-to-PRF boost along the OSL out direction.
 */

#ifndef SMASH_AFTERBURNER_ANALYSIS_HBT_BOOSTS_LCMS_TO_PRF_H
#define SMASH_AFTERBURNER_ANALYSIS_HBT_BOOSTS_LCMS_TO_PRF_H

namespace hbt {

/**
 * @brief Transformed temporal and out components in the pair rest frame.
 */
struct LCMSToPRFResult {
    double delta_t_prf_fm;  ///< Relative emission time in PRF, in fm.
    double r_out_prf_fm;    ///< Out separation in PRF, in fm.
};

/**
 * @brief Boost LCMS temporal/out separation components to the PRF.
 * @param delta_t_lcms_fm Relative emission time in LCMS, in fm.
 * @param r_out_lcms_fm Relative out separation in LCMS, in fm.
 * @param beta_out Pair speed along the LCMS out direction.
 * @return Transformed PRF temporal and out components.
 *
 * @pre beta_out is finite and satisfies 0 <= beta_out < 1.
 *
 * If beta_out is exactly zero, LCMS already is the PRF and the two values are
 * returned unchanged. Any representable non-zero beta uses the complete
 * Lorentz transformation without an epsilon threshold. Side and long do not
 * enter this function because this one-dimensional boost leaves them unchanged.
 */
[[nodiscard]] LCMSToPRFResult boost_lcms_to_prf(
    double delta_t_lcms_fm,
    double r_out_lcms_fm,
    double beta_out
) noexcept;

}  // namespace hbt

#endif  // SMASH_AFTERBURNER_ANALYSIS_HBT_BOOSTS_LCMS_TO_PRF_H
