/**
 * @file lab_to_lcms.cpp
 * @brief HBT Lab-to-LCMS relative-separation transformation.
 */

#include "hbt/boosts/lab_to_lcms.h"

#include <cmath>

namespace hbt {

LabToLCMSResult boost_lab_to_lcms(
    const common::FourVector& relative_separation_lab_fm,
    const PairKinematics& kinematics
) noexcept {
    const common::FourVector& pair = kinematics.pair_four_momentum;
    const double beta_lcms = pair.x3 / pair.x0;
    const double pair_transverse_momentum_gev = 2.0 * kinematics.kt_gev;

    if (beta_lcms == 0.0) {
        return {
            relative_separation_lab_fm,
            pair_transverse_momentum_gev / pair.x0
        };
    }

    const double beta_lcms_squared = beta_lcms * beta_lcms;
    const double gamma_lcms = 1.0 / std::sqrt(
        1.0 - beta_lcms_squared
    );

    const common::FourVector relative_separation_lcms_fm{
        gamma_lcms * (
            relative_separation_lab_fm.x0 -
            beta_lcms * relative_separation_lab_fm.x3
        ),
        relative_separation_lab_fm.x1,
        relative_separation_lab_fm.x2,
        gamma_lcms * (
            relative_separation_lab_fm.x3 -
            beta_lcms * relative_separation_lab_fm.x0
        )
    };

    const double beta_out =
        pair_transverse_momentum_gev * gamma_lcms / pair.x0;

    return {relative_separation_lcms_fm, beta_out};
}

}  // namespace hbt
