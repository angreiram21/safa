/**
 * @file lcms_to_prf.cpp
 * @brief HBT LCMS-to-PRF boost along the OSL out direction.
 */

#include "hbt/boosts/lcms_to_prf.h"

#include <cmath>

namespace hbt {

LCMSToPRFResult boost_lcms_to_prf(
    double delta_t_lcms_fm,
    double r_out_lcms_fm,
    double beta_out
) noexcept {
    if (beta_out == 0.0) {
        return {delta_t_lcms_fm, r_out_lcms_fm};
    }

    const double beta_out_squared = beta_out * beta_out;
    const double gamma_out = 1.0 / std::sqrt(
        1.0 - beta_out_squared
    );

    return {
        gamma_out * (
            delta_t_lcms_fm - beta_out * r_out_lcms_fm
        ),
        gamma_out * (
            r_out_lcms_fm - beta_out * delta_t_lcms_fm
        )
    };
}

}  // namespace hbt
