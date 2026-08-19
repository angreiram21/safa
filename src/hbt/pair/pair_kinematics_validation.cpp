/**
 * @file pair_kinematics_validation.cpp
 * @brief Pure numerical validation for calculated HBT pair kinematics.
 */

#include "hbt/pair/pair_kinematics_validation.h"

#include <cmath>

namespace hbt {

bool is_finite_pair_kt(const PairKinematics& kinematics) noexcept {
    return std::isfinite(kinematics.kt_gev);
}

bool is_finite_pair_mt(const PairKinematics& kinematics) noexcept {
    return std::isfinite(kinematics.mt_gev);
}

}  // namespace hbt
