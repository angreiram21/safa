/**
 * @file pair_kinematics_validation.h
 * @brief Pure numerical validation for calculated HBT pair kinematics.
 */

#ifndef SMASH_AFTERBURNER_ANALYSIS_HBT_PAIR_PAIR_KINEMATICS_VALIDATION_H
#define SMASH_AFTERBURNER_ANALYSIS_HBT_PAIR_PAIR_KINEMATICS_VALIDATION_H

#include "hbt/pair/pair_kinematics.h"

namespace hbt {

/**
 * @brief Check whether the calculated pair kT is finite.
 * @param kinematics Already calculated pair kinematics.
 * @return true when kT is finite, otherwise false.
 *
 * This function performs validation only. It does not reject, clamp, throw,
 * report, route, or modify the supplied pair kinematics.
 */
[[nodiscard]] bool is_finite_pair_kt(
    const PairKinematics& kinematics
) noexcept;

/**
 * @brief Check whether the calculated pair mT is finite.
 * @param kinematics Already calculated pair kinematics.
 * @return true when mT is finite, otherwise false.
 *
 * This function performs validation only. It does not reject, clamp, throw,
 * report, route, or modify the supplied pair kinematics.
 */
[[nodiscard]] bool is_finite_pair_mt(
    const PairKinematics& kinematics
) noexcept;

}  // namespace hbt

#endif  // SMASH_AFTERBURNER_ANALYSIS_HBT_PAIR_PAIR_KINEMATICS_VALIDATION_H
