/**
 * @file kinematics_validation.h
 * @brief Non-fatal validation of common single-particle kinematic data.
 *
 * This file declares reusable validation operations for kinematic inputs and
 * calculated results. Validation functions report validity only. They do not
 * abort, throw, print diagnostics, repair values, or apply analysis cuts.
 */

#ifndef SMASH_AFTERBURNER_ANALYSIS_COMMON_KINEMATICS_VALIDATION_H
#define SMASH_AFTERBURNER_ANALYSIS_COMMON_KINEMATICS_VALIDATION_H

#include "common/four_vector.h"

namespace common {

    /**
     * @brief Test whether every component of a four-momentum is finite.
     *
     * @param momentum Four-momentum to validate.
     * @return true only when all four components satisfy std::isfinite().
     */
    [[nodiscard]] bool is_finite_four_momentum(
        const FourVector& momentum) noexcept;

    /**
     * @brief Test whether a four-momentum has a valid rapidity input domain.
     *
     * The complete four-momentum must be finite. The energy and longitudinal
     * momentum must additionally satisfy
     *
     *     E > 0
     *     |p_z| < E.
     *
     * @param momentum Four-momentum to validate for rapidity calculation.
     * @return true only when the direct rapidity formula has a valid finite
     *         input domain.
     */
    [[nodiscard]] bool is_valid_rapidity_input(
        const FourVector& momentum) noexcept;

    /**
     * @brief Test whether a four-momentum has a valid pseudorapidity domain.
     *
     * The complete four-momentum must be finite and its transverse momentum
     * must be non-zero:
     *
     *     p_x != 0 or p_y != 0.
     *
     * @param momentum Four-momentum to validate for pseudorapidity calculation.
     * @return true only when the direct pseudorapidity formula has a valid
     *         finite input domain.
     */
    [[nodiscard]] bool is_valid_pseudorapidity_input(
        const FourVector& momentum) noexcept;

    /**
     * @brief Test whether a calculated kinematic result is finite.
     *
     * @param result Calculated kinematic value to validate.
     * @return true only when result satisfies std::isfinite().
     */
    [[nodiscard]] bool is_finite_kinematic_result(double result) noexcept;

    /**
     * @brief Test whether invariant mass squared is valid for sqrt evaluation.
     *
     * The value must be finite and strictly positive. Zero, negative values,
     * NaN, and infinities are invalid. No clamp or tolerance is applied.
     *
     * @param mass_squared Calculated invariant mass squared to validate.
     * @return true only when mass_squared is finite and strictly positive.
     */
    [[nodiscard]] bool is_valid_invariant_mass_squared(
        double mass_squared) noexcept;

}  // namespace common

#endif  // SMASH_AFTERBURNER_ANALYSIS_COMMON_KINEMATICS_VALIDATION_H
