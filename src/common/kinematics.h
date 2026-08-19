/**
 * @file kinematics.h
 * @brief Common single-particle kinematic operations.
 *
 * This file declares kinematic operations that are independent of a specific
 * analysis module and can therefore be reused throughout the project.
 */

#ifndef SMASH_AFTERBURNER_ANALYSIS_COMMON_KINEMATICS_H
#define SMASH_AFTERBURNER_ANALYSIS_COMMON_KINEMATICS_H

#include "common/four_vector.h"

namespace common {

    /**
     * @brief Calculate the transverse momentum of a four-momentum.
     *
     * The input FourVector is interpreted as a contravariant four-momentum
     *
     *     p^mu = (E, p_x, p_y, p_z),
     *
     * with x0 = E, x1 = p_x, x2 = p_y, and x3 = p_z.
     *
     * The transverse momentum is defined as
     *
     *     p_T = sqrt(p_x^2 + p_y^2).
     *
     * The operation depends only on the spatial transverse components and does
     * not use the energy or longitudinal-momentum components.
     *
     * @param momentum Four-momentum whose transverse momentum is requested.
     * @return Non-negative transverse-momentum magnitude in the same momentum
     *         units as the input components.
     */
    [[nodiscard]] double transverse_momentum(
        const FourVector& momentum) noexcept;

    /**
     * @brief Calculate the magnitude of the three-momentum.
     *
     * The input FourVector is interpreted as a contravariant four-momentum
     *
     *     p^mu = (E, p_x, p_y, p_z).
     *
     * The three-momentum magnitude is defined as
     *
     *     |p| = sqrt(p_x^2 + p_y^2 + p_z^2).
     *
     * The operation depends only on the spatial momentum components and does
     * not use the energy component.
     *
     * @param momentum Four-momentum whose spatial magnitude is requested.
     * @return Non-negative three-momentum magnitude in the same momentum units
     *         as the input spatial components.
     */
    [[nodiscard]] double momentum_magnitude(
        const FourVector& momentum) noexcept;

    /**
     * @brief Calculate the longitudinal rapidity of a four-momentum.
     *
     * The input FourVector is interpreted as a contravariant four-momentum
     *
     *     p^mu = (E, p_x, p_y, p_z).
     *
     * Rapidity is defined as
     *
     *     y = 0.5 * ln((E + p_z) / (E - p_z)).
     *
     * The operation depends only on the energy and longitudinal-momentum
     * components.
     *
     * @param momentum Four-momentum whose rapidity is requested.
     * @return Dimensionless longitudinal rapidity.
     */
    [[nodiscard]] double rapidity(const FourVector& momentum) noexcept;

    /**
     * @brief Calculate the pseudorapidity of a four-momentum.
     *
     * The input FourVector is interpreted as a contravariant four-momentum
     *
     *     p^mu = (E, p_x, p_y, p_z).
     *
     * Pseudorapidity is defined as
     *
     *     eta = 0.5 * ln((|p| + p_z) / (|p| - p_z)).
     *
     * The operation depends only on the spatial momentum components and does
     * not use the energy component.
     *
     * @param momentum Four-momentum whose pseudorapidity is requested.
     * @return Dimensionless pseudorapidity.
     */
    [[nodiscard]] double pseudorapidity(
        const FourVector& momentum) noexcept;

    /**
     * @brief Calculate invariant mass squared from a four-momentum.
     *
     * The input FourVector is interpreted as a contravariant four-momentum
     *
     *     p^mu = (E, p_x, p_y, p_z).
     *
     * The invariant mass squared is calculated directly as
     *
     *     m^2 = E^2 - p_x^2 - p_y^2 - p_z^2.
     *
     * No validation, clamping, tolerance, or fallback is applied. A negative
     * or non-finite result is returned unchanged for the caller to validate.
     *
     * @param momentum Four-momentum whose invariant mass squared is requested.
     * @return Direct invariant-mass-squared result in squared momentum units.
     */
    [[nodiscard]] double invariant_mass_squared(
        const FourVector& momentum) noexcept;

    /**
     * @brief Calculate invariant mass from validated invariant mass squared.
     *
     * This operation calculates
     *
     *     m = sqrt(m^2).
     *
     * The caller is responsible for validating that mass_squared is finite and
     * strictly positive before calling this function. No validation, clamping,
     * tolerance, or fallback is applied here.
     *
     * @param mass_squared Validated positive invariant mass squared.
     * @return Invariant mass in the same momentum units as the source momentum.
     */
    [[nodiscard]] double invariant_mass(double mass_squared) noexcept;

}  // namespace common

#endif  // SMASH_AFTERBURNER_ANALYSIS_COMMON_KINEMATICS_H
