/**
 * @file four_vector.h
 * @brief Minimal representation of a relativistic four-vector.
 *
 * This file defines the FourVector type used to store four-vector components
 * throughout the project.
 *
 * FourVector is a mathematical data container. It does not assign a specific
 * physical interpretation, unit system, metric convention, or reference frame
 * to its components.
 *
 * The stored components correspond conceptually to the contravariant
 * components
 *
 *     x^mu = (x^0, x^1, x^2, x^3).
 *
 * Depending on context, a FourVector may represent quantities such as a
 * spacetime position or a four-momentum. The physical meaning and units of
 * the components must be documented by the owning type or the operation that
 * uses the vector.
 */

#ifndef SMASH_AFTERBURNER_ANALYSIS_COMMON_FOUR_VECTOR_H
#define SMASH_AFTERBURNER_ANALYSIS_COMMON_FOUR_VECTOR_H

namespace common {

    /**
     * @brief Four-component mathematical representation of a contravariant
     *        relativistic four-vector.
     *
     * The members x0, x1, x2, and x3 identify component numbers 0, 1, 2, and 3,
     * respectively. The member names do not denote lowered tensor indices:
     * x0 represents x^0, not x_0, and analogously for the remaining components.
     *
     * No metric signature is assumed by this type. Operations whose result
     * depends on a spacetime metric must define that convention explicitly in
     * the module that implements them.
     *
     * No units are imposed by this type. Units are determined by the physical
     * quantity represented by a particular FourVector instance.
     */
    struct FourVector {
        double x0{0.0};  ///< Contravariant component x^0.
        double x1{0.0};  ///< Contravariant component x^1.
        double x2{0.0};  ///< Contravariant component x^2.
        double x3{0.0};  ///< Contravariant component x^3.
    };

}  // namespace common

#endif  // SMASH_AFTERBURNER_ANALYSIS_COMMON_FOUR_VECTOR_H
