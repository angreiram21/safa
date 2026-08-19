/**
 * @file pair_frame_observables_validation.h
 * @brief Numerical validation for calculated HBT pair-frame observables.
 */

#ifndef HBT_PAIR_PAIR_FRAME_OBSERVABLES_VALIDATION_H
#define HBT_PAIR_PAIR_FRAME_OBSERVABLES_VALIDATION_H

#include "hbt/pair/pair_frame_observables.h"

#include <optional>

namespace hbt {

/**
 * @brief Exact non-finite reason found in calculated frame observables.
 */
enum class PairFrameObservableNumericalReason {
    NonFiniteDeltaTLab,    ///< delta_t_lab is not finite.
    NonFiniteDeltaTLcms,   ///< delta_t_lcms is not finite.
    NonFiniteDeltaTPrf,    ///< delta_t_prf is not finite.
    NonFiniteROutLcms,     ///< r_out_lcms is not finite.
    NonFiniteROutPrf,      ///< r_out_prf is not finite.
    NonFiniteRSide,        ///< r_side is not finite.
    NonFiniteRLong,        ///< r_long is not finite.
    NonFiniteRRadialLcms,  ///< r_radial_lcms is not finite.
    NonFiniteRRadialPrf    ///< r_radial_prf is not finite.
};

/**
 * @brief Find the first non-finite calculated frame observable.
 * @param observables Already calculated pair-frame observables.
 * @return Exact first non-finite reason, or std::nullopt when all are finite.
 *
 * The check order is stable and follows the stored observable fields. This
 * function validates only; it does not reject, route, report, clamp, or modify
 * the supplied observables.
 */
[[nodiscard]] std::optional<PairFrameObservableNumericalReason>
first_non_finite_pair_frame_observable(
    const PairFrameObservables& observables
) noexcept;

}  // namespace hbt

#endif  // HBT_PAIR_PAIR_FRAME_OBSERVABLES_VALIDATION_H
