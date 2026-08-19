/**
 * @file pair_frame_observables_validation.cpp
 * @brief Numerical validation for calculated HBT pair-frame observables.
 */

#include "hbt/pair/pair_frame_observables_validation.h"

#include <cmath>

namespace hbt {

std::optional<PairFrameObservableNumericalReason>
first_non_finite_pair_frame_observable(
    const PairFrameObservables& observables
) noexcept {
    if (!std::isfinite(observables.delta_t_lab_fm)) {
        return PairFrameObservableNumericalReason::NonFiniteDeltaTLab;
    }
    if (!std::isfinite(observables.delta_t_lcms_fm)) {
        return PairFrameObservableNumericalReason::NonFiniteDeltaTLcms;
    }
    if (!std::isfinite(observables.delta_t_prf_fm)) {
        return PairFrameObservableNumericalReason::NonFiniteDeltaTPrf;
    }
    if (!std::isfinite(observables.r_out_lcms_fm)) {
        return PairFrameObservableNumericalReason::NonFiniteROutLcms;
    }
    if (!std::isfinite(observables.r_out_prf_fm)) {
        return PairFrameObservableNumericalReason::NonFiniteROutPrf;
    }
    if (!std::isfinite(observables.r_side_fm)) {
        return PairFrameObservableNumericalReason::NonFiniteRSide;
    }
    if (!std::isfinite(observables.r_long_fm)) {
        return PairFrameObservableNumericalReason::NonFiniteRLong;
    }
    if (!std::isfinite(observables.r_radial_lcms_fm)) {
        return PairFrameObservableNumericalReason::NonFiniteRRadialLcms;
    }
    if (!std::isfinite(observables.r_radial_prf_fm)) {
        return PairFrameObservableNumericalReason::NonFiniteRRadialPrf;
    }
    return std::nullopt;
}

}  // namespace hbt
