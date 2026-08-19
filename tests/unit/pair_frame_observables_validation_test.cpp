/**
 * @file pair_frame_observables_validation_test.cpp
 * @brief Unit tests for pair-frame observable finiteness validation.
 */

#include "hbt/pair/pair_frame_observables_validation.h"

#include <array>
#include <cstdlib>
#include <iostream>
#include <limits>
#include <optional>

namespace {

/**
 * @brief Return one completely finite frame-observable object.
 * @return Finite representative observables.
 */
hbt::PairFrameObservables finite_observables() {
    return {1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0, 8.0, 9.0};
}

/**
 * @brief Verify all finite observables pass without a numerical reason.
 * @return true when validation returns std::nullopt.
 */
bool verify_all_finite() {
    return !hbt::first_non_finite_pair_frame_observable(
        finite_observables()
    ).has_value();
}

/**
 * @brief Verify every stored field maps to its exact rejection reason.
 * @return true when all nine non-finite fields are identified exactly.
 */
bool verify_exact_non_finite_reasons() {
    using Reason = hbt::PairFrameObservableNumericalReason;
    using Member = double hbt::PairFrameObservables::*;
    const std::array<Member, 9U> members{
        &hbt::PairFrameObservables::delta_t_lab_fm,
        &hbt::PairFrameObservables::delta_t_lcms_fm,
        &hbt::PairFrameObservables::delta_t_prf_fm,
        &hbt::PairFrameObservables::r_out_lcms_fm,
        &hbt::PairFrameObservables::r_out_prf_fm,
        &hbt::PairFrameObservables::r_side_fm,
        &hbt::PairFrameObservables::r_long_fm,
        &hbt::PairFrameObservables::r_radial_lcms_fm,
        &hbt::PairFrameObservables::r_radial_prf_fm
    };
    const std::array<Reason, 9U> reasons{
        Reason::NonFiniteDeltaTLab,
        Reason::NonFiniteDeltaTLcms,
        Reason::NonFiniteDeltaTPrf,
        Reason::NonFiniteROutLcms,
        Reason::NonFiniteROutPrf,
        Reason::NonFiniteRSide,
        Reason::NonFiniteRLong,
        Reason::NonFiniteRRadialLcms,
        Reason::NonFiniteRRadialPrf
    };

    for (std::size_t index = 0U; index < members.size(); ++index) {
        hbt::PairFrameObservables observables = finite_observables();
        observables.*members[index] =
            std::numeric_limits<double>::infinity();
        const std::optional<Reason> reason =
            hbt::first_non_finite_pair_frame_observable(observables);
        if (!reason.has_value() || reason.value() != reasons[index]) {
            std::cerr
                << "pair_frame_observables_validation_test: wrong reason at "
                << index << ".\n";
            return false;
        }
    }
    return true;
}

}  // namespace

/**
 * @brief Run pair-frame observable validation unit tests.
 * @return EXIT_SUCCESS when every test passes, otherwise EXIT_FAILURE.
 */
int main() {
    bool success = true;
    success = verify_all_finite() && success;
    success = verify_exact_non_finite_reasons() && success;
    return success ? EXIT_SUCCESS : EXIT_FAILURE;
}
