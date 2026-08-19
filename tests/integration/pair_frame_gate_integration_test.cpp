/**
 * @file pair_frame_gate_integration_test.cpp
 * @brief Integration coverage for the frame-observable numerical gate.
 */

#include "hbt/pair/pair_frame_observables.h"
#include "hbt/pair/pair_processor.h"

#include <array>
#include <cstddef>
#include <cstdlib>
#include <iostream>
#include <limits>
#include <optional>
#include <vector>

namespace {

/// Number of calls observed through the frame-calculation test double.
std::size_t frame_observable_call_count = 0U;
/// Selected frame field made non-finite by the calculation test double.
std::optional<std::size_t> non_finite_frame_index;

using FrameObservableMember = double hbt::PairFrameObservables::*;

/// Stable member order matching PairFrameObservableNumericalReason.
constexpr std::array<FrameObservableMember, 9U> frame_observable_members{
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

/**
 * @brief Recording consumer that retains only safe copied routing state.
 */
class RecordingPairFrameConsumer final : public hbt::PairFrameConsumer {
public:
    /**
     * @brief Record one synchronous pair delivery without retaining references.
     * @param context Current already resolved routing context.
     * @param kinematics Existing pair kinematics.
     * @param observables Existing finite frame observables.
     */
    void consume(
        const hbt::PairFrameRouteContext& context,
        const hbt::PairKinematics& kinematics,
        const hbt::PairFrameObservables& observables
    ) override {
        static_cast<void>(kinematics);
        static_cast<void>(observables);
        ++call_count;
        channel_index = context.channel_index;
        channel = context.channel;
        routes = context.origin_routes;
        if (context.pair_slice_route != nullptr) {
            slice_route = *context.pair_slice_route;
        } else {
            slice_route.reset();
        }
    }

    /** @brief Reset all copied state before one isolated scenario. */
    void reset() noexcept {
        call_count = 0U;
        channel_index = 0U;
        channel = hbt::PrimitiveChannelId::PiPlusPiPlus;
        routes = {false, false, false};
        slice_route.reset();
    }

    std::size_t call_count = 0U;  ///< Number of consume calls observed.
    std::size_t channel_index = 0U;  ///< Copied required-channel index.
    hbt::PrimitiveChannelId channel =
        hbt::PrimitiveChannelId::PiPlusPiPlus;  ///< Copied channel identity.
    hbt::PairOriginRoutes routes{false, false, false};  ///< Copied routes.
    std::optional<hbt::PairSliceRoute> slice_route;  ///< Copied slice route.
};

/// Consumer instance reused and reset between isolated scenarios.
RecordingPairFrameConsumer frame_consumer;

/**
 * @brief Reset all frame-gate test-double state.
 */
void reset_frame_state() noexcept {
    frame_observable_call_count = 0U;
    non_finite_frame_index.reset();
    frame_consumer.reset();
}

/**
 * @brief Build one synthetic accepted particle for gate integration tests.
 * @param species Canonical particle species.
 * @param px First transverse momentum component in GeV.
 * @param py Second transverse momentum component in GeV.
 * @param mass Stored invariant mass in GeV.
 * @param pdg Raw signed PDG code retained for diagnostics.
 * @param charge Raw electric charge retained for diagnostics.
 * @return Complete synthetic Particle with primordial inclusive origin flags.
 */
hbt::Particle make_particle(
    hbt::SpeciesId species,
    double px,
    double py,
    double mass,
    int pdg,
    int charge
) {
    return {
        species,
        {0.0, 1.0, 2.0, 3.0},
        {2.0, px, py, 0.0},
        mass,
        {true, true, true},
        pdg,
        charge
    };
}

/**
 * @brief Return one kT slice covering [0, 0.3) GeV.
 * @return Valid kT-only slicing configuration.
 */
hbt::PairSlicingConfig active_slicing() {
    return {{true, {0.0, 0.3}}, {false, {}}};
}

/**
 * @brief Return a configuration with both pair-slicing axes disabled.
 * @return Valid disabled-slicing configuration.
 */
hbt::PairSlicingConfig disabled_slicing() {
    return {{false, {}}, {false, {}}};
}

/**
 * @brief Return one test primitive channel as required-channel input.
 * @param channel Primitive channel to process.
 * @return One-element required-channel vector.
 */
std::vector<hbt::PrimitiveChannelId> one_channel(
    hbt::PrimitiveChannelId channel
) {
    return {channel};
}

/**
 * @brief Report one failed frame-gate integration condition.
 * @param message Description of the violated contract.
 * @return Always false.
 */
bool fail(const char* message) {
    std::cerr << "pair_frame_gate_integration_test: " << message << ".\n";
    return false;
}

/**
 * @brief Verify a kinematically rejected pair never calculates frames.
 * @return true when rejection occurs with zero frame and consumer calls.
 */
bool verify_kinematic_rejection_skips_frames() {
    hbt::EventBuffers buffers;
    const double maximum = std::numeric_limits<double>::max();
    buffers.add(make_particle(
        hbt::SpeciesId::PiPlus, maximum, 0.0, 0.14, 211, 1));
    buffers.add(make_particle(
        hbt::SpeciesId::PiPlus, maximum, 0.0, 0.14, 211, 1));

    reset_frame_state();
    const hbt::PairSubeventProcessingResult result =
        hbt::process_subevent_pairs(
            1U,
            0,
            buffers,
            one_channel(hbt::PrimitiveChannelId::PiPlusPiPlus),
            hbt::OriginMode::All,
            active_slicing(),
            frame_consumer
        );

    return frame_observable_call_count == 0U &&
           frame_consumer.call_count == 0U &&
           result.summary.pair_counts.channels[0].pair_count == 1U &&
           result.summary.valid_pair_counts.channels[0].pair_count == 0U &&
           result.summary.numerical_rejection_counts
                   .channels[0].pair_count == 1U;
}

/**
 * @brief Verify a valid pair outside active slicing skips frame calculation.
 * @return true when valid origin routing remains and consumer calls stay zero.
 */
bool verify_out_of_slice_pair_skips_frames() {
    hbt::EventBuffers buffers;
    buffers.add(make_particle(
        hbt::SpeciesId::Proton, 0.8, 0.0, 0.938, 2212, 1));
    buffers.add(make_particle(
        hbt::SpeciesId::Proton, 0.8, 0.0, 0.938, 2212, 1));

    reset_frame_state();
    const hbt::PairSubeventProcessingResult result =
        hbt::process_subevent_pairs(
            1U,
            1,
            buffers,
            one_channel(hbt::PrimitiveChannelId::ProtonProton),
            hbt::OriginMode::All,
            active_slicing(),
            frame_consumer
        );

    const hbt::PairOriginRouteCountSummary& routes =
        result.summary.origin_route_counts;
    return frame_observable_call_count == 0U &&
           frame_consumer.call_count == 0U &&
           result.summary.valid_pair_counts.channels[0].pair_count == 1U &&
           routes.routed_P.channels[0].pair_count == 1U &&
           routes.routed_PR.channels[0].pair_count == 1U &&
           routes.routed_PRD.channels[0].pair_count == 1U;
}

/**
 * @brief Verify one in-slice pair is calculated and consumed exactly once.
 * @return true when one call carries all P/PR/PRD routes and the reused slice.
 */
bool verify_in_slice_pair_has_one_consumer_call() {
    hbt::EventBuffers buffers;
    buffers.add(make_particle(
        hbt::SpeciesId::Proton, 0.2, 0.0, 0.938, 2212, 1));
    buffers.add(make_particle(
        hbt::SpeciesId::Proton, 0.2, 0.0, 0.938, 2212, 1));

    reset_frame_state();
    const hbt::PairSubeventProcessingResult result =
        hbt::process_subevent_pairs(
            1U,
            2,
            buffers,
            one_channel(hbt::PrimitiveChannelId::ProtonProton),
            hbt::OriginMode::All,
            active_slicing(),
            frame_consumer
        );

    if (frame_observable_call_count != 1U ||
        frame_consumer.call_count != 1U) {
        return fail("one in-slice pair was not delivered exactly once");
    }
    if (!frame_consumer.routes.primordial ||
        !frame_consumer.routes.primordial_rescattering ||
        !frame_consumer.routes.primordial_rescattering_decay) {
        return fail("one consumer call did not carry all resolved routes");
    }
    if (!frame_consumer.slice_route.has_value() ||
        frame_consumer.slice_route->flat_slice_index != 0U) {
        return fail("consumer did not receive the resolved flat slice");
    }
    return result.summary.valid_pair_counts.channels[0].pair_count == 1U;
}

/**
 * @brief Verify disabled slicing creates no dummy slice destination.
 * @return true when one valid pair is consumed once with a null slice route.
 */
bool verify_disabled_slicing_has_one_consumer_call() {
    hbt::EventBuffers buffers;
    buffers.add(make_particle(
        hbt::SpeciesId::Proton, 0.5, 0.0, 0.938, 2212, 1));
    buffers.add(make_particle(
        hbt::SpeciesId::Proton, -0.1, 0.0, 0.938, 2212, 1));

    reset_frame_state();
    const hbt::PairSubeventProcessingResult result =
        hbt::process_subevent_pairs(
            1U,
            3,
            buffers,
            one_channel(hbt::PrimitiveChannelId::ProtonProton),
            hbt::OriginMode::All,
            disabled_slicing(),
            frame_consumer
        );

    return frame_observable_call_count == 1U &&
           frame_consumer.call_count == 1U &&
           !frame_consumer.slice_route.has_value() &&
           result.summary.valid_pair_counts.channels[0].pair_count == 1U &&
           result.summary.pair_slice_counts.entries.empty();
}

/**
 * @brief Verify every non-finite frame field maps to its exact report reason.
 * @return true when all nine frame reasons reject before consumer delivery.
 */
bool verify_non_finite_frames_are_rejected_before_consumer() {
    const std::array<hbt::PairRejectionReason, 9U> expected_reasons{
        hbt::PairRejectionReason::NonFiniteDeltaTLab,
        hbt::PairRejectionReason::NonFiniteDeltaTLcms,
        hbt::PairRejectionReason::NonFiniteDeltaTPrf,
        hbt::PairRejectionReason::NonFiniteROutLcms,
        hbt::PairRejectionReason::NonFiniteROutPrf,
        hbt::PairRejectionReason::NonFiniteRSide,
        hbt::PairRejectionReason::NonFiniteRLong,
        hbt::PairRejectionReason::NonFiniteRRadialLcms,
        hbt::PairRejectionReason::NonFiniteRRadialPrf
    };
    hbt::EventBuffers buffers;
    buffers.add(make_particle(
        hbt::SpeciesId::Proton, 0.2, 0.0, 0.938, 2212, 1));
    buffers.add(make_particle(
        hbt::SpeciesId::Proton, 0.2, 0.0, 0.938, 2212, 1));

    for (std::size_t index = 0U; index < expected_reasons.size(); ++index) {
        reset_frame_state();
        non_finite_frame_index = index;
        const hbt::PairSubeventProcessingResult result =
            hbt::process_subevent_pairs(
                1U,
                4,
                buffers,
                one_channel(hbt::PrimitiveChannelId::ProtonProton),
                hbt::OriginMode::All,
                active_slicing(),
                frame_consumer
            );

        if (frame_observable_call_count != 1U ||
            frame_consumer.call_count != 0U ||
            result.summary.pair_counts.channels[0].pair_count != 1U ||
            result.summary.valid_pair_counts.channels[0].pair_count != 0U ||
            result.summary.numerical_rejection_counts
                    .channels[0].pair_count != 1U ||
            result.numerical_rejections.size() != 1U ||
            result.numerical_rejections.count(expected_reasons[index]) != 1U ||
            result.numerical_rejections.records()[0].reason !=
                expected_reasons[index]) {
            return fail("frame rejection reason did not cross pair processing");
        }
    }
    return true;
}

}  // namespace

namespace hbt {

/**
 * @internal
 * @brief Link-time test double counting frame-observable calculation requests.
 * @param particle_a First pair particle, unused by this test double.
 * @param particle_b Second pair particle, unused by this test double.
 * @param kinematics Existing pair kinematics, unused by this test double.
 * @return Finite zero observables or one requested non-finite frame field.
 * @endinternal
 */
PairFrameObservables calculate_pair_frame_observables(
    const Particle& particle_a,
    const Particle& particle_b,
    const PairKinematics& kinematics
) noexcept {
    static_cast<void>(particle_a);
    static_cast<void>(particle_b);
    static_cast<void>(kinematics);
    ++frame_observable_call_count;
    PairFrameObservables observables{};
    if (non_finite_frame_index.has_value()) {
        observables.*frame_observable_members[
            non_finite_frame_index.value()] =
                std::numeric_limits<double>::infinity();
    }
    return observables;
}

}  // namespace hbt

/**
 * @brief Run the frame-observable gate integration-test collection.
 * @return EXIT_SUCCESS when every gate scenario passes, otherwise failure.
 */
int main() {
    bool success = true;
    success = verify_kinematic_rejection_skips_frames() && success;
    success = verify_out_of_slice_pair_skips_frames() && success;
    success = verify_in_slice_pair_has_one_consumer_call() && success;
    success = verify_disabled_slicing_has_one_consumer_call() && success;
    success =
        verify_non_finite_frames_are_rejected_before_consumer() && success;
    return success ? EXIT_SUCCESS : EXIT_FAILURE;
}
