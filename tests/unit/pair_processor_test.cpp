/**
 * @file pair_processor_test.cpp
 * @brief Unit tests for pair kinematics policy and exact pair accounting.
 */

#include "hbt/pair/pair_count_accumulator.h"
#include "hbt/pair/pair_processor.h"
#include "hbt/pair/pair_slice_count_accumulator.h"

#include <cmath>
#include <cstddef>
#include <cstdint>
#include <cstdlib>
#include <iostream>
#include <limits>
#include <stdexcept>
#include <utility>
#include <vector>

namespace {

/**
 * @brief Build one synthetic accepted HBT particle for pair-processing tests.
 * @param species Canonical particle species.
 * @param px First transverse momentum component in GeV.
 * @param py Second transverse momentum component in GeV.
 * @param mass Stored invariant mass in GeV.
 * @param pdg Raw signed PDG code retained for diagnostics.
 * @param charge Raw electric charge retained for diagnostics.
 * @param origin Inclusive origin memberships assigned to the particle.
 * @return Complete synthetic Particle.
 */
hbt::Particle make_particle(
    hbt::SpeciesId species,
    double px,
    double py,
    double mass,
    int pdg,
    int charge,
    hbt::OriginFlags origin = {true, true, true}
) {
    return {
        species,
        {0.0, 0.0, 0.0, 0.0},
        {2.0, px, py, 0.0},
        mass,
        origin,
        pdg,
        charge
    };
}

/**
 * @brief Report one failed pair-processor test condition.
 * @param message Description of the violated contract.
 * @return Always false.
 */
bool fail(const char* message) {
    std::cerr << "pair_processor_test: " << message << ".\n";
    return false;
}

/**
 * @brief Test sink used where frame payloads are irrelevant to the assertion.
 */
class DiscardingPairFrameConsumer final : public hbt::PairFrameConsumer {
public:
    /** @brief Consume one frame route without retaining test state. */
    void consume(
        const hbt::PairFrameRouteContext& context,
        const hbt::PairKinematics& kinematics,
        const hbt::PairFrameObservables& observables
    ) override {
        static_cast<void>(context);
        static_cast<void>(kinematics);
        static_cast<void>(observables);
    }
};

/**
 * @brief Process one test subevent through an explicit frame consumer.
 */
hbt::PairSubeventProcessingResult process_subevent_pairs_for_test(
    std::size_t outer_event_number,
    int subevent_id,
    const hbt::EventBuffers& buffers,
    const std::vector<hbt::PrimitiveChannelId>& required_channels,
    hbt::OriginMode origin_mode,
    const hbt::PairSlicingConfig& pair_slicing
) {
    DiscardingPairFrameConsumer frame_consumer;
    return hbt::process_subevent_pairs(
        outer_event_number,
        subevent_id,
        buffers,
        required_channels,
        origin_mode,
        pair_slicing,
        frame_consumer
    );
}

/**
 * @brief Verify one channel count in a three-channel summary.
 * @param summary Ordered three-channel summary to inspect.
 * @param index Expected channel position.
 * @param channel Expected primitive channel identifier.
 * @param count Expected count at that position.
 * @return true when channel and count both match.
 */
bool has_count(
    const hbt::PairCountSummary& summary,
    std::size_t index,
    hbt::PrimitiveChannelId channel,
    std::uint64_t count
) {
    return summary.channels.size() == 3U &&
           summary.channels[index].channel == channel &&
           summary.channels[index].pair_count == count;
}

/**
 * @brief Return a configuration with both kinetic slicing axes disabled.
 * @return Valid no-slicing configuration.
 */
hbt::PairSlicingConfig disabled_slicing() {
    return {{false, {}}, {false, {}}};
}

/**
 * @brief Return one kT slice that contains the valid policy-test pair.
 * @return Valid kT-only slicing configuration.
 */
hbt::PairSlicingConfig policy_slicing() {
    return {{true, {0.0, 0.3}}, {false, {}}};
}

/**
 * @brief Return three kT slices used to test origin/channel routing.
 * @return Valid kT-only slicing configuration with three bins.
 */
hbt::PairSlicingConfig three_slice_slicing() {
    return {{true, {0.0, 0.3, 0.45, 0.6}}, {false, {}}};
}

/**
 * @brief Build buffers producing one kT reject, one mT reject, and one valid.
 * @return Synthetic same-subevent particle buffers.
 */
hbt::EventBuffers make_policy_buffers() {
    hbt::EventBuffers buffers;
    const double maximum = std::numeric_limits<double>::max();

    buffers.add(make_particle(
        hbt::SpeciesId::PiPlus, maximum, 0.0, 0.14, 211, 1));
    buffers.add(make_particle(
        hbt::SpeciesId::PiPlus, maximum, 0.0, 0.14, 211, 1));

    buffers.add(make_particle(
        hbt::SpeciesId::KPlus, 0.5, 0.0, maximum, 321, 1));
    buffers.add(make_particle(
        hbt::SpeciesId::KPlus, 0.5, 0.0, maximum, 321, 1));

    buffers.add(make_particle(
        hbt::SpeciesId::Proton, 0.5, 0.0, 0.938, 2212, 1));
    buffers.add(make_particle(
        hbt::SpeciesId::Proton, -0.1, 0.0, 0.938, 2212, 1));

    return buffers;
}

/**
 * @brief Return the fixed three-channel policy-test order.
 * @return Required primitive channels in deterministic startup order.
 */
std::vector<hbt::PrimitiveChannelId> policy_channels() {
    return {
        hbt::PrimitiveChannelId::PiPlusPiPlus,
        hbt::PrimitiveChannelId::KPlusKPlus,
        hbt::PrimitiveChannelId::ProtonProton
    };
}

/**
 * @brief Build zeroed run-total pair processing state for one origin mode.
 * @param channels Ordered primitive channels.
 * @param mode Configured pair-origin routing mode.
 * @param slicing Validated kinetic slicing configuration.
 * @return Empty run-total summary with matching route structures.
 */
hbt::HBTPairProcessingSummary make_total_summary(
    const std::vector<hbt::PrimitiveChannelId>& channels,
    hbt::OriginMode mode,
    const hbt::PairSlicingConfig& slicing
) {
    return {
        hbt::make_zero_pair_count_summary(channels),
        hbt::make_zero_pair_count_summary(channels),
        hbt::make_zero_pair_count_summary(channels),
        {
            mode,
            hbt::make_zero_pair_count_summary(channels),
            hbt::make_zero_pair_count_summary(channels),
            hbt::make_zero_pair_count_summary(channels)
        },
        hbt::make_zero_pair_slice_count_summary(
            slicing,
            mode,
            channels
        ),
        {},
        {}
    };
}

/**
 * @brief Build zeroed local origin-route counts for one channel list.
 * @param channels Ordered primitive channels.
 * @param mode Configured pair-origin routing mode.
 * @return Empty local route summary.
 */
hbt::PairOriginRouteCountSummary make_zero_route_counts(
    const std::vector<hbt::PrimitiveChannelId>& channels,
    hbt::OriginMode mode
) {
    return {
        mode,
        hbt::make_zero_pair_count_summary(channels),
        hbt::make_zero_pair_count_summary(channels),
        hbt::make_zero_pair_count_summary(channels)
    };
}

/**
 * @brief Verify formed pairs are exactly partitioned into valid and rejected.
 * @return true when counts, reasons, ordinals, and snapshots are exact.
 */
bool verify_local_processing_and_rejection_reporting() {
    const hbt::PairSubeventProcessingResult result =
        process_subevent_pairs_for_test(
            7U,
            19,
            make_policy_buffers(),
            policy_channels(),
            hbt::OriginMode::All,
            policy_slicing()
        );

    const hbt::HBTPairSubeventSummary& summary = result.summary;

    if (summary.outer_event_number != 7U || summary.subevent_id != 19) {
        return fail("subevent identity changed");
    }

    if (!has_count(
            summary.pair_counts,
            0U,
            hbt::PrimitiveChannelId::PiPlusPiPlus,
            1U) ||
        !has_count(
            summary.pair_counts,
            1U,
            hbt::PrimitiveChannelId::KPlusKPlus,
            1U) ||
        !has_count(
            summary.pair_counts,
            2U,
            hbt::PrimitiveChannelId::ProtonProton,
            1U)) {
        return fail("formed-pair counts differ");
    }

    if (!has_count(
            summary.valid_pair_counts,
            0U,
            hbt::PrimitiveChannelId::PiPlusPiPlus,
            0U) ||
        !has_count(
            summary.valid_pair_counts,
            1U,
            hbt::PrimitiveChannelId::KPlusKPlus,
            0U) ||
        !has_count(
            summary.valid_pair_counts,
            2U,
            hbt::PrimitiveChannelId::ProtonProton,
            1U)) {
        return fail("valid-pair counts differ");
    }

    if (!has_count(
            summary.numerical_rejection_counts,
            0U,
            hbt::PrimitiveChannelId::PiPlusPiPlus,
            1U) ||
        !has_count(
            summary.numerical_rejection_counts,
            1U,
            hbt::PrimitiveChannelId::KPlusKPlus,
            1U) ||
        !has_count(
            summary.numerical_rejection_counts,
            2U,
            hbt::PrimitiveChannelId::ProtonProton,
            0U)) {
        return fail("numerical-rejection counts differ");
    }

    const hbt::PairOriginRouteCountSummary& routes =
        summary.origin_route_counts;

    if (routes.origin_mode != hbt::OriginMode::All ||
        !has_count(
            routes.routed_P,
            2U,
            hbt::PrimitiveChannelId::ProtonProton,
            1U) ||
        !has_count(
            routes.routed_PR,
            2U,
            hbt::PrimitiveChannelId::ProtonProton,
            1U) ||
        !has_count(
            routes.routed_PRD,
            2U,
            hbt::PrimitiveChannelId::ProtonProton,
            1U) ||
        routes.routed_P.channels[0].pair_count != 0U ||
        routes.routed_P.channels[1].pair_count != 0U ||
        routes.routed_PR.channels[0].pair_count != 0U ||
        routes.routed_PR.channels[1].pair_count != 0U ||
        routes.routed_PRD.channels[0].pair_count != 0U ||
        routes.routed_PRD.channels[1].pair_count != 0U) {
        return fail("local origin-route counts differ");
    }

    const hbt::PairSliceCountSummary& slices =
        summary.pair_slice_counts;
    if (slices.entries.size() != 1U ||
        slices.entries[0].kt_slice_index != 0U ||
        slices.entries[0].mt_slice_index.has_value() ||
        !has_count(
            slices.entries[0].origin_counts.routed_P,
            2U,
            hbt::PrimitiveChannelId::ProtonProton,
            1U) ||
        !has_count(
            slices.entries[0].origin_counts.routed_PR,
            2U,
            hbt::PrimitiveChannelId::ProtonProton,
            1U) ||
        !has_count(
            slices.entries[0].origin_counts.routed_PRD,
            2U,
            hbt::PrimitiveChannelId::ProtonProton,
            1U)) {
        return fail("rejected pairs entered slice accounting");
    }

    const hbt::RejectedPairReport& report = result.numerical_rejections;

    if (report.size() != 2U ||
        report.count(hbt::PairRejectionReason::NonFiniteKt) != 1U ||
        report.count(hbt::PairRejectionReason::NonFiniteMt) != 1U) {
        return fail("rejection-reason counts differ");
    }

    const hbt::RejectedPairRecord& kt_record = report.records()[0];
    const hbt::RejectedPairRecord& mt_record = report.records()[1];

    if (kt_record.outer_event_number != 7U ||
        kt_record.subevent_id != 19 ||
        kt_record.channel != hbt::PrimitiveChannelId::PiPlusPiPlus ||
        kt_record.pair_ordinal_in_channel != 1U ||
        kt_record.reason != hbt::PairRejectionReason::NonFiniteKt ||
        std::isfinite(kt_record.kinematics.kt_gev) ||
        kt_record.particle_a.raw_pdg != 211 ||
        kt_record.particle_b.raw_charge != 1) {
        return fail("non-finite-kT rejection record differs");
    }

    if (mt_record.outer_event_number != 7U ||
        mt_record.subevent_id != 19 ||
        mt_record.channel != hbt::PrimitiveChannelId::KPlusKPlus ||
        mt_record.pair_ordinal_in_channel != 1U ||
        mt_record.reason != hbt::PairRejectionReason::NonFiniteMt ||
        !std::isfinite(mt_record.kinematics.kt_gev) ||
        std::isfinite(mt_record.kinematics.mt_gev) ||
        mt_record.particle_a.invariant_mass_gev !=
            std::numeric_limits<double>::max()) {
        return fail("non-finite-mT rejection record differs");
    }

    return true;
}

/**
 * @brief Verify run totals retain the exact accounting partition and records.
 * @return true when two subevents accumulate without losing any rejection.
 */
bool verify_run_total_accumulation() {
    const std::vector<hbt::PrimitiveChannelId> channels = policy_channels();
    hbt::HBTPairProcessingSummary total =
        make_total_summary(
            channels,
            hbt::OriginMode::All,
            policy_slicing()
        );

    hbt::accumulate_pair_processing_result(
        total,
        process_subevent_pairs_for_test(
            1U,
            3,
            make_policy_buffers(),
            channels,
            hbt::OriginMode::All,
            policy_slicing()
        )
    );
    hbt::accumulate_pair_processing_result(
        total,
        process_subevent_pairs_for_test(
            1U,
            4,
            make_policy_buffers(),
            channels,
            hbt::OriginMode::All,
            policy_slicing()
        )
    );

    return
        has_count(
            total.total_pair_counts,
            0U,
            hbt::PrimitiveChannelId::PiPlusPiPlus,
            2U) &&
        has_count(
            total.total_valid_pair_counts,
            2U,
            hbt::PrimitiveChannelId::ProtonProton,
            2U) &&
        has_count(
            total.total_numerical_rejection_counts,
            0U,
            hbt::PrimitiveChannelId::PiPlusPiPlus,
            2U) &&
        has_count(
            total.total_numerical_rejection_counts,
            1U,
            hbt::PrimitiveChannelId::KPlusKPlus,
            2U) &&
        has_count(
            total.total_origin_route_counts.routed_P,
            2U,
            hbt::PrimitiveChannelId::ProtonProton,
            2U) &&
        has_count(
            total.total_origin_route_counts.routed_PR,
            2U,
            hbt::PrimitiveChannelId::ProtonProton,
            2U) &&
        has_count(
            total.total_origin_route_counts.routed_PRD,
            2U,
            hbt::PrimitiveChannelId::ProtonProton,
            2U) &&
        total.total_pair_slice_counts.entries.size() == 1U &&
        has_count(
            total.total_pair_slice_counts.entries[0]
                .origin_counts.routed_P,
            2U,
            hbt::PrimitiveChannelId::ProtonProton,
            2U) &&
        has_count(
            total.total_pair_slice_counts.entries[0]
                .origin_counts.routed_PR,
            2U,
            hbt::PrimitiveChannelId::ProtonProton,
            2U) &&
        has_count(
            total.total_pair_slice_counts.entries[0]
                .origin_counts.routed_PRD,
            2U,
            hbt::PrimitiveChannelId::ProtonProton,
            2U) &&
        total.numerical_rejections.size() == 4U &&
        total.numerical_rejections.count(
            hbt::PairRejectionReason::NonFiniteKt) == 2U &&
        total.numerical_rejections.count(
            hbt::PairRejectionReason::NonFiniteMt) == 2U &&
        total.subevents.size() == 2U;
}

/**
 * @brief Verify an inconsistent local partition fails before mutating totals.
 * @return true when std::logic_error is raised and totals remain untouched.
 */
bool verify_inconsistent_partition_is_structural_failure() {
    const std::vector<hbt::PrimitiveChannelId> channels{
        hbt::PrimitiveChannelId::PiPlusPiPlus
    };
    hbt::HBTPairProcessingSummary total =
        make_total_summary(
            channels,
            hbt::OriginMode::All,
            disabled_slicing()
        );

    hbt::PairSubeventProcessingResult invalid{
        {
            1U,
            0,
            {{{hbt::PrimitiveChannelId::PiPlusPiPlus, 1U}}},
            {{{hbt::PrimitiveChannelId::PiPlusPiPlus, 1U}}},
            {{{hbt::PrimitiveChannelId::PiPlusPiPlus, 1U}}},
            {
                hbt::OriginMode::All,
                {{{hbt::PrimitiveChannelId::PiPlusPiPlus, 1U}}},
                {{{hbt::PrimitiveChannelId::PiPlusPiPlus, 1U}}},
                {{{hbt::PrimitiveChannelId::PiPlusPiPlus, 1U}}}
            },
            hbt::make_zero_pair_slice_count_summary(
                disabled_slicing(),
                hbt::OriginMode::All,
                channels
            )
        },
        {}
    };

    try {
        hbt::accumulate_pair_processing_result(total, std::move(invalid));
    } catch (const std::logic_error&) {
        return total.total_pair_counts.channels[0].pair_count == 0U &&
               total.total_valid_pair_counts.channels[0].pair_count == 0U &&
               total.total_numerical_rejection_counts.channels[0].pair_count ==
                   0U &&
               total.numerical_rejections.empty() &&
               total.subevents.empty();
    } catch (...) {
        return fail("inconsistent partition changed exception type");
    }

    return fail("inconsistent partition was not rejected");
}

/**
 * @brief Verify duplicate rejected-pair identities are structural failures.
 * @return true when duplicate channel/ordinal records are rejected pre-commit.
 */
bool verify_duplicate_rejection_identity_is_structural_failure() {
    const std::vector<hbt::PrimitiveChannelId> channels{
        hbt::PrimitiveChannelId::PiPlusPiPlus
    };
    hbt::HBTPairProcessingSummary total =
        make_total_summary(
            channels,
            hbt::OriginMode::All,
            disabled_slicing()
        );

    const hbt::RejectedPairParticleSnapshot particle{
        hbt::SpeciesId::PiPlus,
        {2.0, 1.0, 0.0, 0.0},
        0.139,
        211,
        1
    };
    const hbt::RejectedPairRecord duplicate{
        1U,
        0,
        hbt::PrimitiveChannelId::PiPlusPiPlus,
        1U,
        particle,
        particle,
        {
            {std::numeric_limits<double>::infinity(), 0.0, 0.0, 0.0},
            0.0,
            0.0,
            std::numeric_limits<double>::infinity(),
            1.0
        },
        hbt::PairRejectionReason::NonFiniteKt
    };

    hbt::RejectedPairReport report;
    report.add(duplicate);
    report.add(duplicate);

    hbt::PairSubeventProcessingResult invalid{
        {
            1U,
            0,
            {{{hbt::PrimitiveChannelId::PiPlusPiPlus, 2U}}},
            {{{hbt::PrimitiveChannelId::PiPlusPiPlus, 0U}}},
            {{{hbt::PrimitiveChannelId::PiPlusPiPlus, 2U}}},
            make_zero_route_counts(channels, hbt::OriginMode::All),
            hbt::make_zero_pair_slice_count_summary(
                disabled_slicing(),
                hbt::OriginMode::All,
                channels
            )
        },
        std::move(report)
    };

    try {
        hbt::accumulate_pair_processing_result(total, std::move(invalid));
    } catch (const std::logic_error&) {
        return total.total_pair_counts.channels[0].pair_count == 0U &&
               total.numerical_rejections.empty() &&
               total.subevents.empty();
    } catch (...) {
        return fail("duplicate rejection identity changed exception type");
    }

    return fail("duplicate rejection identity was not rejected");
}


/**
 * @brief Verify every individual origin mode routes only its requested slice.
 * @return true when one valid primordial pair follows each exact mode contract.
 */
bool verify_individual_origin_mode_routing() {
    const std::vector<hbt::PrimitiveChannelId> channels{
        hbt::PrimitiveChannelId::ProtonProton
    };
    hbt::EventBuffers buffers;
    buffers.add(make_particle(
        hbt::SpeciesId::Proton, 0.5, 0.0, 0.938, 2212, 1));
    buffers.add(make_particle(
        hbt::SpeciesId::Proton, -0.1, 0.0, 0.938, 2212, 1));

    const struct {
        hbt::OriginMode mode;
        std::uint64_t p;
        std::uint64_t pr;
        std::uint64_t prd;
    } cases[]{
        {hbt::OriginMode::Primordial, 1U, 0U, 0U},
        {hbt::OriginMode::PrimordialRescattering, 0U, 1U, 0U},
        {hbt::OriginMode::PrimordialRescatteringDecay, 0U, 0U, 1U},
        {hbt::OriginMode::All, 1U, 1U, 1U}
    };

    for (const auto& test_case : cases) {
        const hbt::PairSubeventProcessingResult result =
            process_subevent_pairs_for_test(
                1U,
                0,
                buffers,
                channels,
                test_case.mode,
                policy_slicing()
            );
        const hbt::PairOriginRouteCountSummary& routes =
            result.summary.origin_route_counts;

        const hbt::PairSliceCountSummary& slices =
            result.summary.pair_slice_counts;

        if (routes.origin_mode != test_case.mode ||
            routes.routed_P.channels[0].pair_count != test_case.p ||
            routes.routed_PR.channels[0].pair_count != test_case.pr ||
            routes.routed_PRD.channels[0].pair_count != test_case.prd ||
            slices.entries.size() != 1U ||
            slices.entries[0].origin_counts.routed_P
                    .channels[0].pair_count != test_case.p ||
            slices.entries[0].origin_counts.routed_PR
                    .channels[0].pair_count != test_case.pr ||
            slices.entries[0].origin_counts.routed_PRD
                    .channels[0].pair_count != test_case.prd) {
            return fail("individual origin/slice routing differs");
        }
    }

    return true;
}

/**
 * @brief Verify one slice lookup is reused across compatible origin routes.
 * @return true when slice/origin/channel counts match the routed pairs exactly.
 */
bool verify_slice_origin_channel_routing() {
    const std::vector<hbt::PrimitiveChannelId> channels = policy_channels();
    const hbt::PairSlicingConfig slicing = three_slice_slicing();
    const hbt::OriginFlags primordial{true, true, true};
    const hbt::OriginFlags rescattering{false, true, true};
    const hbt::OriginFlags decay{false, false, true};

    hbt::EventBuffers buffers;
    buffers.add(make_particle(
        hbt::SpeciesId::PiPlus,
        0.2,
        0.0,
        0.14,
        211,
        1,
        primordial));
    buffers.add(make_particle(
        hbt::SpeciesId::PiPlus,
        0.2,
        0.0,
        0.14,
        211,
        1,
        primordial));
    buffers.add(make_particle(
        hbt::SpeciesId::KPlus,
        0.4,
        0.0,
        0.494,
        321,
        1,
        primordial));
    buffers.add(make_particle(
        hbt::SpeciesId::KPlus,
        0.4,
        0.0,
        0.494,
        321,
        1,
        rescattering));
    buffers.add(make_particle(
        hbt::SpeciesId::Proton,
        0.5,
        0.0,
        0.938,
        2212,
        1,
        decay));
    buffers.add(make_particle(
        hbt::SpeciesId::Proton,
        0.5,
        0.0,
        0.938,
        2212,
        1,
        decay));

    const hbt::PairSubeventProcessingResult result =
        process_subevent_pairs_for_test(
            2U,
            5,
            buffers,
            channels,
            hbt::OriginMode::All,
            slicing
        );

    const hbt::PairOriginRouteCountSummary& routes =
        result.summary.origin_route_counts;
    const hbt::PairSliceCountSummary& slices =
        result.summary.pair_slice_counts;

    if (!has_count(
            routes.routed_P,
            0U,
            hbt::PrimitiveChannelId::PiPlusPiPlus,
            1U) ||
        !has_count(
            routes.routed_PR,
            1U,
            hbt::PrimitiveChannelId::KPlusKPlus,
            1U) ||
        !has_count(
            routes.routed_PRD,
            2U,
            hbt::PrimitiveChannelId::ProtonProton,
            1U) ||
        slices.entries.size() != 3U) {
        return fail("origin routes before slicing differ");
    }

    const hbt::PairSliceCountEntry& p_slice = slices.entries[0];
    const hbt::PairSliceCountEntry& pr_slice = slices.entries[1];
    const hbt::PairSliceCountEntry& prd_slice = slices.entries[2];

    if (p_slice.kt_slice_index != 0U ||
        !has_count(
            p_slice.origin_counts.routed_P,
            0U,
            hbt::PrimitiveChannelId::PiPlusPiPlus,
            1U) ||
        !has_count(
            p_slice.origin_counts.routed_PR,
            0U,
            hbt::PrimitiveChannelId::PiPlusPiPlus,
            1U) ||
        !has_count(
            p_slice.origin_counts.routed_PRD,
            0U,
            hbt::PrimitiveChannelId::PiPlusPiPlus,
            1U)) {
        return fail("P-P pair did not reuse one slice for P/PR/PRD");
    }

    if (pr_slice.kt_slice_index != 1U ||
        !has_count(
            pr_slice.origin_counts.routed_P,
            1U,
            hbt::PrimitiveChannelId::KPlusKPlus,
            0U) ||
        !has_count(
            pr_slice.origin_counts.routed_PR,
            1U,
            hbt::PrimitiveChannelId::KPlusKPlus,
            1U) ||
        !has_count(
            pr_slice.origin_counts.routed_PRD,
            1U,
            hbt::PrimitiveChannelId::KPlusKPlus,
            1U)) {
        return fail("P-R pair slice routes differ");
    }

    if (prd_slice.kt_slice_index != 2U ||
        !has_count(
            prd_slice.origin_counts.routed_P,
            2U,
            hbt::PrimitiveChannelId::ProtonProton,
            0U) ||
        !has_count(
            prd_slice.origin_counts.routed_PR,
            2U,
            hbt::PrimitiveChannelId::ProtonProton,
            0U) ||
        !has_count(
            prd_slice.origin_counts.routed_PRD,
            2U,
            hbt::PrimitiveChannelId::ProtonProton,
            1U)) {
        return fail("D-D pair slice routes differ");
    }

    hbt::HBTPairProcessingSummary total =
        make_total_summary(channels, hbt::OriginMode::All, slicing);
    hbt::accumulate_pair_processing_result(total, result);
    hbt::accumulate_pair_processing_result(
        total,
        process_subevent_pairs_for_test(
            2U,
            6,
            buffers,
            channels,
            hbt::OriginMode::All,
            slicing
        )
    );

    return total.subevents.size() == 2U &&
           has_count(
               total.total_pair_slice_counts.entries[0]
                   .origin_counts.routed_P,
               0U,
               hbt::PrimitiveChannelId::PiPlusPiPlus,
               2U) &&
           has_count(
               total.total_pair_slice_counts.entries[1]
                   .origin_counts.routed_PR,
               1U,
               hbt::PrimitiveChannelId::KPlusKPlus,
               2U) &&
           has_count(
               total.total_pair_slice_counts.entries[2]
                   .origin_counts.routed_PRD,
               2U,
               hbt::PrimitiveChannelId::ProtonProton,
               2U);
}

/**
 * @brief Verify out-of-range pairs remain valid/routed but enter no slice.
 * @return true when kinetic coverage does not change pair validity.
 */
bool verify_out_of_slice_pair_remains_valid() {
    const std::vector<hbt::PrimitiveChannelId> channels{
        hbt::PrimitiveChannelId::ProtonProton
    };
    const hbt::OriginFlags decay{false, false, true};
    hbt::EventBuffers buffers;
    buffers.add(make_particle(
        hbt::SpeciesId::Proton, 0.8, 0.0, 0.938, 2212, 1, decay));
    buffers.add(make_particle(
        hbt::SpeciesId::Proton, 0.8, 0.0, 0.938, 2212, 1, decay));

    const hbt::PairSubeventProcessingResult result =
        process_subevent_pairs_for_test(
            3U,
            1,
            buffers,
            channels,
            hbt::OriginMode::All,
            three_slice_slicing()
        );

    if (result.summary.valid_pair_counts.channels[0].pair_count != 1U ||
        result.summary.origin_route_counts.routed_PRD
                .channels[0].pair_count != 1U) {
        return fail("out-of-slice pair stopped being valid/routed");
    }

    for (const hbt::PairSliceCountEntry& entry :
         result.summary.pair_slice_counts.entries) {
        if (entry.origin_counts.routed_P.channels[0].pair_count != 0U ||
            entry.origin_counts.routed_PR.channels[0].pair_count != 0U ||
            entry.origin_counts.routed_PRD.channels[0].pair_count != 0U) {
            return fail("out-of-slice pair incremented slice counts");
        }
    }

    return true;
}

/**
 * @brief Verify disabled slicing creates no production dummy slice.
 * @return true when valid origin-routed pairs leave slice entries empty.
 */
bool verify_disabled_slicing_has_no_dummy_entries() {
    const std::vector<hbt::PrimitiveChannelId> channels{
        hbt::PrimitiveChannelId::ProtonProton
    };
    hbt::EventBuffers buffers;
    buffers.add(make_particle(
        hbt::SpeciesId::Proton, 0.5, 0.0, 0.938, 2212, 1));
    buffers.add(make_particle(
        hbt::SpeciesId::Proton, -0.1, 0.0, 0.938, 2212, 1));

    const hbt::PairSubeventProcessingResult result =
        process_subevent_pairs_for_test(
            4U,
            2,
            buffers,
            channels,
            hbt::OriginMode::All,
            disabled_slicing()
        );

    return result.summary.valid_pair_counts.channels[0].pair_count == 1U &&
           result.summary.origin_route_counts.routed_PRD
                   .channels[0].pair_count == 1U &&
           result.summary.pair_slice_counts.entries.empty();
}

/**
 * @brief Verify slice counts cannot exceed their pre-slicing origin routes.
 * @return true when a double-counted local slice fails before total mutation.
 */
bool verify_slice_overcount_is_structural_failure() {
    const std::vector<hbt::PrimitiveChannelId> channels{
        hbt::PrimitiveChannelId::ProtonProton
    };
    const hbt::PairSlicingConfig slicing = policy_slicing();
    hbt::EventBuffers buffers;
    buffers.add(make_particle(
        hbt::SpeciesId::Proton, 0.5, 0.0, 0.938, 2212, 1));
    buffers.add(make_particle(
        hbt::SpeciesId::Proton, -0.1, 0.0, 0.938, 2212, 1));

    hbt::PairSubeventProcessingResult local =
        process_subevent_pairs_for_test(
            1U,
            0,
            buffers,
            channels,
            hbt::OriginMode::All,
            slicing
        );
    local.summary.pair_slice_counts.entries[0]
        .origin_counts.routed_PRD.channels[0].pair_count = 2U;

    hbt::HBTPairProcessingSummary total =
        make_total_summary(channels, hbt::OriginMode::All, slicing);

    try {
        hbt::accumulate_pair_processing_result(total, std::move(local));
    } catch (const std::logic_error&) {
        return total.total_pair_counts.channels[0].pair_count == 0U &&
               total.total_pair_slice_counts.entries[0]
                       .origin_counts.routed_PRD
                       .channels[0].pair_count == 0U &&
               total.subevents.empty();
    } catch (...) {
        return fail("slice overcount changed exception type");
    }

    return fail("slice overcount was not rejected");
}

/**
 * @brief Verify event-local pair summaries reduce in caller-provided order.
 * @return true when totals and subevent identity are preserved exactly.
 */
bool verify_event_summary_reduction() {
    const std::vector<hbt::PrimitiveChannelId> channels{
        hbt::PrimitiveChannelId::ProtonProton
    };
    hbt::EventBuffers first_buffers;
    first_buffers.add(make_particle(
        hbt::SpeciesId::Proton, 0.5, 0.0, 0.938, 2212, 1));
    first_buffers.add(make_particle(
        hbt::SpeciesId::Proton, -0.1, 0.0, 0.938, 2212, 1));
    hbt::EventBuffers second_buffers = first_buffers;

    hbt::HBTPairProcessingSummary first = make_total_summary(
        channels,
        hbt::OriginMode::All,
        disabled_slicing()
    );
    hbt::HBTPairProcessingSummary second = make_total_summary(
        channels,
        hbt::OriginMode::All,
        disabled_slicing()
    );
    hbt::accumulate_pair_processing_result(
        first,
        process_subevent_pairs_for_test(
            1U,
            0,
            first_buffers,
            channels,
            hbt::OriginMode::All,
            disabled_slicing()
        )
    );
    hbt::accumulate_pair_processing_result(
        second,
        process_subevent_pairs_for_test(
            2U,
            0,
            second_buffers,
            channels,
            hbt::OriginMode::All,
            disabled_slicing()
        )
    );

    hbt::HBTPairProcessingSummary total = make_total_summary(
        channels,
        hbt::OriginMode::All,
        disabled_slicing()
    );
    hbt::accumulate_pair_processing_summary(total, std::move(first));
    hbt::accumulate_pair_processing_summary(total, std::move(second));

    return total.total_pair_counts.channels.size() == 1U &&
           total.total_pair_counts.channels[0].pair_count == 2U &&
           total.subevents.size() == 2U &&
           total.subevents[0].outer_event_number == 1U &&
           total.subevents[1].outer_event_number == 2U;
}

/**
 * @brief Verify a valid pair with no requested route is structural failure.
 * @return true when incompatible particle origins do not disappear silently.
 */
bool verify_missing_requested_route_is_structural_failure() {
    const std::vector<hbt::PrimitiveChannelId> channels{
        hbt::PrimitiveChannelId::ProtonProton
    };
    hbt::EventBuffers buffers;
    hbt::Particle decay = make_particle(
        hbt::SpeciesId::Proton, 0.5, 0.0, 0.938, 2212, 1);
    decay.origin = {false, false, true};
    buffers.add(decay);
    buffers.add(decay);

    try {
        static_cast<void>(process_subevent_pairs_for_test(
            1U,
            0,
            buffers,
            channels,
            hbt::OriginMode::Primordial,
            disabled_slicing()
        ));
    } catch (const std::logic_error&) {
        return true;
    } catch (...) {
        return fail("missing-route failure changed exception type");
    }

    return fail("valid pair without requested route was hidden");
}

/**
 * @brief Verify local/run-total OriginMode mismatch is structural failure.
 * @return true when totals remain unchanged after the mismatch is rejected.
 */
bool verify_origin_mode_mismatch_is_structural_failure() {
    const std::vector<hbt::PrimitiveChannelId> channels{
        hbt::PrimitiveChannelId::ProtonProton
    };
    hbt::EventBuffers buffers;
    buffers.add(make_particle(
        hbt::SpeciesId::Proton, 0.5, 0.0, 0.938, 2212, 1));
    buffers.add(make_particle(
        hbt::SpeciesId::Proton, -0.1, 0.0, 0.938, 2212, 1));

    hbt::HBTPairProcessingSummary total =
        make_total_summary(
            channels,
            hbt::OriginMode::All,
            disabled_slicing()
        );
    hbt::PairSubeventProcessingResult local =
        process_subevent_pairs_for_test(
            1U,
            0,
            buffers,
            channels,
            hbt::OriginMode::Primordial,
            disabled_slicing()
        );

    try {
        hbt::accumulate_pair_processing_result(total, std::move(local));
    } catch (const std::invalid_argument&) {
        return total.total_pair_counts.channels[0].pair_count == 0U &&
               total.total_origin_route_counts.routed_P
                       .channels[0].pair_count == 0U &&
               total.subevents.empty();
    } catch (...) {
        return fail("origin-mode mismatch changed exception type");
    }

    return fail("origin-mode mismatch was not rejected");
}

}  // namespace

/**
 * @brief Run the complete pair-processor unit-test collection.
 * @return EXIT_SUCCESS when every test passes, otherwise EXIT_FAILURE.
 */
int main() {
    bool success = true;
    success = verify_local_processing_and_rejection_reporting() && success;
    success = verify_run_total_accumulation() && success;
    success = verify_inconsistent_partition_is_structural_failure() && success;
    success =
        verify_duplicate_rejection_identity_is_structural_failure() && success;
    success = verify_individual_origin_mode_routing() && success;
    success = verify_slice_origin_channel_routing() && success;
    success = verify_out_of_slice_pair_remains_valid() && success;
    success = verify_disabled_slicing_has_no_dummy_entries() && success;
    success = verify_slice_overcount_is_structural_failure() && success;
    success = verify_event_summary_reduction() && success;
    success = verify_missing_requested_route_is_structural_failure() && success;
    success = verify_origin_mode_mismatch_is_structural_failure() && success;
    return success ? EXIT_SUCCESS : EXIT_FAILURE;
}
