/**
 * @file pair_processor.cpp
 * @brief Pair-level kinematics, validation, origin routing, and accounting.
 */

#include "hbt/pair/pair_processor.h"

#include "hbt/pair/pair_count_accumulator.h"
#include "hbt/pair/pair_iterator.h"
#include "hbt/pair/pair_frame_observables.h"
#include "hbt/pair/pair_frame_observables_validation.h"
#include "hbt/pair/pair_kinematics.h"
#include "hbt/pair/pair_kinematics_validation.h"
#include "hbt/pair/pair_origin_routing.h"
#include "hbt/pair/pair_slice_count_accumulator.h"
#include "hbt/pair/pair_slice_routing.h"

#include <cstddef>
#include <cstdint>
#include <limits>
#include <optional>
#include <set>
#include <stdexcept>
#include <string>
#include <utility>

namespace hbt {
namespace {

/**
 * @brief Return the mutable count entry for one channel.
 * @param summary Ordered channel-count summary to search.
 * @param channel Channel whose entry is required.
 * @return Mutable matching entry.
 * @throws std::logic_error If the channel is absent from the summary.
 */
PairChannelCount& count_entry_for_channel(
    PairCountSummary& summary,
    PrimitiveChannelId channel
) {
    for (PairChannelCount& entry : summary.channels) {
        if (entry.channel == channel) {
            return entry;
        }
    }

    throw std::logic_error(
        "pair processing: iterator returned an unexpected channel"
    );
}

/**
 * @brief Return the immutable count entry for one channel.
 * @param summary Ordered channel-count summary to search.
 * @param channel Channel whose entry is required.
 * @return Immutable matching entry.
 * @throws std::logic_error If the channel is absent from the summary.
 */
const PairChannelCount& count_entry_for_channel(
    const PairCountSummary& summary,
    PrimitiveChannelId channel
) {
    for (const PairChannelCount& entry : summary.channels) {
        if (entry.channel == channel) {
            return entry;
        }
    }

    throw std::logic_error(
        "pair processing: summary does not contain rejection channel"
    );
}

/**
 * @brief Increment one uint64 pair counter with explicit overflow checking.
 * @param count Counter to increment.
 * @throws std::overflow_error If the increment would overflow.
 */
void increment_pair_count(std::uint64_t& count) {
    if (count == std::numeric_limits<std::uint64_t>::max()) {
        throw std::overflow_error("pair processing: local count overflow");
    }

    ++count;
}

/**
 * @brief Verify formed = valid + numerical rejected for every channel.
 * @param formed All formed physical pair counts.
 * @param valid Numerically valid pair counts.
 * @param rejected Numerical pair-rejection counts.
 * @throws std::logic_error If any channel structure or count identity differs.
 */
void require_pair_count_partition(
    const PairCountSummary& formed,
    const PairCountSummary& valid,
    const PairCountSummary& rejected
) {
    if (formed.channels.size() != valid.channels.size() ||
        formed.channels.size() != rejected.channels.size()) {
        throw std::logic_error(
            "pair processing: formed/valid/rejected size mismatch"
        );
    }

    for (std::size_t index = 0U; index < formed.channels.size(); ++index) {
        const PairChannelCount& formed_entry = formed.channels[index];
        const PairChannelCount& valid_entry = valid.channels[index];
        const PairChannelCount& rejected_entry = rejected.channels[index];

        if (formed_entry.channel != valid_entry.channel ||
            formed_entry.channel != rejected_entry.channel) {
            throw std::logic_error(
                "pair processing: formed/valid/rejected channel mismatch"
            );
        }

        if (valid_entry.pair_count > formed_entry.pair_count ||
            rejected_entry.pair_count >
                formed_entry.pair_count - valid_entry.pair_count ||
            valid_entry.pair_count + rejected_entry.pair_count !=
                formed_entry.pair_count) {
            throw std::logic_error(
                "pair processing: formed != valid + numerical rejected"
            );
        }
    }
}

/**
 * @brief Verify rejection records exactly reproduce one rejected-count summary.
 * @param rejected Numerical rejection counts by primitive channel.
 * @param report Complete rejection records to compare.
 * @throws std::logic_error If record channels do not reproduce the counts.
 */
void require_rejection_report_matches_counts(
    const PairCountSummary& rejected,
    const RejectedPairReport& report
) {
    PairCountSummary report_counts;
    report_counts.channels.reserve(rejected.channels.size());

    for (const PairChannelCount& entry : rejected.channels) {
        report_counts.channels.push_back({entry.channel, 0U});
    }

    for (const RejectedPairRecord& record : report.records()) {
        PairChannelCount& entry =
            count_entry_for_channel(report_counts, record.channel);
        increment_pair_count(entry.pair_count);
    }

    if (report_counts.channels.size() != rejected.channels.size()) {
        throw std::logic_error(
            "pair processing: rejection report channel-size mismatch"
        );
    }

    for (std::size_t index = 0U; index < rejected.channels.size(); ++index) {
        if (report_counts.channels[index].channel !=
                rejected.channels[index].channel ||
            report_counts.channels[index].pair_count !=
                rejected.channels[index].pair_count) {
            throw std::logic_error(
                "pair processing: rejection report/count mismatch"
            );
        }
    }
}

/**
 * @brief Verify local rejection identities belong to the local subevent.
 * @param summary Local subevent summary owning the records.
 * @param report Local rejection records to inspect.
 * @throws std::logic_error If any record carries different event identity.
 */
void require_local_rejection_identity(
    const HBTPairSubeventSummary& summary,
    const RejectedPairReport& report
) {
    std::set<std::pair<PrimitiveChannelId, std::uint64_t>> identities;

    for (const RejectedPairRecord& record : report.records()) {
        if (record.outer_event_number != summary.outer_event_number ||
            record.subevent_id != summary.subevent_id) {
            throw std::logic_error(
                "pair processing: rejection record subevent mismatch"
            );
        }

        const PairChannelCount& formed = count_entry_for_channel(
            summary.pair_counts,
            record.channel
        );

        if (record.pair_ordinal_in_channel == 0U ||
            record.pair_ordinal_in_channel > formed.pair_count) {
            throw std::logic_error(
                "pair processing: rejection pair ordinal out of range"
            );
        }

        const auto identity = std::make_pair(
            record.channel,
            record.pair_ordinal_in_channel
        );

        if (!identities.insert(identity).second) {
            throw std::logic_error(
                "pair processing: duplicate rejected-pair identity"
            );
        }
    }
}

/**
 * @brief Verify one route summary has the same ordered channel structure.
 * @param valid Numerically valid pair counts.
 * @param routed Routed pair counts to compare.
 * @param label Human-readable route label for diagnostics.
 * @throws std::logic_error If size or channel order differs.
 */
void require_route_structure(
    const PairCountSummary& valid,
    const PairCountSummary& routed,
    const char* label
) {
    if (valid.channels.size() != routed.channels.size()) {
        throw std::logic_error(
            std::string("pair processing: ") + label +
            " route size mismatch"
        );
    }

    for (std::size_t index = 0U; index < valid.channels.size(); ++index) {
        if (valid.channels[index].channel != routed.channels[index].channel) {
            throw std::logic_error(
                std::string("pair processing: ") + label +
                " route channel mismatch"
            );
        }
    }
}

/**
 * @brief Verify requested-origin routing counts against valid pairs.
 * @param valid Numerically valid pair counts.
 * @param routes Routed pair counts and the mode that produced them.
 * @throws std::logic_error If route counts violate the configured contract.
 * @throws std::invalid_argument If the stored OriginMode value is invalid.
 *
 * Individual origin modes route every valid pair to exactly their one requested
 * selection. OriginMode::All keeps nested overlapping routes and therefore
 * requires routed_P <= routed_PR <= routed_PRD = valid for each channel.
 */
void require_origin_route_contract(
    const PairCountSummary& valid,
    const PairOriginRouteCountSummary& routes
) {
    require_route_structure(valid, routes.routed_P, "P");
    require_route_structure(valid, routes.routed_PR, "PR");
    require_route_structure(valid, routes.routed_PRD, "PRD");

    for (std::size_t index = 0U; index < valid.channels.size(); ++index) {
        const std::uint64_t valid_count = valid.channels[index].pair_count;
        const std::uint64_t routed_P =
            routes.routed_P.channels[index].pair_count;
        const std::uint64_t routed_PR =
            routes.routed_PR.channels[index].pair_count;
        const std::uint64_t routed_PRD =
            routes.routed_PRD.channels[index].pair_count;

        switch (routes.origin_mode) {
        case OriginMode::Primordial:
            if (routed_P != valid_count ||
                routed_PR != 0U ||
                routed_PRD != 0U) {
                throw std::logic_error(
                    "pair processing: invalid Primordial route counts"
                );
            }
            break;
        case OriginMode::PrimordialRescattering:
            if (routed_P != 0U ||
                routed_PR != valid_count ||
                routed_PRD != 0U) {
                throw std::logic_error(
                    "pair processing: invalid PR route counts"
                );
            }
            break;
        case OriginMode::PrimordialRescatteringDecay:
            if (routed_P != 0U ||
                routed_PR != 0U ||
                routed_PRD != valid_count) {
                throw std::logic_error(
                    "pair processing: invalid PRD route counts"
                );
            }
            break;
        case OriginMode::All:
            if (routed_P > routed_PR ||
                routed_PR > routed_PRD ||
                routed_PRD != valid_count) {
                throw std::logic_error(
                    "pair processing: invalid All origin route counts"
                );
            }
            break;
        default:
            throw std::invalid_argument(
                "pair processing: invalid origin mode in route summary"
            );
        }
    }
}

/**
 * @brief Verify slice counts are a routed subset with exact structure.
 * @param routes Run/local origin-route counts available before slicing.
 * @param slices Pair counts grouped by kinetic slice and origin route.
 * @throws std::invalid_argument If stored routing policy is invalid.
 * @throws std::logic_error If slice layout or counts violate routing policy.
 */
void require_pair_slice_count_contract(
    const PairOriginRouteCountSummary& routes,
    const PairSliceCountSummary& slices
) {
    if (routes.origin_mode != slices.origin_mode) {
        throw std::invalid_argument(
            "pair processing: origin/slice mode mismatch"
        );
    }

    std::vector<PrimitiveChannelId> channels;
    channels.reserve(routes.routed_P.channels.size());
    for (const PairChannelCount& entry : routes.routed_P.channels) {
        channels.push_back(entry.channel);
    }

    require_route_structure(routes.routed_P, routes.routed_PR, "PR");
    require_route_structure(routes.routed_P, routes.routed_PRD, "PRD");

    const PairSliceCountSummary expected =
        make_zero_pair_slice_count_summary(
            slices.pair_slicing,
            slices.origin_mode,
            channels
        );

    if (slices.entries.size() != expected.entries.size()) {
        throw std::logic_error(
            "pair processing: slice entry-count mismatch"
        );
    }

    std::vector<std::uint64_t> sum_P(channels.size(), 0U);
    std::vector<std::uint64_t> sum_PR(channels.size(), 0U);
    std::vector<std::uint64_t> sum_PRD(channels.size(), 0U);

    for (std::size_t slice_index = 0U;
         slice_index < slices.entries.size();
         ++slice_index) {
        const PairSliceCountEntry& entry = slices.entries[slice_index];
        const PairSliceCountEntry& expected_entry =
            expected.entries[slice_index];

        if (entry.kt_slice_index != expected_entry.kt_slice_index ||
            entry.mt_slice_index != expected_entry.mt_slice_index) {
            throw std::logic_error(
                "pair processing: slice entry identity mismatch"
            );
        }

        require_route_structure(
            routes.routed_P,
            entry.origin_counts.routed_P,
            "slice P"
        );
        require_route_structure(
            routes.routed_PR,
            entry.origin_counts.routed_PR,
            "slice PR"
        );
        require_route_structure(
            routes.routed_PRD,
            entry.origin_counts.routed_PRD,
            "slice PRD"
        );

        for (std::size_t channel_index = 0U;
             channel_index < channels.size();
             ++channel_index) {
            const std::uint64_t count_P =
                entry.origin_counts.routed_P.channels[channel_index]
                    .pair_count;
            const std::uint64_t count_PR =
                entry.origin_counts.routed_PR.channels[channel_index]
                    .pair_count;
            const std::uint64_t count_PRD =
                entry.origin_counts.routed_PRD.channels[channel_index]
                    .pair_count;

            switch (slices.origin_mode) {
            case OriginMode::Primordial:
                if (count_PR != 0U || count_PRD != 0U) {
                    throw std::logic_error(
                        "pair processing: invalid primordial slice routes"
                    );
                }
                break;
            case OriginMode::PrimordialRescattering:
                if (count_P != 0U || count_PRD != 0U) {
                    throw std::logic_error(
                        "pair processing: invalid PR slice routes"
                    );
                }
                break;
            case OriginMode::PrimordialRescatteringDecay:
                if (count_P != 0U || count_PR != 0U) {
                    throw std::logic_error(
                        "pair processing: invalid PRD slice routes"
                    );
                }
                break;
            case OriginMode::All:
                if (count_P > count_PR || count_PR > count_PRD) {
                    throw std::logic_error(
                        "pair processing: non-nested All slice routes"
                    );
                }
                break;
            default:
                throw std::invalid_argument(
                    "pair processing: invalid slice origin mode"
                );
            }

            const std::uint64_t values[] = {
                count_P,
                count_PR,
                count_PRD
            };
            std::vector<std::uint64_t>* sums[] = {
                &sum_P,
                &sum_PR,
                &sum_PRD
            };

            for (std::size_t route_index = 0U;
                 route_index < 3U;
                 ++route_index) {
                std::uint64_t& sum =
                    (*sums[route_index])[channel_index];
                if (values[route_index] >
                    std::numeric_limits<std::uint64_t>::max() - sum) {
                    throw std::logic_error(
                        "pair processing: slice route sum overflow"
                    );
                }
                sum += values[route_index];
            }
        }
    }

    for (std::size_t index = 0U; index < channels.size(); ++index) {
        if (sum_P[index] > routes.routed_P.channels[index].pair_count ||
            sum_PR[index] > routes.routed_PR.channels[index].pair_count ||
            sum_PRD[index] > routes.routed_PRD.channels[index].pair_count) {
            throw std::logic_error(
                "pair processing: slice counts exceed origin routes"
            );
        }
    }
}

/**
 * @brief Return whether at least one requested origin route is active.
 * @param routes Requested routes for one valid physical pair.
 * @return true when at least one route is active.
 */
bool has_requested_origin_route(const PairOriginRoutes& routes) noexcept {
    return routes.primordial ||
           routes.primordial_rescattering ||
           routes.primordial_rescattering_decay;
}

/**
 * @brief Return whether at least one kinetic slicing axis is active.
 * @param slicing Validated pair-slicing configuration.
 * @return true when kT or mT slicing is enabled.
 */
bool has_active_pair_slicing(
    const PairSlicingConfig& slicing
) noexcept {
    return slicing.kt.enabled || slicing.mt.enabled;
}

/**
 * @brief Increment all requested route counters for one valid pair.
 * @param counts Local origin-route counts to update.
 * @param index Channel index in the ordered route summaries.
 * @param routes Requested routes activated for the pair.
 * @throws std::overflow_error If any route counter would overflow.
 */
void increment_origin_route_counts(
    PairOriginRouteCountSummary& counts,
    std::size_t index,
    const PairOriginRoutes& routes
) {
    if (routes.primordial) {
        increment_pair_count(counts.routed_P.channels[index].pair_count);
    }
    if (routes.primordial_rescattering) {
        increment_pair_count(counts.routed_PR.channels[index].pair_count);
    }
    if (routes.primordial_rescattering_decay) {
        increment_pair_count(counts.routed_PRD.channels[index].pair_count);
    }
}

/**
 * @brief Add one complete numerical pair-rejection record.
 * @param report Local report to update.
 * @param outer_event_number One-based outer-event number.
 * @param subevent_id Current subevent identifier.
 * @param channel Primitive channel owning the pair.
 * @param ordinal One-based pair ordinal within the channel traversal.
 * @param particle_a Particle in canonical role A.
 * @param particle_b Particle in canonical role B.
 * @param kinematics Already calculated pair kinematics.
 * @param reason Exact numerical rejection reason.
 */
void record_pair_rejection(
    RejectedPairReport& report,
    std::size_t outer_event_number,
    int subevent_id,
    PrimitiveChannelId channel,
    std::uint64_t ordinal,
    const Particle& particle_a,
    const Particle& particle_b,
    const PairKinematics& kinematics,
    PairRejectionReason reason
) {
    report.add({
        outer_event_number,
        subevent_id,
        channel,
        ordinal,
        make_rejected_pair_particle_snapshot(particle_a),
        make_rejected_pair_particle_snapshot(particle_b),
        kinematics,
        reason
    });
}

/**
 * @brief Map a frame-observable numerical reason to report terminology.
 * @param reason Validation reason returned after frame calculation.
 * @return Matching stable pair-rejection reason.
 * @throws std::invalid_argument If @p reason is invalid.
 */
PairRejectionReason report_reason(
    PairFrameObservableNumericalReason reason
) {
    switch (reason) {
        case PairFrameObservableNumericalReason::NonFiniteDeltaTLab:
            return PairRejectionReason::NonFiniteDeltaTLab;
        case PairFrameObservableNumericalReason::NonFiniteDeltaTLcms:
            return PairRejectionReason::NonFiniteDeltaTLcms;
        case PairFrameObservableNumericalReason::NonFiniteDeltaTPrf:
            return PairRejectionReason::NonFiniteDeltaTPrf;
        case PairFrameObservableNumericalReason::NonFiniteROutLcms:
            return PairRejectionReason::NonFiniteROutLcms;
        case PairFrameObservableNumericalReason::NonFiniteROutPrf:
            return PairRejectionReason::NonFiniteROutPrf;
        case PairFrameObservableNumericalReason::NonFiniteRSide:
            return PairRejectionReason::NonFiniteRSide;
        case PairFrameObservableNumericalReason::NonFiniteRLong:
            return PairRejectionReason::NonFiniteRLong;
        case PairFrameObservableNumericalReason::NonFiniteRRadialLcms:
            return PairRejectionReason::NonFiniteRRadialLcms;
        case PairFrameObservableNumericalReason::NonFiniteRRadialPrf:
            return PairRejectionReason::NonFiniteRRadialPrf;
    }
    throw std::invalid_argument(
        "pair processing: invalid frame-observable rejection reason"
    );
}

}  // namespace

PairSubeventProcessingResult process_subevent_pairs(
    std::size_t outer_event_number,
    int subevent_id,
    const EventBuffers& buffers,
    const std::vector<PrimitiveChannelId>& required_channels,
    OriginMode origin_mode,
    const PairSlicingConfig& pair_slicing,
    PairFrameConsumer& frame_consumer
) {
    PairCountSummary valid_counts =
        make_zero_pair_count_summary(required_channels);
    PairCountSummary rejected_counts =
        make_zero_pair_count_summary(required_channels);
    PairCountSummary ordinals =
        make_zero_pair_count_summary(required_channels);
    PairOriginRouteCountSummary origin_route_counts{
        origin_mode,
        make_zero_pair_count_summary(required_channels),
        make_zero_pair_count_summary(required_channels),
        make_zero_pair_count_summary(required_channels)
    };
    PairSliceCountSummary pair_slice_counts =
        make_zero_pair_slice_count_summary(
            pair_slicing,
            origin_mode,
            required_channels
        );
    RejectedPairReport rejections;

    PairCountSummary formed_counts = for_each_pair(
        buffers,
        required_channels,
        [&](std::size_t channel_index,
            PrimitiveChannelId channel,
            const Particle& particle_a,
            const Particle& particle_b) {
            const std::size_t index = channel_index;
            PairChannelCount& ordinal_entry = ordinals.channels[index];
            increment_pair_count(ordinal_entry.pair_count);

            const PairKinematics kinematics =
                calculate_pair_kinematics(particle_a, particle_b);

            if (!is_finite_pair_kt(kinematics)) {
                increment_pair_count(
                    rejected_counts.channels[index].pair_count
                );
                record_pair_rejection(
                    rejections,
                    outer_event_number,
                    subevent_id,
                    channel,
                    ordinal_entry.pair_count,
                    particle_a,
                    particle_b,
                    kinematics,
                    PairRejectionReason::NonFiniteKt
                );
                return;
            }

            if (!is_finite_pair_mt(kinematics)) {
                increment_pair_count(
                    rejected_counts.channels[index].pair_count
                );
                record_pair_rejection(
                    rejections,
                    outer_event_number,
                    subevent_id,
                    channel,
                    ordinal_entry.pair_count,
                    particle_a,
                    particle_b,
                    kinematics,
                    PairRejectionReason::NonFiniteMt
                );
                return;
            }

            const PairOriginMemberships memberships =
                calculate_pair_origin_memberships(
                    particle_a.origin,
                    particle_b.origin
                );
            const PairOriginRoutes routes =
                route_pair_origin_memberships(memberships, origin_mode);

            if (!has_requested_origin_route(routes)) {
                throw std::logic_error(
                    "pair processing: valid pair has no requested origin route"
                );
            }

            const std::optional<PairSliceRoute> slice_route =
                route_pair_to_slices(kinematics, pair_slicing);
            const bool frame_required =
                !has_active_pair_slicing(pair_slicing) ||
                slice_route.has_value();
            std::optional<PairFrameObservables> observables;

            if (frame_required) {
                observables = calculate_pair_frame_observables(
                    particle_a,
                    particle_b,
                    kinematics
                );
                const auto numerical_reason =
                    first_non_finite_pair_frame_observable(
                        observables.value()
                    );
                if (numerical_reason.has_value()) {
                    increment_pair_count(
                        rejected_counts.channels[index].pair_count
                    );
                    record_pair_rejection(
                        rejections,
                        outer_event_number,
                        subevent_id,
                        channel,
                        ordinal_entry.pair_count,
                        particle_a,
                        particle_b,
                        kinematics,
                        report_reason(numerical_reason.value())
                    );
                    return;
                }
            }

            increment_pair_count(valid_counts.channels[index].pair_count);
            increment_origin_route_counts(
                origin_route_counts,
                index,
                routes
            );

            if (slice_route.has_value()) {
                increment_pair_slice_count(
                    pair_slice_counts,
                    slice_route.value(),
                    routes,
                    index,
                    channel
                );
            }

            if (frame_required) {
                const PairSliceRoute* const route =
                    slice_route.has_value() ? &slice_route.value() : nullptr;
                const PairFrameRouteContext context{
                    index,
                    channel,
                    routes,
                    route
                };
                frame_consumer.consume(
                    context,
                    kinematics,
                    observables.value()
                );
            }
        }
    );

    require_pair_count_partition(
        formed_counts,
        valid_counts,
        rejected_counts
    );

    require_rejection_report_matches_counts(
        rejected_counts,
        rejections
    );
    require_origin_route_contract(valid_counts, origin_route_counts);
    require_pair_slice_count_contract(
        origin_route_counts,
        pair_slice_counts
    );

    return {
        {
            outer_event_number,
            subevent_id,
            std::move(formed_counts),
            std::move(valid_counts),
            std::move(rejected_counts),
            std::move(origin_route_counts),
            std::move(pair_slice_counts)
        },
        std::move(rejections)
    };
}

void accumulate_pair_processing_result(
    HBTPairProcessingSummary& total,
    PairSubeventProcessingResult local
) {
    require_pair_count_partition(
        local.summary.pair_counts,
        local.summary.valid_pair_counts,
        local.summary.numerical_rejection_counts
    );
    require_local_rejection_identity(
        local.summary,
        local.numerical_rejections
    );
    require_rejection_report_matches_counts(
        local.summary.numerical_rejection_counts,
        local.numerical_rejections
    );
    require_origin_route_contract(
        local.summary.valid_pair_counts,
        local.summary.origin_route_counts
    );
    require_pair_slice_count_contract(
        local.summary.origin_route_counts,
        local.summary.pair_slice_counts
    );
    require_pair_count_partition(
        total.total_pair_counts,
        total.total_valid_pair_counts,
        total.total_numerical_rejection_counts
    );
    require_rejection_report_matches_counts(
        total.total_numerical_rejection_counts,
        total.numerical_rejections
    );
    require_origin_route_contract(
        total.total_valid_pair_counts,
        total.total_origin_route_counts
    );
    require_pair_slice_count_contract(
        total.total_origin_route_counts,
        total.total_pair_slice_counts
    );

    if (local.summary.origin_route_counts.origin_mode !=
        total.total_origin_route_counts.origin_mode) {
        throw std::invalid_argument(
            "pair processing: local/run-total origin mode mismatch"
        );
    }

    PairCountSummary formed = total.total_pair_counts;
    PairCountSummary valid = total.total_valid_pair_counts;
    PairCountSummary rejected = total.total_numerical_rejection_counts;
    PairOriginRouteCountSummary origin_routes =
        total.total_origin_route_counts;
    PairSliceCountSummary slice_counts =
        total.total_pair_slice_counts;

    accumulate_pair_counts(formed, local.summary.pair_counts);
    accumulate_pair_counts(valid, local.summary.valid_pair_counts);
    accumulate_pair_counts(
        rejected,
        local.summary.numerical_rejection_counts
    );
    accumulate_pair_counts(
        origin_routes.routed_P,
        local.summary.origin_route_counts.routed_P
    );
    accumulate_pair_counts(
        origin_routes.routed_PR,
        local.summary.origin_route_counts.routed_PR
    );
    accumulate_pair_counts(
        origin_routes.routed_PRD,
        local.summary.origin_route_counts.routed_PRD
    );
    accumulate_pair_slice_counts(
        slice_counts,
        local.summary.pair_slice_counts
    );

    require_pair_count_partition(formed, valid, rejected);
    require_origin_route_contract(valid, origin_routes);
    require_pair_slice_count_contract(origin_routes, slice_counts);

    total.total_pair_counts = std::move(formed);
    total.total_valid_pair_counts = std::move(valid);
    total.total_numerical_rejection_counts = std::move(rejected);
    total.total_origin_route_counts = std::move(origin_routes);
    total.total_pair_slice_counts = std::move(slice_counts);

    for (const RejectedPairRecord& record :
         local.numerical_rejections.records()) {
        total.numerical_rejections.add(record);
    }

    require_rejection_report_matches_counts(
        total.total_numerical_rejection_counts,
        total.numerical_rejections
    );

    total.subevents.push_back(std::move(local.summary));
}

void accumulate_pair_processing_summary(
    HBTPairProcessingSummary& total,
    HBTPairProcessingSummary local
) {
    require_pair_count_partition(
        local.total_pair_counts,
        local.total_valid_pair_counts,
        local.total_numerical_rejection_counts
    );
    require_rejection_report_matches_counts(
        local.total_numerical_rejection_counts,
        local.numerical_rejections
    );
    require_origin_route_contract(
        local.total_valid_pair_counts,
        local.total_origin_route_counts
    );
    require_pair_slice_count_contract(
        local.total_origin_route_counts,
        local.total_pair_slice_counts
    );
    require_pair_count_partition(
        total.total_pair_counts,
        total.total_valid_pair_counts,
        total.total_numerical_rejection_counts
    );
    require_rejection_report_matches_counts(
        total.total_numerical_rejection_counts,
        total.numerical_rejections
    );
    require_origin_route_contract(
        total.total_valid_pair_counts,
        total.total_origin_route_counts
    );
    require_pair_slice_count_contract(
        total.total_origin_route_counts,
        total.total_pair_slice_counts
    );

    if (local.total_origin_route_counts.origin_mode !=
        total.total_origin_route_counts.origin_mode) {
        throw std::invalid_argument(
            "pair processing: event/run-total origin mode mismatch"
        );
    }

    PairCountSummary formed = total.total_pair_counts;
    PairCountSummary valid = total.total_valid_pair_counts;
    PairCountSummary rejected = total.total_numerical_rejection_counts;
    PairOriginRouteCountSummary origin_routes =
        total.total_origin_route_counts;
    PairSliceCountSummary slice_counts = total.total_pair_slice_counts;

    accumulate_pair_counts(formed, local.total_pair_counts);
    accumulate_pair_counts(valid, local.total_valid_pair_counts);
    accumulate_pair_counts(
        rejected,
        local.total_numerical_rejection_counts
    );
    accumulate_pair_counts(
        origin_routes.routed_P,
        local.total_origin_route_counts.routed_P
    );
    accumulate_pair_counts(
        origin_routes.routed_PR,
        local.total_origin_route_counts.routed_PR
    );
    accumulate_pair_counts(
        origin_routes.routed_PRD,
        local.total_origin_route_counts.routed_PRD
    );
    accumulate_pair_slice_counts(
        slice_counts,
        local.total_pair_slice_counts
    );

    require_pair_count_partition(formed, valid, rejected);
    require_origin_route_contract(valid, origin_routes);
    require_pair_slice_count_contract(origin_routes, slice_counts);

    total.total_pair_counts = std::move(formed);
    total.total_valid_pair_counts = std::move(valid);
    total.total_numerical_rejection_counts = std::move(rejected);
    total.total_origin_route_counts = std::move(origin_routes);
    total.total_pair_slice_counts = std::move(slice_counts);

    for (const RejectedPairRecord& record :
         local.numerical_rejections.records()) {
        total.numerical_rejections.add(record);
    }
    for (HBTPairSubeventSummary& subevent : local.subevents) {
        total.subevents.push_back(std::move(subevent));
    }

    require_rejection_report_matches_counts(
        total.total_numerical_rejection_counts,
        total.numerical_rejections
    );
}

}  // namespace hbt
