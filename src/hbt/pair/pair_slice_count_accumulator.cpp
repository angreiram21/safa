/**
 * @file pair_slice_count_accumulator.cpp
 * @brief Pair slice/origin/channel count accumulation implementation.
 */

#include "hbt/pair/pair_slice_count_accumulator.h"

#include "hbt/pair/pair_count_accumulator.h"

#include <cstddef>
#include <cstdint>
#include <limits>
#include <stdexcept>
#include <utility>

namespace hbt {
namespace {

/**
 * @brief Validate that an OriginMode enumerator is supported.
 * @param mode Origin mode to validate.
 * @throws std::invalid_argument If the enumerator is invalid.
 */
void validate_origin_mode(OriginMode mode) {
    switch (mode) {
    case OriginMode::Primordial:
    case OriginMode::PrimordialRescattering:
    case OriginMode::PrimordialRescatteringDecay:
    case OriginMode::All:
        return;
    }

    throw std::invalid_argument(
        "pair slice counts: invalid OriginMode"
    );
}

/**
 * @brief Return the number of configured slices on one axis.
 * @param axis Validated slicing-axis configuration.
 * @return Zero for a disabled axis, otherwise edge count minus one.
 * @throws std::invalid_argument If an enabled axis has fewer than two edges.
 */
std::size_t axis_slice_count(const PairSlicingAxisConfig& axis) {
    if (!axis.enabled) {
        return 0U;
    }

    if (axis.bin_edges_gev.size() < 2U) {
        throw std::invalid_argument(
            "pair slice counts: enabled axis has fewer than two edges"
        );
    }

    return axis.bin_edges_gev.size() - 1U;
}

/**
 * @brief Compare two slicing axes exactly.
 * @param lhs First axis.
 * @param rhs Second axis.
 * @return `true` when enabled state and retained edges are identical.
 */
bool same_axis(
    const PairSlicingAxisConfig& lhs,
    const PairSlicingAxisConfig& rhs
) {
    return lhs.enabled == rhs.enabled &&
           lhs.bin_edges_gev == rhs.bin_edges_gev;
}

/**
 * @brief Compare two complete slicing configurations exactly.
 * @param lhs First configuration.
 * @param rhs Second configuration.
 * @return `true` when both kT and mT axis definitions match.
 */
bool same_slicing(
    const PairSlicingConfig& lhs,
    const PairSlicingConfig& rhs
) {
    return same_axis(lhs.kt, rhs.kt) && same_axis(lhs.mt, rhs.mt);
}

/**
 * @brief Build zeroed nested-origin counts for one slice.
 * @param channels Ordered primitive channels.
 * @return Three zeroed per-channel summaries.
 */
PairSliceOriginCounts make_zero_origin_counts(
    const std::vector<PrimitiveChannelId>& channels
) {
    return {
        make_zero_pair_count_summary(channels),
        make_zero_pair_count_summary(channels),
        make_zero_pair_count_summary(channels)
    };
}

/**
 * @brief Validate requested origin-route booleans against one OriginMode.
 * @param routes Requested compatible routes for one pair.
 * @param mode Origin mode that owns the count summary.
 * @throws std::invalid_argument If routes cannot be produced by @p mode.
 */
void validate_routes_for_mode(
    const PairOriginRoutes& routes,
    OriginMode mode
) {
    validate_origin_mode(mode);

    switch (mode) {
    case OriginMode::Primordial:
        if (routes.primordial_rescattering ||
            routes.primordial_rescattering_decay) {
            throw std::invalid_argument(
                "pair slice counts: routes contradict primordial mode"
            );
        }
        return;
    case OriginMode::PrimordialRescattering:
        if (routes.primordial ||
            routes.primordial_rescattering_decay) {
            throw std::invalid_argument(
                "pair slice counts: routes contradict PR mode"
            );
        }
        return;
    case OriginMode::PrimordialRescatteringDecay:
        if (routes.primordial || routes.primordial_rescattering) {
            throw std::invalid_argument(
                "pair slice counts: routes contradict PRD mode"
            );
        }
        return;
    case OriginMode::All:
        if ((routes.primordial &&
             !routes.primordial_rescattering) ||
            (routes.primordial_rescattering &&
             !routes.primordial_rescattering_decay)) {
            throw std::invalid_argument(
                "pair slice counts: non-nested routes in All mode"
            );
        }
        return;
    }

    throw std::invalid_argument(
        "pair slice counts: invalid OriginMode"
    );
}

/**
 * @brief Validate and reuse the flat entry index carried by one slice route.
 * @param summary Slice-count summary owning the configured layout.
 * @param route Unique already resolved kinetic route.
 * @return The route's already calculated flat slice index.
 * @throws std::invalid_argument If route identity contradicts the summary.
 */
std::size_t entry_index_for_route(
    const PairSliceCountSummary& summary,
    const PairSliceRoute& route
) {
    const bool kt_enabled = summary.pair_slicing.kt.enabled;
    const bool mt_enabled = summary.pair_slicing.mt.enabled;

    if (!kt_enabled && !mt_enabled) {
        throw std::invalid_argument(
            "pair slice counts: no slicing axis is enabled"
        );
    }

    if (route.kt_slice_index.has_value() != kt_enabled ||
        route.mt_slice_index.has_value() != mt_enabled) {
        throw std::invalid_argument(
            "pair slice counts: route axis presence mismatch"
        );
    }

    if (route.flat_slice_index >= summary.entries.size()) {
        throw std::invalid_argument(
            "pair slice counts: flat slice index out of range"
        );
    }

    const PairSliceCountEntry& entry =
        summary.entries[route.flat_slice_index];
    if (entry.kt_slice_index != route.kt_slice_index ||
        entry.mt_slice_index != route.mt_slice_index) {
        throw std::invalid_argument(
            "pair slice counts: route identity/flat index mismatch"
        );
    }

    return route.flat_slice_index;
}

/**
 * @brief Validate one entry's channel identity at a known ordered index.
 * @param counts Per-origin channel counts in one kinetic slice.
 * @param channel_index Expected ordered channel index.
 * @param channel Expected primitive channel identifier.
 * @throws std::invalid_argument If any origin route disagrees structurally.
 */
void validate_channel_identity(
    const PairSliceOriginCounts& counts,
    std::size_t channel_index,
    PrimitiveChannelId channel
) {
    const PairCountSummary* summaries[] = {
        &counts.routed_P,
        &counts.routed_PR,
        &counts.routed_PRD
    };

    for (const PairCountSummary* summary : summaries) {
        if (channel_index >= summary->channels.size() ||
            summary->channels[channel_index].channel != channel) {
            throw std::invalid_argument(
                "pair slice counts: channel index/identity mismatch"
            );
        }
    }
}

/**
 * @brief Throw if incrementing one count would overflow uint64_t.
 * @param count Count to inspect before mutation.
 * @throws std::overflow_error If @p count is already at its maximum.
 */
void require_increment_capacity(std::uint64_t count) {
    if (count == std::numeric_limits<std::uint64_t>::max()) {
        throw std::overflow_error(
            "pair slice counts: pair-count overflow"
        );
    }
}

/**
 * @brief Validate deterministic entry identity between two summaries.
 * @param total Existing total summary.
 * @param subevent Local summary to accumulate.
 * @throws std::invalid_argument If summary-level slice structure differs.
 */
void validate_summary_identity(
    const PairSliceCountSummary& total,
    const PairSliceCountSummary& subevent
) {
    validate_origin_mode(total.origin_mode);
    validate_origin_mode(subevent.origin_mode);

    if (total.origin_mode != subevent.origin_mode) {
        throw std::invalid_argument(
            "pair slice counts: origin-mode mismatch"
        );
    }
    if (!same_slicing(total.pair_slicing, subevent.pair_slicing)) {
        throw std::invalid_argument(
            "pair slice counts: slicing-config mismatch"
        );
    }
    if (total.entries.size() != subevent.entries.size()) {
        throw std::invalid_argument(
            "pair slice counts: slice-entry-count mismatch"
        );
    }

    for (std::size_t index = 0U; index < total.entries.size(); ++index) {
        if (total.entries[index].kt_slice_index !=
                subevent.entries[index].kt_slice_index ||
            total.entries[index].mt_slice_index !=
                subevent.entries[index].mt_slice_index) {
            throw std::invalid_argument(
                "pair slice counts: slice-entry identity mismatch"
            );
        }
    }
}

}  // namespace

PairSliceCountSummary make_zero_pair_slice_count_summary(
    const PairSlicingConfig& slicing,
    OriginMode origin_mode,
    const std::vector<PrimitiveChannelId>& channels
) {
    validate_origin_mode(origin_mode);

    const std::size_t kt_count = axis_slice_count(slicing.kt);
    const std::size_t mt_count = axis_slice_count(slicing.mt);

    PairSliceCountSummary summary{slicing, origin_mode, {}};

    if (!slicing.kt.enabled && !slicing.mt.enabled) {
        return summary;
    }

    if (slicing.kt.enabled && slicing.mt.enabled) {
        if (mt_count != 0U &&
            kt_count > std::numeric_limits<std::size_t>::max() / mt_count) {
            throw std::overflow_error(
                "pair slice counts: Cartesian slice-count overflow"
            );
        }

        summary.entries.reserve(kt_count * mt_count);
        for (std::size_t kt_index = 0U;
             kt_index < kt_count;
             ++kt_index) {
            for (std::size_t mt_index = 0U;
                 mt_index < mt_count;
                 ++mt_index) {
                summary.entries.push_back({
                    kt_index,
                    mt_index,
                    make_zero_origin_counts(channels)
                });
            }
        }
        return summary;
    }

    if (slicing.kt.enabled) {
        summary.entries.reserve(kt_count);
        for (std::size_t kt_index = 0U;
             kt_index < kt_count;
             ++kt_index) {
            summary.entries.push_back({
                kt_index,
                std::nullopt,
                make_zero_origin_counts(channels)
            });
        }
        return summary;
    }

    summary.entries.reserve(mt_count);
    for (std::size_t mt_index = 0U;
         mt_index < mt_count;
         ++mt_index) {
        summary.entries.push_back({
            std::nullopt,
            mt_index,
            make_zero_origin_counts(channels)
        });
    }

    return summary;
}

void increment_pair_slice_count(
    PairSliceCountSummary& summary,
    const PairSliceRoute& slice_route,
    const PairOriginRoutes& origin_routes,
    std::size_t channel_index,
    PrimitiveChannelId channel
) {
    validate_routes_for_mode(origin_routes, summary.origin_mode);

    const std::size_t entry_index =
        entry_index_for_route(summary, slice_route);
    if (entry_index >= summary.entries.size()) {
        throw std::invalid_argument(
            "pair slice counts: missing configured slice entry"
        );
    }

    PairSliceCountEntry& entry = summary.entries[entry_index];
    if (entry.kt_slice_index != slice_route.kt_slice_index ||
        entry.mt_slice_index != slice_route.mt_slice_index) {
        throw std::invalid_argument(
            "pair slice counts: stored slice layout mismatch"
        );
    }

    validate_channel_identity(
        entry.origin_counts,
        channel_index,
        channel
    );

    PairChannelCount& p =
        entry.origin_counts.routed_P.channels[channel_index];
    PairChannelCount& pr =
        entry.origin_counts.routed_PR.channels[channel_index];
    PairChannelCount& prd =
        entry.origin_counts.routed_PRD.channels[channel_index];

    if (origin_routes.primordial) {
        require_increment_capacity(p.pair_count);
    }
    if (origin_routes.primordial_rescattering) {
        require_increment_capacity(pr.pair_count);
    }
    if (origin_routes.primordial_rescattering_decay) {
        require_increment_capacity(prd.pair_count);
    }

    if (origin_routes.primordial) {
        ++p.pair_count;
    }
    if (origin_routes.primordial_rescattering) {
        ++pr.pair_count;
    }
    if (origin_routes.primordial_rescattering_decay) {
        ++prd.pair_count;
    }
}

void accumulate_pair_slice_counts(
    PairSliceCountSummary& total,
    const PairSliceCountSummary& subevent
) {
    validate_summary_identity(total, subevent);

    PairSliceCountSummary candidate = total;
    for (std::size_t index = 0U; index < candidate.entries.size(); ++index) {
        PairSliceOriginCounts& candidate_counts =
            candidate.entries[index].origin_counts;
        const PairSliceOriginCounts& local_counts =
            subevent.entries[index].origin_counts;

        accumulate_pair_counts(
            candidate_counts.routed_P,
            local_counts.routed_P
        );
        accumulate_pair_counts(
            candidate_counts.routed_PR,
            local_counts.routed_PR
        );
        accumulate_pair_counts(
            candidate_counts.routed_PRD,
            local_counts.routed_PRD
        );
    }

    total = std::move(candidate);
}

}  // namespace hbt
