/**
 * @file pair_slice_count_accumulator_test.cpp
 * @brief Unit tests for pair counts by kinetic slice, origin, and channel.
 */

#include "hbt/pair/pair_slice_count_accumulator.h"

#include <cstdint>
#include <cstdlib>
#include <iostream>
#include <limits>
#include <stdexcept>
#include <utility>
#include <vector>

namespace {

/**
 * @brief Build one slicing axis for concise count-test setup.
 * @param enabled Whether the axis participates in slicing.
 * @param edges Configured bin edges.
 * @return Axis configuration containing the supplied values.
 */
hbt::PairSlicingAxisConfig axis(
    bool enabled,
    std::vector<double> edges
) {
    return {enabled, std::move(edges)};
}

/**
 * @brief Return the two-channel test ordering used by all count checks.
 * @return Stable unique primitive-channel list.
 */
std::vector<hbt::PrimitiveChannelId> test_channels() {
    return {
        hbt::PrimitiveChannelId::PiPlusPiPlus,
        hbt::PrimitiveChannelId::PiPlusProton
    };
}

/**
 * @brief Read one count from one slice entry and selected origin route.
 * @param entry Slice entry containing nested origin counts.
 * @param route One-character route selector: P, R for PR, or D for PRD.
 * @param channel_index Ordered primitive-channel index.
 * @return Stored pair count.
 * @throws std::invalid_argument If @p route is unsupported.
 */
std::uint64_t count_for(
    const hbt::PairSliceCountEntry& entry,
    char route,
    std::size_t channel_index
) {
    switch (route) {
    case 'P':
        return entry.origin_counts.routed_P
            .channels[channel_index].pair_count;
    case 'R':
        return entry.origin_counts.routed_PR
            .channels[channel_index].pair_count;
    case 'D':
        return entry.origin_counts.routed_PRD
            .channels[channel_index].pair_count;
    }

    throw std::invalid_argument("invalid test route selector");
}

/**
 * @brief Verify disabled slicing creates no count entry or dummy slice.
 * @return `true` when the zero summary has no entries.
 */
bool verify_disabled_slicing_has_no_entries() {
    const hbt::PairSlicingConfig slicing{
        axis(false, {0.2, 0.4}),
        axis(false, {0.5, 0.8})
    };
    const hbt::PairSliceCountSummary summary =
        hbt::make_zero_pair_slice_count_summary(
            slicing,
            hbt::OriginMode::All,
            test_channels()
        );

    if (!summary.entries.empty()) {
        std::cerr
            << "pair_slice_count_accumulator_test: dummy slice created.\n";
        return false;
    }

    return true;
}

/**
 * @brief Verify deterministic kT-major Cartesian slice layout.
 * @return `true` when all four expected cells are present in exact order.
 */
bool verify_cartesian_layout() {
    const hbt::PairSlicingConfig slicing{
        axis(true, {0.2, 0.4, 0.6}),
        axis(true, {0.5, 0.8, 1.1})
    };
    const hbt::PairSliceCountSummary summary =
        hbt::make_zero_pair_slice_count_summary(
            slicing,
            hbt::OriginMode::All,
            test_channels()
        );

    if (summary.entries.size() != 4U) {
        std::cerr
            << "pair_slice_count_accumulator_test: wrong cell count.\n";
        return false;
    }

    const std::size_t expected[][2] = {
        {0U, 0U},
        {0U, 1U},
        {1U, 0U},
        {1U, 1U}
    };

    for (std::size_t index = 0U; index < summary.entries.size(); ++index) {
        const auto& entry = summary.entries[index];
        if (entry.kt_slice_index != expected[index][0] ||
            entry.mt_slice_index != expected[index][1]) {
            std::cerr
                << "pair_slice_count_accumulator_test: wrong cell order.\n";
            return false;
        }
    }

    return true;
}

/**
 * @brief Verify one resolved slice is reused for nested All-origin counts.
 * @return `true` when P-P, P-R, and D-like routes update only one cell.
 */
bool verify_all_mode_counts_one_resolved_cell() {
    const hbt::PairSlicingConfig slicing{
        axis(true, {0.2, 0.4, 0.6}),
        axis(true, {0.5, 0.8, 1.1})
    };
    hbt::PairSliceCountSummary summary =
        hbt::make_zero_pair_slice_count_summary(
            slicing,
            hbt::OriginMode::All,
            test_channels()
        );
    const hbt::PairSliceRoute route{1U, 0U, 2U};
    constexpr std::size_t channel_index = 1U;
    constexpr hbt::PrimitiveChannelId channel =
        hbt::PrimitiveChannelId::PiPlusProton;

    hbt::increment_pair_slice_count(
        summary,
        route,
        {true, true, true},
        channel_index,
        channel
    );
    hbt::increment_pair_slice_count(
        summary,
        route,
        {false, true, true},
        channel_index,
        channel
    );
    hbt::increment_pair_slice_count(
        summary,
        route,
        {false, false, true},
        channel_index,
        channel
    );

    const hbt::PairSliceCountEntry& target = summary.entries[2U];
    if (count_for(target, 'P', channel_index) != 1U ||
        count_for(target, 'R', channel_index) != 2U ||
        count_for(target, 'D', channel_index) != 3U) {
        std::cerr
            << "pair_slice_count_accumulator_test: wrong All counts.\n";
        return false;
    }

    for (std::size_t index = 0U; index < summary.entries.size(); ++index) {
        if (index == 2U) {
            continue;
        }
        if (count_for(summary.entries[index], 'P', channel_index) != 0U ||
            count_for(summary.entries[index], 'R', channel_index) != 0U ||
            count_for(summary.entries[index], 'D', channel_index) != 0U) {
            std::cerr
                << "pair_slice_count_accumulator_test: other cell changed.\n";
            return false;
        }
    }

    return true;
}

/**
 * @brief Verify an individual mode increments only its requested origin route.
 * @return `true` when Primordial updates P and leaves PR/PRD untouched.
 */
bool verify_individual_mode_counts_only_requested_route() {
    const hbt::PairSlicingConfig slicing{
        axis(true, {0.2, 0.4}),
        axis(false, {})
    };
    hbt::PairSliceCountSummary summary =
        hbt::make_zero_pair_slice_count_summary(
            slicing,
            hbt::OriginMode::Primordial,
            test_channels()
        );

    hbt::increment_pair_slice_count(
        summary,
        {0U, std::nullopt, 0U},
        {true, false, false},
        0U,
        hbt::PrimitiveChannelId::PiPlusPiPlus
    );

    const hbt::PairSliceCountEntry& entry = summary.entries[0U];
    if (count_for(entry, 'P', 0U) != 1U ||
        count_for(entry, 'R', 0U) != 0U ||
        count_for(entry, 'D', 0U) != 0U) {
        std::cerr
            << "pair_slice_count_accumulator_test: individual mode leaked.\n";
        return false;
    }

    return true;
}

/**
 * @brief Verify contradictory origin routes are structural errors.
 * @return `true` when Primordial mode refuses a PR route.
 */
bool verify_mode_route_mismatch_is_rejected() {
    const hbt::PairSlicingConfig slicing{
        axis(true, {0.2, 0.4}),
        axis(false, {})
    };
    hbt::PairSliceCountSummary summary =
        hbt::make_zero_pair_slice_count_summary(
            slicing,
            hbt::OriginMode::Primordial,
            test_channels()
        );

    try {
        hbt::increment_pair_slice_count(
            summary,
            {0U, std::nullopt, 0U},
            {false, true, false},
            0U,
            hbt::PrimitiveChannelId::PiPlusPiPlus
        );
    } catch (const std::invalid_argument&) {
        return true;
    }

    std::cerr
        << "pair_slice_count_accumulator_test: mode mismatch accepted.\n";
    return false;
}

/**
 * @brief Verify a flat index must identify the same routed slice indices.
 * @return `true` when an inconsistent precomputed flat index is rejected.
 */
bool verify_flat_index_identity_mismatch_is_rejected() {
    const hbt::PairSlicingConfig slicing{
        axis(true, {0.2, 0.4, 0.6}),
        axis(true, {0.5, 0.8, 1.1})
    };
    hbt::PairSliceCountSummary summary =
        hbt::make_zero_pair_slice_count_summary(
            slicing,
            hbt::OriginMode::All,
            test_channels()
        );

    try {
        hbt::increment_pair_slice_count(
            summary,
            {1U, 0U, 0U},
            {true, true, true},
            0U,
            hbt::PrimitiveChannelId::PiPlusPiPlus
        );
    } catch (const std::invalid_argument&) {
        return true;
    }

    std::cerr
        << "pair_slice_count_accumulator_test: bad flat index accepted.\n";
    return false;
}

/**
 * @brief Verify channel index and channel identity must agree.
 * @return `true` when a mismatched channel is rejected.
 */
bool verify_channel_identity_mismatch_is_rejected() {
    const hbt::PairSlicingConfig slicing{
        axis(true, {0.2, 0.4}),
        axis(false, {})
    };
    hbt::PairSliceCountSummary summary =
        hbt::make_zero_pair_slice_count_summary(
            slicing,
            hbt::OriginMode::All,
            test_channels()
        );

    try {
        hbt::increment_pair_slice_count(
            summary,
            {0U, std::nullopt, 0U},
            {true, true, true},
            0U,
            hbt::PrimitiveChannelId::PiPlusProton
        );
    } catch (const std::invalid_argument&) {
        return true;
    }

    std::cerr
        << "pair_slice_count_accumulator_test: channel mismatch accepted.\n";
    return false;
}

/**
 * @brief Verify overflow is detected before any nested route mutates.
 * @return `true` when a failed P-P-like increment is atomic.
 */
bool verify_increment_overflow_is_atomic() {
    const hbt::PairSlicingConfig slicing{
        axis(true, {0.2, 0.4}),
        axis(false, {})
    };
    hbt::PairSliceCountSummary summary =
        hbt::make_zero_pair_slice_count_summary(
            slicing,
            hbt::OriginMode::All,
            test_channels()
        );
    hbt::PairSliceCountEntry& entry = summary.entries[0U];
    entry.origin_counts.routed_P.channels[0U].pair_count =
        std::numeric_limits<std::uint64_t>::max();

    try {
        hbt::increment_pair_slice_count(
            summary,
            {0U, std::nullopt, 0U},
            {true, true, true},
            0U,
            hbt::PrimitiveChannelId::PiPlusPiPlus
        );
    } catch (const std::overflow_error&) {
        if (count_for(entry, 'R', 0U) == 0U &&
            count_for(entry, 'D', 0U) == 0U) {
            return true;
        }
    }

    std::cerr
        << "pair_slice_count_accumulator_test: overflow was not atomic.\n";
    return false;
}

/**
 * @brief Verify compatible local slice counts accumulate into run totals.
 * @return `true` when all nested route counts are added exactly.
 */
bool verify_summary_accumulation() {
    const hbt::PairSlicingConfig slicing{
        axis(false, {}),
        axis(true, {0.5, 0.8, 1.1})
    };
    hbt::PairSliceCountSummary total =
        hbt::make_zero_pair_slice_count_summary(
            slicing,
            hbt::OriginMode::All,
            test_channels()
        );
    hbt::PairSliceCountSummary local =
        hbt::make_zero_pair_slice_count_summary(
            slicing,
            hbt::OriginMode::All,
            test_channels()
        );

    hbt::increment_pair_slice_count(
        local,
        {std::nullopt, 1U, 1U},
        {false, true, true},
        1U,
        hbt::PrimitiveChannelId::PiPlusProton
    );
    hbt::accumulate_pair_slice_counts(total, local);

    const hbt::PairSliceCountEntry& entry = total.entries[1U];
    if (count_for(entry, 'P', 1U) != 0U ||
        count_for(entry, 'R', 1U) != 1U ||
        count_for(entry, 'D', 1U) != 1U) {
        std::cerr
            << "pair_slice_count_accumulator_test: bad accumulation.\n";
        return false;
    }

    return true;
}

/**
 * @brief Verify incompatible slice configurations cannot be accumulated.
 * @return `true` when a different edge sequence is rejected.
 */
bool verify_accumulation_mismatch_is_rejected() {
    hbt::PairSliceCountSummary total =
        hbt::make_zero_pair_slice_count_summary(
            {axis(true, {0.2, 0.4}), axis(false, {})},
            hbt::OriginMode::All,
            test_channels()
        );
    const hbt::PairSliceCountSummary local =
        hbt::make_zero_pair_slice_count_summary(
            {axis(true, {0.2, 0.5}), axis(false, {})},
            hbt::OriginMode::All,
            test_channels()
        );

    try {
        hbt::accumulate_pair_slice_counts(total, local);
    } catch (const std::invalid_argument&) {
        return true;
    }

    std::cerr
        << "pair_slice_count_accumulator_test: config mismatch accepted.\n";
    return false;
}

}  // namespace

/**
 * @brief Run all slice/origin/channel count-accumulator unit tests.
 * @return `EXIT_SUCCESS` when every test passes, otherwise `EXIT_FAILURE`.
 */
int main() {
    if (!verify_disabled_slicing_has_no_entries() ||
        !verify_cartesian_layout() ||
        !verify_all_mode_counts_one_resolved_cell() ||
        !verify_individual_mode_counts_only_requested_route() ||
        !verify_mode_route_mismatch_is_rejected() ||
        !verify_flat_index_identity_mismatch_is_rejected() ||
        !verify_channel_identity_mismatch_is_rejected() ||
        !verify_increment_overflow_is_atomic() ||
        !verify_summary_accumulation() ||
        !verify_accumulation_mismatch_is_rejected()) {
        return EXIT_FAILURE;
    }

    return EXIT_SUCCESS;
}
