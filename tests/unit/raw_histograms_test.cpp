/**
 * @file raw_histograms_test.cpp
 * @brief Unit tests for Phase-6 raw histogram allocation and accumulation.
 */

#include "hbt/histograms/raw_histograms.h"
#include "hbt/startup/hbt_startup_builder.h"

#include <cstddef>
#include <cstdint>
#include <cstdlib>
#include <iostream>
#include <limits>
#include <optional>
#include <stdexcept>
#include <utility>
#include <vector>

namespace hbt::test {

/**
 * @brief Reset the test-only raw-histogram classification counter.
 */
void reset_raw_histogram_bin_classification_count() noexcept;

/**
 * @brief Return the test-only raw-histogram classification count.
 * @return Number of bin classifications observed since the last reset.
 */
std::size_t raw_histogram_bin_classification_count() noexcept;

}  // namespace hbt::test

namespace {

/**
 * @brief Report one failed raw-histogram condition.
 * @param message Description of the violated contract.
 * @return Always false.
 */
bool fail(const char* message) {
    std::cerr << "raw_histograms_test: " << message << ".\n";
    return false;
}

/**
 * @brief Construct one complete compact histogram test configuration.
 * @param selection Final products to configure.
 * @param origin_mode Requested origin mode.
 * @param slicing Optional pair-slicing configuration.
 * @return Complete validated-style HBT configuration.
 */
hbt::HBTConfig make_config(
    hbt::HBTSelection selection,
    hbt::OriginMode origin_mode,
    hbt::PairSlicingConfig slicing
) {
    return {
        std::move(selection),
        {
            hbt::LongitudinalVariable::Pseudorapidity,
            {0.8, 0.14, 4.0},
            {0.8, 0.4, 1.4},
            {0.8, 0.5, 4.05},
            {0.8, 1.0, 10000.0},
            {0.8, 0.3, 10000.0}
        },
        std::move(slicing),
        {
            {4U, 0.0, 4.0, 1.0},
            {4U, 0.0, 4.0, 1.0},
            {4U, -2.0, 2.0, 1.0}
        },
        origin_mode,
        hbt::FitEstimatorMode::All
    };
}

/**
 * @brief Return two overlapping products sharing PiPlusPiPlus.
 * @return Product zero A+B and product one A+C.
 */
hbt::HBTSelection overlapping_selection() {
    return {{
        hbt::AnalysisProduct{{
            hbt::PrimitiveChannelId::PiPlusPiPlus,
            hbt::PrimitiveChannelId::PiMinusPiMinus
        }},
        hbt::AnalysisProduct{{
            hbt::PrimitiveChannelId::PiPlusPiPlus,
            hbt::PrimitiveChannelId::KPlusKPlus
        }}
    }};
}

/**
 * @brief Return frame observables whose nine bins are easy to identify.
 * @return Finite observables spanning all configured logical slots.
 */
hbt::PairFrameObservables in_range_observables() {
    return {
        -1.5,
        -0.5,
        0.5,
        -0.5,
        1.5,
        -2.5,
        3.5,
        0.5,
        1.5
    };
}

/**
 * @brief Return one dummy finite kinematics object for consumer tests.
 * @return Complete finite pair kinematics.
 */
hbt::PairKinematics finite_kinematics() {
    return {{4.0, 1.0, 0.0, 0.0}, 0.5, 1.0, 3.0, 1.0};
}

/**
 * @brief Return one flattened counter for a logical histogram slot and bin.
 * @param bins Slot-major flattened counter block.
 * @param nbins Number of bins owned by every logical slot.
 * @param slot Logical observable slot.
 * @param bin Zero-based bin within the logical slot.
 * @return Stored raw count.
 */
std::uint64_t count_at(
    const std::vector<std::uint64_t>& bins,
    std::size_t nbins,
    std::size_t slot,
    std::size_t bin
) {
    return bins[slot * nbins + bin];
}

/**
 * @brief Verify eager allocation matches product, origin, and slice layout.
 * @return true when all exact dimensions and zero state are correct.
 */
bool verify_eager_layout() {
    const hbt::HBTConfig config = make_config(
        overlapping_selection(),
        hbt::OriginMode::All,
        {{true, {0.0, 1.0, 2.0}}, {false, {}}}
    );
    const hbt::RawHistogramState state =
        hbt::make_zero_raw_histogram_state(config);

    if (state.products.size() != 2U) {
        return fail("one state was not allocated per final product");
    }
    for (const hbt::ProductRawHistogramState& product : state.products) {
        if (product.origins.size() != 3U) {
            return fail("OriginMode::All did not allocate P/PR/PRD");
        }
        for (const hbt::OriginRawHistogramState& origin : product.origins) {
            if (origin.slices.size() != 2U ||
                origin.global.osl.bins.size() != 16U ||
                origin.global.radial.bins.size() != 8U ||
                origin.global.delta_t.bins.size() != 12U) {
                return fail("raw histogram destination dimensions differ");
            }
            for (const std::uint64_t count : origin.global.osl.bins) {
                if (count != 0U) {
                    return fail("eagerly allocated OSL state is not zero");
                }
            }
        }
    }
    return true;
}

/**
 * @brief Verify one channel result fans out without intermediate histograms.
 * @return true when both final products, all origins, global, and one slice
 *         receive the same already classified observables exactly once.
 */
bool verify_shared_product_origin_and_slice_fanout() {
    const hbt::HBTConfig config = make_config(
        overlapping_selection(),
        hbt::OriginMode::All,
        {{true, {0.0, 1.0, 2.0}}, {false, {}}}
    );
    const hbt::HBTStartupState startup =
        hbt::build_hbt_startup_state(config);
    hbt::RawHistogramState state =
        hbt::make_zero_raw_histogram_state(config);
    hbt::RawHistogramAccumulator accumulator(config, startup, state);
    const hbt::PairSliceRoute slice{{0U}, std::nullopt, 0U};
    const hbt::PairFrameRouteContext context{
        0U,
        hbt::PrimitiveChannelId::PiPlusPiPlus,
        {true, true, true},
        &slice
    };

    hbt::test::reset_raw_histogram_bin_classification_count();
    accumulator.consume(
        context,
        finite_kinematics(),
        in_range_observables()
    );

    if (hbt::test::raw_histogram_bin_classification_count() != 9U) {
        return fail(
            "one consumed pair did not execute exactly nine classifications"
        );
    }

    for (const hbt::ProductRawHistogramState& product : state.products) {
        for (const hbt::OriginRawHistogramState& origin : product.origins) {
            if (count_at(origin.global.osl.bins, 4U, 0U, 0U) != 1U ||
                count_at(origin.global.osl.bins, 4U, 1U, 1U) != 1U ||
                count_at(origin.global.osl.bins, 4U, 2U, 2U) != 1U ||
                count_at(origin.global.osl.bins, 4U, 3U, 3U) != 1U ||
                count_at(origin.global.radial.bins, 4U, 0U, 0U) != 1U ||
                count_at(origin.global.radial.bins, 4U, 1U, 1U) != 1U ||
                count_at(origin.global.delta_t.bins, 4U, 0U, 0U) != 1U ||
                count_at(origin.global.delta_t.bins, 4U, 1U, 1U) != 1U ||
                count_at(origin.global.delta_t.bins, 4U, 2U, 2U) != 1U) {
                return fail("global fan-out counts are incorrect");
            }
            if (origin.slices[0].osl.bins != origin.global.osl.bins ||
                origin.slices[0].radial.bins != origin.global.radial.bins ||
                origin.slices[0].delta_t.bins != origin.global.delta_t.bins) {
                return fail("slice did not receive the global pair counts");
            }
            for (const std::uint64_t count : origin.slices[1].osl.bins) {
                if (count != 0U) {
                    return fail("unrouted slice received an OSL count");
                }
            }
        }
    }
    return true;
}

/**
 * @brief Verify marginal range diagnostics do not reject other observables.
 * @return true when one OSL overflow coexists with other filled marginals.
 */
bool verify_marginal_overflow_is_independent() {
    const hbt::HBTConfig config = make_config(
        {{hbt::AnalysisProduct{{
            hbt::PrimitiveChannelId::PiPlusPiPlus
        }}}},
        hbt::OriginMode::PrimordialRescatteringDecay,
        {{false, {}}, {false, {}}}
    );
    const hbt::HBTStartupState startup =
        hbt::build_hbt_startup_state(config);
    hbt::RawHistogramState state =
        hbt::make_zero_raw_histogram_state(config);
    hbt::RawHistogramAccumulator accumulator(config, startup, state);
    hbt::PairFrameObservables observables = in_range_observables();
    observables.r_out_lcms_fm = 5.0;
    observables.delta_t_lab_fm = -3.0;
    const hbt::PairFrameRouteContext context{
        0U,
        hbt::PrimitiveChannelId::PiPlusPiPlus,
        {false, false, true},
        nullptr
    };

    accumulator.consume(context, finite_kinematics(), observables);
    const hbt::RawHistogramSet& set =
        state.products[0].origins[0].global;

    if (set.osl.overflow_counts[0] != 1U ||
        set.delta_t.underflow_counts[0] != 1U) {
        return fail("range diagnostics were not recorded exactly");
    }
    if (count_at(set.osl.bins, 4U, 2U, 2U) != 1U ||
        count_at(set.delta_t.bins, 4U, 1U, 1U) != 1U) {
        return fail("one marginal overflow suppressed other observables");
    }
    return true;
}

/**
 * @brief Verify nested OriginMode::All routes fill only compatible origins.
 * @return true when PR and PRD fan-out does not fabricate primordial counts.
 */
bool verify_nested_origin_subset_fanout() {
    const hbt::HBTConfig config = make_config(
        {{hbt::AnalysisProduct{{
            hbt::PrimitiveChannelId::PiPlusPiPlus
        }}}},
        hbt::OriginMode::All,
        {{false, {}}, {false, {}}}
    );
    const hbt::HBTStartupState startup =
        hbt::build_hbt_startup_state(config);
    hbt::RawHistogramState state =
        hbt::make_zero_raw_histogram_state(config);
    hbt::RawHistogramAccumulator accumulator(config, startup, state);

    const hbt::PairFrameRouteContext pr_context{
        0U,
        hbt::PrimitiveChannelId::PiPlusPiPlus,
        {false, true, true},
        nullptr
    };
    accumulator.consume(
        pr_context,
        finite_kinematics(),
        in_range_observables()
    );

    const hbt::PairFrameRouteContext prd_context{
        0U,
        hbt::PrimitiveChannelId::PiPlusPiPlus,
        {false, false, true},
        nullptr
    };
    accumulator.consume(
        prd_context,
        finite_kinematics(),
        in_range_observables()
    );

    const auto first_osl_count = [](
        const hbt::OriginRawHistogramState& origin
    ) {
        return count_at(origin.global.osl.bins, 4U, 0U, 0U);
    };
    const hbt::ProductRawHistogramState& product = state.products[0];
    if (!product.origins[0].slices.empty() ||
        !product.origins[1].slices.empty() ||
        !product.origins[2].slices.empty()) {
        return fail("disabled slicing allocated slice histogram state");
    }
    if (first_osl_count(product.origins[0]) != 0U ||
        first_osl_count(product.origins[1]) != 1U ||
        first_osl_count(product.origins[2]) != 2U) {
        return fail("nested origin routes filled incompatible destinations");
    }
    return true;
}

/**
 * @brief Verify configured histogram ranges are exactly half-open.
 * @return true when minimum is included and maximum becomes overflow.
 */
bool verify_half_open_bin_boundaries() {
    const hbt::HBTConfig config = make_config(
        {{hbt::AnalysisProduct{{
            hbt::PrimitiveChannelId::PiPlusPiPlus
        }}}},
        hbt::OriginMode::Primordial,
        {{false, {}}, {false, {}}}
    );
    const hbt::HBTStartupState startup =
        hbt::build_hbt_startup_state(config);
    hbt::RawHistogramState state =
        hbt::make_zero_raw_histogram_state(config);
    hbt::RawHistogramAccumulator accumulator(config, startup, state);
    hbt::PairFrameObservables observables = in_range_observables();
    observables.r_out_lcms_fm = 0.0;
    observables.r_out_prf_fm = -4.0;
    observables.delta_t_lab_fm = -2.0;
    observables.delta_t_lcms_fm = 2.0;
    const hbt::PairFrameRouteContext context{
        0U,
        hbt::PrimitiveChannelId::PiPlusPiPlus,
        {true, false, false},
        nullptr
    };

    accumulator.consume(context, finite_kinematics(), observables);
    const hbt::RawHistogramSet& set =
        state.products[0].origins[0].global;
    if (count_at(set.osl.bins, 4U, 0U, 0U) != 1U ||
        set.osl.overflow_counts[1] != 1U ||
        count_at(set.delta_t.bins, 4U, 0U, 0U) != 1U ||
        set.delta_t.overflow_counts[1] != 1U) {
        return fail("histogram range is not exactly half-open");
    }
    return true;
}

/**
 * @brief Verify total histogram allocation arithmetic rejects overflow early.
 * @return true when an impossible family size fails before allocation.
 */
bool verify_allocation_size_overflow_is_explicit() {
    hbt::HBTConfig config = make_config(
        {{hbt::AnalysisProduct{{
            hbt::PrimitiveChannelId::PiPlusPiPlus
        }}}},
        hbt::OriginMode::Primordial,
        {{false, {}}, {false, {}}}
    );
    config.histogram_config.osl.nbins =
        std::numeric_limits<std::size_t>::max();

    try {
        static_cast<void>(hbt::make_zero_raw_histogram_state(config));
    } catch (const std::overflow_error&) {
        return true;
    } catch (...) {
        return fail("allocation overflow changed exception type");
    }
    return fail("allocation-size overflow was not rejected");
}

/**
 * @brief Verify worker-private raw states reduce into exact homologous bins.
 * @return true when bins and range counters are added exactly once.
 */
bool verify_worker_state_merge_is_exact() {
    const hbt::HBTConfig config = make_config(
        {{hbt::AnalysisProduct{{
            hbt::PrimitiveChannelId::PiPlusPiPlus
        }}}},
        hbt::OriginMode::Primordial,
        {{false, {}}, {true, {0.0, 1.0, 2.0}}}
    );
    hbt::RawHistogramState total =
        hbt::make_zero_raw_histogram_state(config);
    hbt::RawHistogramState first =
        hbt::make_zero_raw_histogram_state(config);
    hbt::RawHistogramState second =
        hbt::make_zero_raw_histogram_state(config);

    first.products[0].origins[0].global.osl.bins[0] = 3U;
    first.products[0].origins[0].global.radial.underflow_counts[1] = 2U;
    first.products[0].origins[0].slices[0].osl.bins[1] = 11U;
    first.products[0].origins[0].slices[1].delta_t.bins[2] = 17U;
    second.products[0].origins[0].global.osl.bins[0] = 5U;
    second.products[0].origins[0].global.delta_t.overflow_counts[2] = 7U;
    second.products[0].origins[0].slices[0].osl.bins[1] = 13U;
    second.products[0].origins[0].slices[1].delta_t.bins[2] = 19U;

    hbt::accumulate_raw_histogram_state(config, total, first);
    hbt::accumulate_raw_histogram_state(config, total, second);

    const hbt::OriginRawHistogramState& origin =
        total.products[0].origins[0];
    return origin.global.osl.bins[0] == 8U &&
           origin.global.radial.underflow_counts[1] == 2U &&
           origin.global.delta_t.overflow_counts[2] == 7U &&
           origin.slices[0].osl.bins[1] == 24U &&
           origin.slices[1].delta_t.bins[2] == 36U &&
           origin.slices[1].osl.bins[1] == 0U &&
           origin.slices[0].delta_t.bins[2] == 0U;
}

/**
 * @brief Verify a failed worker-state merge leaves the destination unchanged.
 * @return true when overflow is detected before any destination mutation.
 */
bool verify_worker_state_merge_overflow_is_atomic() {
    const hbt::HBTConfig config = make_config(
        {{hbt::AnalysisProduct{{
            hbt::PrimitiveChannelId::PiPlusPiPlus
        }}}},
        hbt::OriginMode::Primordial,
        {{false, {}}, {false, {}}}
    );
    hbt::RawHistogramState total =
        hbt::make_zero_raw_histogram_state(config);
    hbt::RawHistogramState worker =
        hbt::make_zero_raw_histogram_state(config);
    auto& total_bins = total.products[0].origins[0].global.osl.bins;
    auto& worker_bins = worker.products[0].origins[0].global.osl.bins;
    total_bins[0] = 11U;
    worker_bins[0] = 13U;
    total_bins[1] = std::numeric_limits<std::uint64_t>::max();
    worker_bins[1] = 1U;

    try {
        hbt::accumulate_raw_histogram_state(config, total, worker);
    } catch (const std::overflow_error&) {
        return total_bins[0] == 11U &&
               total_bins[1] ==
                   std::numeric_limits<std::uint64_t>::max();
    } catch (...) {
        return fail("worker merge overflow changed exception type");
    }
    return fail("worker merge overflow was not rejected");
}

/**
 * @brief Verify accumulator startup channels must match configuration.
 * @return true when inconsistent startup channel state is rejected.
 */
bool verify_startup_channel_mismatch_is_rejected() {
    const hbt::HBTConfig config = make_config(
        {{hbt::AnalysisProduct{{
            hbt::PrimitiveChannelId::PiPlusPiPlus
        }}}},
        hbt::OriginMode::Primordial,
        {{false, {}}, {false, {}}}
    );
    hbt::HBTStartupState startup =
        hbt::build_hbt_startup_state(config);
    startup.required_primitive_channels.clear();
    hbt::RawHistogramState state =
        hbt::make_zero_raw_histogram_state(config);

    try {
        hbt::RawHistogramAccumulator accumulator(config, startup, state);
        static_cast<void>(accumulator);
    } catch (const std::invalid_argument&) {
        return true;
    } catch (...) {
        return fail("startup mismatch changed exception type");
    }
    return fail("startup channel mismatch was not rejected");
}

/**
 * @brief Verify raw counter overflow fails before uint64_t wraparound.
 * @return true when incrementing a saturated first OSL bin throws.
 */
bool verify_counter_overflow_is_explicit() {
    const hbt::HBTConfig config = make_config(
        {{hbt::AnalysisProduct{{
            hbt::PrimitiveChannelId::PiPlusPiPlus
        }}}},
        hbt::OriginMode::Primordial,
        {{false, {}}, {false, {}}}
    );
    const hbt::HBTStartupState startup =
        hbt::build_hbt_startup_state(config);
    hbt::RawHistogramState state =
        hbt::make_zero_raw_histogram_state(config);
    state.products[0].origins[0].global.osl.bins[0] =
        std::numeric_limits<std::uint64_t>::max();
    hbt::RawHistogramAccumulator accumulator(config, startup, state);
    const hbt::PairFrameRouteContext context{
        0U,
        hbt::PrimitiveChannelId::PiPlusPiPlus,
        {true, false, false},
        nullptr
    };

    try {
        accumulator.consume(
            context,
            finite_kinematics(),
            in_range_observables()
        );
    } catch (const std::overflow_error&) {
        return true;
    } catch (...) {
        return fail("counter overflow changed exception type");
    }
    return fail("saturated raw counter wrapped silently");
}

}  // namespace

/**
 * @brief Run all Phase-6 raw histogram unit tests.
 * @return EXIT_SUCCESS when every test passes, otherwise EXIT_FAILURE.
 */
int main() {
    bool success = true;
    success = verify_eager_layout() && success;
    success = verify_shared_product_origin_and_slice_fanout() && success;
    success = verify_marginal_overflow_is_independent() && success;
    success = verify_nested_origin_subset_fanout() && success;
    success = verify_half_open_bin_boundaries() && success;
    success = verify_allocation_size_overflow_is_explicit() && success;
    success = verify_worker_state_merge_is_exact() && success;
    success = verify_worker_state_merge_overflow_is_atomic() && success;
    success = verify_startup_channel_mismatch_is_rejected() && success;
    success = verify_counter_overflow_is_explicit() && success;
    return success ? EXIT_SUCCESS : EXIT_FAILURE;
}
