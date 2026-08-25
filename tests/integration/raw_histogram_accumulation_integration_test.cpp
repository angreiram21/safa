/**
 * @file raw_histogram_accumulation_integration_test.cpp
 * @brief Integration coverage for pair processing into Phase-6 raw histograms.
 */

#include "hbt/histograms/raw_histograms.h"
#include "hbt/pair/pair_processor.h"
#include "hbt/startup/hbt_startup_builder.h"

#include <cstdint>
#include <cstdlib>
#include <iostream>
#include <numeric>
#include <utility>

namespace {

/**
 * @brief Report one failed raw-histogram integration condition.
 * @param message Description of the violated contract.
 * @return Always false.
 */
bool fail(const char* message) {
    std::cerr << "raw_histogram_accumulation_integration_test: "
              << message << ".\n";
    return false;
}

/**
 * @brief Return two final products sharing the proton-proton primitive channel.
 * @return Product zero A+B and product one A+C with shared A.
 */
hbt::HBTSelection overlapping_selection() {
    return {{
        hbt::AnalysisProduct{{
            hbt::PrimitiveChannelId::ProtonProton,
            hbt::PrimitiveChannelId::PiPlusPiPlus
        }},
        hbt::AnalysisProduct{{
            hbt::PrimitiveChannelId::ProtonProton,
            hbt::PrimitiveChannelId::KPlusKPlus
        }}
    }};
}

/**
 * @brief Build the complete integration configuration.
 * @return Valid HBT configuration with all origins and one kT slice.
 */
hbt::HBTConfig make_config() {
    return {
        overlapping_selection(),
        {
            hbt::LongitudinalVariable::Pseudorapidity,
            {0.8, 0.14, 4.0},
            {0.8, 0.4, 1.4},
            {0.8, 0.5, 4.05},
            {0.8, 1.0, 10000.0},
            {0.8, 0.3, 10000.0}
        },
        {{true, {0.0, 0.3}}, {false, {}}},
        {
            {20U, 0.0, 10.0, 2.0},
            {20U, 0.0, 10.0, 2.0},
            {40U, -10.0, 10.0, 2.0}
        },
        hbt::OriginMode::All,
        hbt::FitEstimatorMode::All
    };
}

/**
 * @brief Build one primordial proton with finite emission and momentum data.
 * @param x Emission x coordinate in fm.
 * @return Complete synthetic proton accepted into EventBuffers.
 */
hbt::Particle make_proton(double x) {
    return {
        hbt::SpeciesId::Proton,
        {0.0, x, 0.0, 0.0},
        {2.0, 0.2, 0.0, 0.0},
        0.938,
        {true, true, true},
        2212,
        1
    };
}

/**
 * @brief Build one independent synthetic subevent with exactly one A pair.
 * @return Buffers containing two protons and no B/C-channel particles.
 */
hbt::EventBuffers make_subevent() {
    hbt::EventBuffers buffers;
    buffers.add(make_proton(1.0));
    buffers.add(make_proton(-1.0));
    return buffers;
}

/**
 * @brief Sum all counters in one flattened histogram family.
 * @param bins Raw flattened counter block.
 * @return Exact sum of all counters in @p bins.
 */
std::uint64_t sum_bins(const std::vector<std::uint64_t>& bins) {
    return std::accumulate(
        bins.begin(),
        bins.end(),
        std::uint64_t{0U}
    );
}

/**
 * @brief Verify one destination contains two complete pair contributions.
 * @param set Raw histogram destination accumulated from two subevents.
 * @return true when all nine marginals received exactly two in-range entries.
 */
bool has_two_complete_pair_contributions(const hbt::RawHistogramSet& set) {
    const std::uint64_t osl_total = sum_bins(set.osl.bins);
    const std::uint64_t radial_total = sum_bins(set.radial.bins);
    const std::uint64_t delta_t_total = sum_bins(set.delta_t.bins);

    if (osl_total != 8U || radial_total != 4U || delta_t_total != 6U) {
        return false;
    }
    for (const std::uint64_t count : set.osl.underflow_counts) {
        if (count != 0U) {
            return false;
        }
    }
    for (const std::uint64_t count : set.osl.overflow_counts) {
        if (count != 0U) {
            return false;
        }
    }
    for (const std::uint64_t count : set.radial.underflow_counts) {
        if (count != 0U) {
            return false;
        }
    }
    for (const std::uint64_t count : set.radial.overflow_counts) {
        if (count != 0U) {
            return false;
        }
    }
    for (const std::uint64_t count : set.delta_t.underflow_counts) {
        if (count != 0U) {
            return false;
        }
    }
    for (const std::uint64_t count : set.delta_t.overflow_counts) {
        if (count != 0U) {
            return false;
        }
    }
    return true;
}

/**
 * @brief Verify real pair processing fans out once into persistent raw state.
 * @return true when two independent subevents feed both products and all
 *         origin/global/slice destinations without intermediate products.
 */
bool verify_pair_processor_to_raw_histograms() {
    const hbt::HBTConfig config = make_config();
    const hbt::HBTStartupState startup =
        hbt::build_hbt_startup_state(config);
    hbt::RawHistogramState state =
        hbt::make_zero_raw_histogram_state(config);
    hbt::RawHistogramAccumulator accumulator(config, startup, state);

    for (int subevent_id = 0; subevent_id < 2; ++subevent_id) {
        const hbt::PairSubeventProcessingResult result =
            hbt::process_subevent_pairs(
                1U,
                subevent_id,
                make_subevent(),
                startup.required_primitive_channels,
                config.origin_mode,
                config.pair_slicing,
                accumulator
            );
        if (result.summary.pair_counts.channels[0].pair_count != 1U ||
            result.summary.valid_pair_counts.channels[0].pair_count != 1U ||
            result.summary.numerical_rejection_counts
                    .channels[0].pair_count != 0U) {
            return fail(
                "one synthetic subevent did not produce one valid A pair"
            );
        }
    }

    if (state.products.size() != 2U) {
        return fail("shared channel did not preserve two final products");
    }
    for (const hbt::ProductRawHistogramState& product : state.products) {
        if (product.origins.size() != 3U) {
            return fail("OriginMode::All state does not contain three origins");
        }
        for (const hbt::OriginRawHistogramState& origin : product.origins) {
            if (origin.slices.size() != 1U) {
                return fail("configured kT slicing did not allocate one slice");
            }
            if (!has_two_complete_pair_contributions(origin.global) ||
                !has_two_complete_pair_contributions(origin.slices[0])) {
                return fail(
                    "global and slice states did not reuse both pair results"
                );
            }
        }
    }
    return true;
}

}  // namespace

/**
 * @brief Run the raw-histogram pair-processing integration test.
 * @return EXIT_SUCCESS when the complete integration contract holds.
 */
int main() {
    return verify_pair_processor_to_raw_histograms()
               ? EXIT_SUCCESS
               : EXIT_FAILURE;
}
