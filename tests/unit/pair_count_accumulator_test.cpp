/**
 * @file pair_count_accumulator_test.cpp
 * @brief Unit tests for primitive-channel pair-count accumulation.
 */

#include "hbt/pair/pair_count_accumulator.h"

#include <cstdint>
#include <cstdlib>
#include <iostream>
#include <limits>
#include <stdexcept>
#include <vector>

namespace {

/**
 * @brief Report one failed pair-count-accumulator condition.
 * @param message Description of the violated contract.
 * @return Always false.
 */
bool fail(const char* message) {
    std::cerr << "pair_count_accumulator_test: " << message << ".\n";
    return false;
}

/**
 * @brief Verify ordered zero-summary construction.
 * @return true when channel order and zero counts are preserved exactly.
 */
bool verify_zero_summary() {
    const std::vector<hbt::PrimitiveChannelId> channels{
        hbt::PrimitiveChannelId::PiPlusPiPlus,
        hbt::PrimitiveChannelId::KMinusProton
    };

    const hbt::PairCountSummary summary =
        hbt::make_zero_pair_count_summary(channels);

    return summary.channels.size() == 2U &&
           summary.channels[0].channel == channels[0] &&
           summary.channels[0].pair_count == 0U &&
           summary.channels[1].channel == channels[1] &&
           summary.channels[1].pair_count == 0U;
}

/**
 * @brief Verify invalid or duplicate initialization channels are rejected.
 * @return true when both invalid inputs throw std::invalid_argument.
 */
bool verify_zero_summary_rejects_invalid_structure() {
    try {
        static_cast<void>(hbt::make_zero_pair_count_summary({
            hbt::PrimitiveChannelId::PiPlusPiPlus,
            hbt::PrimitiveChannelId::PiPlusPiPlus
        }));
        return fail("duplicate initialization channel was accepted");
    } catch (const std::invalid_argument&) {
    } catch (...) {
        return fail("duplicate initialization used wrong exception type");
    }

    try {
        static_cast<void>(hbt::make_zero_pair_count_summary({
            static_cast<hbt::PrimitiveChannelId>(999)
        }));
        return fail("invalid initialization channel was accepted");
    } catch (const std::invalid_argument&) {
        return true;
    } catch (...) {
        return fail("invalid initialization used wrong exception type");
    }
}

/**
 * @brief Verify matching summaries accumulate counter by counter.
 * @return true when all resulting counts equal their exact sums.
 */
bool verify_accumulation() {
    hbt::PairCountSummary total{{
        {hbt::PrimitiveChannelId::PiPlusPiPlus, 3U},
        {hbt::PrimitiveChannelId::KMinusProton, 5U}
    }};
    const hbt::PairCountSummary local{{
        {hbt::PrimitiveChannelId::PiPlusPiPlus, 7U},
        {hbt::PrimitiveChannelId::KMinusProton, 11U}
    }};

    hbt::accumulate_pair_counts(total, local);

    return total.channels[0].pair_count == 10U &&
           total.channels[1].pair_count == 16U;
}

/**
 * @brief Verify structural mismatches fail without modifying the total.
 * @return true when size and order mismatches both preserve the original total.
 */
bool verify_structure_mismatch_preserves_total() {
    const hbt::PairCountSummary original{{
        {hbt::PrimitiveChannelId::PiPlusPiPlus, 3U},
        {hbt::PrimitiveChannelId::KMinusProton, 5U}
    }};

    hbt::PairCountSummary total = original;
    try {
        hbt::accumulate_pair_counts(
            total,
            {{{hbt::PrimitiveChannelId::PiPlusPiPlus, 1U}}}
        );
        return fail("channel-count mismatch was accepted");
    } catch (const std::invalid_argument&) {
    } catch (...) {
        return fail("channel-count mismatch used wrong exception type");
    }

    if (total.channels[0].pair_count != 3U ||
        total.channels[1].pair_count != 5U) {
        return fail("size mismatch modified total counts");
    }

    try {
        hbt::accumulate_pair_counts(
            total,
            {{
                {hbt::PrimitiveChannelId::KMinusProton, 1U},
                {hbt::PrimitiveChannelId::PiPlusPiPlus, 1U}
            }}
        );
        return fail("channel-order mismatch was accepted");
    } catch (const std::invalid_argument&) {
    } catch (...) {
        return fail("channel-order mismatch used wrong exception type");
    }

    return total.channels[0].pair_count == 3U &&
           total.channels[1].pair_count == 5U;
}

/**
 * @brief Verify duplicate summary channels are rejected during accumulation.
 * @return true when duplicate structure throws before any total modification.
 */
bool verify_duplicate_summary_rejected() {
    hbt::PairCountSummary total{{
        {hbt::PrimitiveChannelId::PiPlusPiPlus, 3U},
        {hbt::PrimitiveChannelId::PiPlusPiPlus, 5U}
    }};
    const hbt::PairCountSummary local = total;

    try {
        hbt::accumulate_pair_counts(total, local);
    } catch (const std::invalid_argument&) {
        return total.channels[0].pair_count == 3U &&
               total.channels[1].pair_count == 5U;
    } catch (...) {
        return fail("duplicate summary used wrong exception type");
    }

    return fail("duplicate summary channels were accepted");
}

/**
 * @brief Verify uint64_t overflow is rejected before any mutation.
 * @return true when overflow throws and all total counters remain unchanged.
 */
bool verify_overflow_preserves_total() {
    const std::uint64_t maximum =
        std::numeric_limits<std::uint64_t>::max();
    hbt::PairCountSummary total{{
        {hbt::PrimitiveChannelId::PiPlusPiPlus, 10U},
        {hbt::PrimitiveChannelId::KMinusProton, maximum}
    }};
    const hbt::PairCountSummary local{{
        {hbt::PrimitiveChannelId::PiPlusPiPlus, 20U},
        {hbt::PrimitiveChannelId::KMinusProton, 1U}
    }};

    try {
        hbt::accumulate_pair_counts(total, local);
    } catch (const std::overflow_error&) {
        return total.channels[0].pair_count == 10U &&
               total.channels[1].pair_count == maximum;
    } catch (...) {
        return fail("overflow used wrong exception type");
    }

    return fail("pair-count overflow was accepted");
}

}  // namespace

/**
 * @brief Run all pair-count-accumulator unit tests.
 * @return EXIT_SUCCESS when every test passes, otherwise EXIT_FAILURE.
 */
int main() {
    bool success = true;
    success = verify_zero_summary() && success;
    success = verify_zero_summary_rejects_invalid_structure() && success;
    success = verify_accumulation() && success;
    success = verify_structure_mismatch_preserves_total() && success;
    success = verify_duplicate_summary_rejected() && success;
    success = verify_overflow_preserves_total() && success;
    return success ? EXIT_SUCCESS : EXIT_FAILURE;
}
