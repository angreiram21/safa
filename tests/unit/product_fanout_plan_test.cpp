/**
 * @file product_fanout_plan_test.cpp
 * @brief Unit tests for startup-resolved channel-to-product fan-out.
 */

#include "hbt/histograms/product_fanout_plan.h"

#include <cstdlib>
#include <iostream>
#include <vector>

namespace {

/**
 * @brief Report one failed product-fan-out condition.
 * @param message Description of the violated contract.
 * @return Always false.
 */
bool fail(const char* message) {
    std::cerr << "product_fanout_plan_test: " << message << ".\n";
    return false;
}

/**
 * @brief Verify shared channels fan out to all and only dependent products.
 * @return true when CSR offsets and product indices are exact.
 */
bool verify_shared_channel_fanout() {
    const hbt::HBTSelection selection{{
        hbt::AnalysisProduct{{
            hbt::PrimitiveChannelId::PiPlusPiPlus,
            hbt::PrimitiveChannelId::PiMinusPiMinus
        }},
        hbt::AnalysisProduct{{
            hbt::PrimitiveChannelId::PiPlusPiPlus,
            hbt::PrimitiveChannelId::KPlusKPlus
        }},
        hbt::AnalysisProduct{{
            hbt::PrimitiveChannelId::ProtonProton
        }}
    }};
    const std::vector<hbt::PrimitiveChannelId> channels{
        hbt::PrimitiveChannelId::PiPlusPiPlus,
        hbt::PrimitiveChannelId::PiMinusPiMinus,
        hbt::PrimitiveChannelId::KPlusKPlus,
        hbt::PrimitiveChannelId::ProtonProton
    };

    const hbt::ProductFanoutPlan plan =
        hbt::build_product_fanout_plan(selection, channels);
    const std::vector<std::size_t> expected_offsets{0U, 2U, 3U, 4U, 5U};
    const std::vector<std::size_t> expected_products{0U, 1U, 0U, 1U, 2U};

    if (plan.offsets != expected_offsets) {
        return fail("channel offsets are incorrect");
    }
    if (plan.product_indices != expected_products) {
        return fail("product destinations are incorrect");
    }
    return true;
}

/**
 * @brief Verify one product creates no implicit per-channel product states.
 * @return true when both primitive channels point to product index zero.
 */
bool verify_aggregate_product_has_one_destination() {
    const hbt::HBTSelection selection{{
        hbt::AnalysisProduct{{
            hbt::PrimitiveChannelId::PiPlusPiPlus,
            hbt::PrimitiveChannelId::PiMinusPiMinus
        }}
    }};
    const std::vector<hbt::PrimitiveChannelId> channels{
        hbt::PrimitiveChannelId::PiPlusPiPlus,
        hbt::PrimitiveChannelId::PiMinusPiMinus
    };

    const hbt::ProductFanoutPlan plan =
        hbt::build_product_fanout_plan(selection, channels);
    return plan.offsets == std::vector<std::size_t>{0U, 1U, 2U} &&
           plan.product_indices == std::vector<std::size_t>{0U, 0U};
}

}  // namespace

/**
 * @brief Run all product-fan-out unit tests.
 * @return EXIT_SUCCESS when every test passes, otherwise EXIT_FAILURE.
 */
int main() {
    bool success = true;
    success = verify_shared_channel_fanout() && success;
    success = verify_aggregate_product_has_one_destination() && success;
    return success ? EXIT_SUCCESS : EXIT_FAILURE;
}
