/**
 * @file product_fanout_plan.cpp
 * @brief Startup-resolved primitive-channel to final-product fan-out.
 */

#include "hbt/histograms/product_fanout_plan.h"

#include "hbt/channels/channel_catalog.h"

#include <stdexcept>

namespace hbt {

ProductFanoutPlan build_product_fanout_plan(
    const HBTSelection& selection,
    const std::vector<PrimitiveChannelId>& required_channels
) {
    ProductFanoutPlan plan;
    plan.offsets.reserve(required_channels.size() + 1U);
    plan.offsets.push_back(0U);

    for (const PrimitiveChannelId channel : required_channels) {
        static_cast<void>(primitive_channel_definition(channel));
        bool channel_used = false;

        for (std::size_t product_index = 0U;
             product_index < selection.products.size();
             ++product_index) {
            const AnalysisProduct& product =
                selection.products[product_index];
            bool found = false;

            for (const PrimitiveChannelId product_channel :
                 product.primitive_channels) {
                if (product_channel != channel) {
                    continue;
                }
                if (found) {
                    throw std::logic_error(
                        "product fan-out: duplicate channel in product"
                    );
                }
                found = true;
            }

            if (found) {
                plan.product_indices.push_back(product_index);
                channel_used = true;
            }
        }

        if (!channel_used) {
            throw std::invalid_argument(
                "product fan-out: required channel has no final product"
            );
        }
        plan.offsets.push_back(plan.product_indices.size());
    }

    return plan;
}

}  // namespace hbt
