/**
 * @file product_fanout_plan.h
 * @brief Startup-resolved primitive-channel to final-product fan-out.
 */

#ifndef HBT_HISTOGRAMS_PRODUCT_FANOUT_PLAN_H
#define HBT_HISTOGRAMS_PRODUCT_FANOUT_PLAN_H

#include "hbt/channels/primitive_channel.h"
#include "hbt/selection/hbt_selection.h"

#include <cstddef>
#include <vector>

namespace hbt {

/**
 * @brief Compact startup-resolved product destinations for every channel.
 *
 * For required channel i, product_indices in the half-open range
 * [offsets[i], offsets[i + 1]) identify every final product fed by that
 * channel. The plan contains no scientific state and is immutable after
 * startup construction.
 */
struct ProductFanoutPlan {
    /// CSR-style offsets, with required-channel count plus one entries.
    std::vector<std::size_t> offsets;
    /// Final product indices grouped by required primitive channel.
    std::vector<std::size_t> product_indices;
};

/**
 * @brief Resolve channel-to-product fan-out once before event processing.
 * @param selection Validated final-product selection.
 * @param required_channels Unique required channels in startup order.
 * @return Compact fan-out plan indexed directly by channel_index.
 * @throws std::invalid_argument If a required channel is invalid or absent
 *         from every selected product.
 * @throws std::logic_error If one product contains a duplicate channel.
 *
 * The operation may scan the small startup collections. No lookup or
 * allocation performed here is repeated in the pair-processing hot path.
 */
[[nodiscard]] ProductFanoutPlan build_product_fanout_plan(
    const HBTSelection& selection,
    const std::vector<PrimitiveChannelId>& required_channels
);

}  // namespace hbt

#endif  // HBT_HISTOGRAMS_PRODUCT_FANOUT_PLAN_H
