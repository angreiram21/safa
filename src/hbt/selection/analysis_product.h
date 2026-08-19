/**
 * @file analysis_product.h
 * @brief Definition of an HBT analysis product.
 *
 * An analysis product represents one requested final HBT result assembled
 * from one or more primitive HBT channels.
 *
 * This type describes the resolved composition and configured user-facing
 * identity of the product. Validation, pair construction, observable
 * calculation, histogramming, fitting, and output generation remain outside
 * its responsibility.
 */

#ifndef HBT_SELECTION_ANALYSIS_PRODUCT_H
#define HBT_SELECTION_ANALYSIS_PRODUCT_H

#include "hbt/channels/primitive_channel.h"

#include <string>
#include <vector>

namespace hbt {

    /**
     * @brief Composition of one requested final HBT analysis product.
     *
     * A product contains the primitive channels whose contributions are to be
     * accumulated into the same final result.
     *
     * A product containing one PrimitiveChannelId represents a primitive
     * result.
     * A product containing multiple PrimitiveChannelId values represents an
     * aggregate result.
     *
     * The primitive-channel list describes product composition only. It does
     * not imply that corresponding particle pairs should be recomputed
     * separately for each product.
     *
     * Product parsing preserves the primitive channels written by the user.
     * Validation of product structure belongs to the configuration layer.
     */
    struct AnalysisProduct {
        /**
         * @brief Primitive channels contributing to this final product.
         *
         * Each entry identifies one canonical primitive HBT channel.
         */
        std::vector<PrimitiveChannelId> primitive_channels;

        /**
         * @brief ASCII product expression preserved from configuration.
         *
         * Parsing removes surrounding horizontal whitespace from the complete
         * product expression once. Channel order and internal spacing remain
         * exactly as configured for output metadata and traceability.
         */
        std::string configured_expression{};
    };

}  // namespace hbt

#endif  // HBT_SELECTION_ANALYSIS_PRODUCT_H
