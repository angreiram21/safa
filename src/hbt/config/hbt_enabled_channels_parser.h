/**
 * @file hbt_enabled_channels_parser.h
 * @brief Parsing of the complete HBT enabled-channels configuration expression.
 *
 * This file declares the operation used to convert the textual value of
 * hbt_enabled_channels into an HBTSelection.
 *
 * Commas separate final analysis products. Parsing of each individual product
 * expression is delegated to the single-product parser, where '+' separates
 * primitive channels belonging to the same final product.
 *
 * YAML loading and semantic validation other than final-product uniqueness are
 * outside this operation's responsibility.
 */

#ifndef HBT_CONFIG_HBT_ENABLED_CHANNELS_PARSER_H
#define HBT_CONFIG_HBT_ENABLED_CHANNELS_PARSER_H

#include "hbt/selection/hbt_selection.h"

#include <string_view>

namespace hbt {

    /**
     * @brief Parse the complete hbt_enabled_channels expression.
     *
     * Final analysis products are separated by ','. Each resulting product
     * expression is parsed according to the single-product configuration
     * contract.
     *
     * Product order is preserved exactly as written in the expression. Each
     * parsed product also retains its ASCII configured expression as metadata.
     * Final products are compared canonically without changing their visible
     * order; repeated products and order-equivalent products are rejected.
     *
     * Surrounding ASCII spaces and horizontal tabs are accepted according to
     * the single-product parser rules.
     *
     * @param expression Textual value of hbt_enabled_channels.
     *
     * @return HBTSelection containing the parsed final products in expression
     *         order.
     *
     * @throws std::invalid_argument if the expression is empty, contains an
     *         empty product expression, contains a duplicate final product, or
     *         contains any product expression that is invalid according to the
     *         single-product parser contract.
     */
    HBTSelection parse_hbt_enabled_channels(
        std::string_view expression
    );

}  // namespace hbt

#endif  // HBT_CONFIG_HBT_ENABLED_CHANNELS_PARSER_H
