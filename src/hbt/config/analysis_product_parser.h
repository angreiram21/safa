/**
 * @file analysis_product_parser.h
 * @brief Parsing of one HBT analysis-product configuration expression.
 *
 * This file declares the operation used to convert one textual analysis-product
 * expression into an AnalysisProduct.
 *
 * Within one product expression, the '+' character separates primitive-channel
 * tokens. Parsing of multiple products separated by ',', YAML loading, complete
 * HBT selection construction and cross-product semantic validation are outside
 * this operation's responsibility. Duplicate primitive channels inside this one
 * product are rejected here because they would duplicate statistical weight.
 */

#ifndef HBT_CONFIG_ANALYSIS_PRODUCT_PARSER_H
#define HBT_CONFIG_ANALYSIS_PRODUCT_PARSER_H

#include "hbt/selection/analysis_product.h"

#include <string_view>

namespace hbt {

    /**
     * @brief Parse one HBT analysis-product configuration expression.
     *
     * Primitive-channel tokens are separated by '+'. Surrounding ASCII spaces
     * and horizontal tabs around each token are ignored before exact token
     * resolution.
     *
     * Channel order is preserved exactly as written in the expression. A
     * primitive channel may occur at most once in the product. Repetition is a
     * configuration error rather than an instruction to duplicate pair weight.
     *
     * The expression must describe exactly one product. Therefore ',' is not
     * accepted by this operation.
     *
     * @param expression Textual expression describing one HBT analysis product.
     *
     * @return AnalysisProduct containing the resolved primitive channels and
     *         the ASCII configured expression with surrounding horizontal
     *         whitespace removed.
     *
     * @throws std::invalid_argument if the expression is empty, contains an
     *         empty or duplicate primitive-channel token, contains ',', or
     *         contains an unrecognized primitive-channel token.
     */
    AnalysisProduct parse_analysis_product(
        std::string_view expression
    );

}  // namespace hbt

#endif  // HBT_CONFIG_ANALYSIS_PRODUCT_PARSER_H
