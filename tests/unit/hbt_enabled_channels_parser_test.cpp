/**
 * @file hbt_enabled_channels_parser_test.cpp
 * @brief Unit tests for parsing the complete HBT enabled-channels expression.
 *
 * This test verifies conversion of the textual hbt_enabled_channels value into
 * an HBTSelection. It also verifies the defined rejection cases for invalid
 * complete expressions.
 */

#include "hbt/config/hbt_enabled_channels_parser.h"

#include <cstdlib>
#include <exception>
#include <iostream>
#include <stdexcept>
#include <string_view>
#include <vector>

namespace {

    /**
     * @brief Verify the products and primitive-channel sequences of one
     * selection.
     *
     * @param test_name Human-readable test-case name.
     * @param expression Complete enabled-channels expression to parse.
     * @param expected Expected primitive channels for every product, in exact
     *        product and channel order.
     *
     * @return true when parsing succeeds and the complete selection matches
     *         @p expected, otherwise false.
     */
    bool verify_selection(
        const char* test_name,
        std::string_view expression,
        const std::vector<std::vector<hbt::PrimitiveChannelId>>& expected
    ) {
        try {
            const hbt::HBTSelection selection =
                hbt::parse_hbt_enabled_channels(expression);

            if (selection.products.size() != expected.size()) {
                std::cerr
                    << "hbt_enabled_channels_parser_test: "
                    << test_name
                    << " produced "
                    << selection.products.size()
                    << " products; expected "
                    << expected.size()
                    << ".\n";
                return false;
            }

            for (
                std::size_t product_index = 0;
                product_index < expected.size();
                ++product_index
            ) {
                const std::vector<hbt::PrimitiveChannelId>& actual_channels =
                    selection.products[product_index].primitive_channels;

                const std::vector<hbt::PrimitiveChannelId>& expected_channels =
                    expected[product_index];

                if (actual_channels.size() != expected_channels.size()) {
                    std::cerr
                        << "hbt_enabled_channels_parser_test: "
                        << test_name
                        << " produced "
                        << actual_channels.size()
                        << " channels in product "
                        << product_index
                        << "; expected "
                        << expected_channels.size()
                        << ".\n";
                    return false;
                }

                for (
                    std::size_t channel_index = 0;
                    channel_index < expected_channels.size();
                    ++channel_index
                ) {
                    if (
                        actual_channels[channel_index] !=
                        expected_channels[channel_index]
                    ) {
                        std::cerr
                            << "hbt_enabled_channels_parser_test: "
                            << test_name
                            << " produced the wrong channel at product "
                            << product_index
                            << ", channel "
                            << channel_index
                            << ".\n";
                        return false;
                    }
                }
            }

            return true;
        } catch (const std::exception& exception) {
            std::cerr
                << "hbt_enabled_channels_parser_test: "
                << test_name
                << " unexpectedly threw: "
                << exception.what()
                << ".\n";
            return false;
        }
    }

    /**
     * @brief Verify that one invalid complete expression throws
     *        std::invalid_argument.
     *
     * @param test_name Human-readable test-case name.
     * @param expression Invalid enabled-channels expression.
     *
     * @return true when std::invalid_argument is thrown, otherwise false.
     */
    bool verify_invalid_expression(
        const char* test_name,
        std::string_view expression
    ) {
        try {
            static_cast<void>(
                hbt::parse_hbt_enabled_channels(expression)
            );
        } catch (const std::invalid_argument&) {
            return true;
        } catch (const std::exception& exception) {
            std::cerr
                << "hbt_enabled_channels_parser_test: "
                << test_name
                << " threw an unexpected exception type: "
                << exception.what()
                << ".\n";
            return false;
        }

        std::cerr
            << "hbt_enabled_channels_parser_test: "
            << test_name
            << " did not throw std::invalid_argument.\n";
        return false;
    }

    /**
     * @brief Verify parsing of one final product.
     *
     * @return true when exactly one product with the expected channel is
     *         produced.
     */
    bool verify_single_product() {
        return verify_selection(
            "single product",
            "pi_plus_pi_plus",
            {
                {
                    hbt::PrimitiveChannelId::PiPlusPiPlus
                }
            }
        );
    }

    /**
     * @brief Verify parsing of multiple final products separated by commas.
     *
     * @return true when all products are preserved in expression order.
     */
    bool verify_multiple_products() {
        return verify_selection(
            "multiple products",
            "pi_plus_pi_plus,k_plus_p,p_p",
            {
                {
                    hbt::PrimitiveChannelId::PiPlusPiPlus
                },
                {
                    hbt::PrimitiveChannelId::KPlusProton
                },
                {
                    hbt::PrimitiveChannelId::ProtonProton
                }
            }
        );
    }

    /**
     * @brief Verify a product containing multiple primitive channels.
     *
     * @return true when '+' composition remains inside one final product.
     */
    bool verify_composed_product() {
        return verify_selection(
            "composed product",
            "pi_plus_pi_plus,k_plus_p+k_minus_p_bar",
            {
                {
                    hbt::PrimitiveChannelId::PiPlusPiPlus
                },
                {
                    hbt::PrimitiveChannelId::KPlusProton,
                    hbt::PrimitiveChannelId::KMinusProtonBar
                }
            }
        );
    }

    /**
     * @brief Verify ASCII spaces around product separators and channel tokens.
     *
     * @return true when surrounding spaces are accepted through the
     *         single-product parser contract.
     */
    bool verify_spaces_around_products() {
        return verify_selection(
            "spaces around products",
            "  pi_plus_pi_plus  ,  k_plus_p + k_minus_p_bar  ",
            {
                {
                    hbt::PrimitiveChannelId::PiPlusPiPlus
                },
                {
                    hbt::PrimitiveChannelId::KPlusProton,
                    hbt::PrimitiveChannelId::KMinusProtonBar
                }
            }
        );
    }

    /**
     * @brief Verify product-level configured expressions survive list parsing.
     * @return true when each product retains its trimmed ASCII source text.
     */
    bool verify_product_expression_metadata() {
        const hbt::HBTSelection selection =
            hbt::parse_hbt_enabled_channels(
                "  pi_plus_pi_plus  ,  k_plus_p + k_minus_p_bar  "
            );
        if (selection.products.size() != 2U ||
            selection.products[0U].configured_expression !=
                "pi_plus_pi_plus" ||
            selection.products[1U].configured_expression !=
                "k_plus_p + k_minus_p_bar") {
            std::cerr
                << "hbt_enabled_channels_parser_test: product expression "
                << "metadata was not preserved.\n";
            return false;
        }
        return true;
    }

    /**
     * @brief Verify horizontal tabs around product separators and channel
     * tokens.
     *
     * @return true when surrounding tabs are accepted through the
     *         single-product parser contract.
     */
    bool verify_tabs_around_products() {
        return verify_selection(
            "tabs around products",
            "\tpi_plus_pi_plus\t,\tk_plus_p\t+\tk_minus_p_bar\t",
            {
                {
                    hbt::PrimitiveChannelId::PiPlusPiPlus
                },
                {
                    hbt::PrimitiveChannelId::KPlusProton,
                    hbt::PrimitiveChannelId::KMinusProtonBar
                }
            }
        );
    }

    /**
     * @brief Verify preservation of final-product order.
     *
     * @return true when product order exactly matches the input expression.
     */
    bool verify_product_order() {
        return verify_selection(
            "product order",
            "lambda_lambda,pi_plus_p,k_minus_p",
            {
                {
                    hbt::PrimitiveChannelId::LambdaLambda
                },
                {
                    hbt::PrimitiveChannelId::PiPlusProton
                },
                {
                    hbt::PrimitiveChannelId::KMinusProton
                }
            }
        );
    }

    /**
     * @brief Verify rejection of a repeated primitive channel in one product.
     *
     * @return true when A+A is rejected by complete-expression parsing.
     */
    bool verify_duplicate_channel_in_product() {
        return verify_invalid_expression(
            "duplicate channel in product",
            "k_plus_p+k_plus_p"
        );
    }

    /**
     * @brief Verify rejection of an exactly repeated final product.
     *
     * @return true when the repeated product is rejected.
     */
    bool verify_duplicate_products() {
        return verify_invalid_expression(
            "duplicate products",
            "k_plus_p,k_plus_p"
        );
    }

    /**
     * @brief Verify rejection of an order-equivalent final product.
     *
     * @return true when A+B and B+A are treated as the same final product.
     */
    bool verify_order_equivalent_duplicate_products() {
        return verify_invalid_expression(
            "order-equivalent duplicate products",
            "pi_plus_p+k_minus_p_bar,k_minus_p_bar+pi_plus_p"
        );
    }

    /**
     * @brief Verify rejection of an empty complete expression.
     *
     * @return true when the expression is rejected.
     */
    bool verify_empty_expression() {
        return verify_invalid_expression(
            "empty expression",
            ""
        );
    }

    /**
     * @brief Verify rejection of a whitespace-only complete expression.
     *
     * @return true when the expression is rejected.
     */
    bool verify_whitespace_only_expression() {
        return verify_invalid_expression(
            "whitespace-only expression",
            " \t "
        );
    }

    /**
     * @brief Verify rejection of a comma-only expression.
     *
     * @return true when the expression is rejected.
     */
    bool verify_comma_only_expression() {
        return verify_invalid_expression(
            "comma-only expression",
            ","
        );
    }

    /**
     * @brief Verify rejection of a leading empty product.
     *
     * @return true when the expression is rejected.
     */
    bool verify_leading_comma() {
        return verify_invalid_expression(
            "leading comma",
            ",k_plus_p"
        );
    }

    /**
     * @brief Verify rejection of a trailing empty product.
     *
     * @return true when the expression is rejected.
     */
    bool verify_trailing_comma() {
        return verify_invalid_expression(
            "trailing comma",
            "k_plus_p,"
        );
    }

    /**
     * @brief Verify rejection of an empty product between two commas.
     *
     * @return true when the expression is rejected.
     */
    bool verify_double_comma() {
        return verify_invalid_expression(
            "double comma",
            "k_plus_p,,p_p"
        );
    }

    /**
     * @brief Verify rejection of a whitespace-only product between commas.
     *
     * @return true when the expression is rejected.
     */
    bool verify_whitespace_only_product() {
        return verify_invalid_expression(
            "whitespace-only product",
            "k_plus_p, \t ,p_p"
        );
    }

    /**
     * @brief Verify propagation of an invalid primitive-channel token.
     *
     * @return true when the complete expression is rejected.
     */
    bool verify_unknown_token() {
        return verify_invalid_expression(
            "unknown token",
            "pi_plus_pi_plus,unknown_channel"
        );
    }

    /**
     * @brief Verify propagation of an invalid single-product expression.
     *
     * @return true when an empty primitive-channel token is rejected.
     */
    bool verify_invalid_composed_product() {
        return verify_invalid_expression(
            "invalid composed product",
            "pi_plus_pi_plus,k_plus_p++p_p"
        );
    }

    /**
     * @brief Verify that non-horizontal whitespace is not accepted implicitly.
     *
     * @return true when a product containing a newline is rejected.
     */
    bool verify_newline_not_trimmed() {
        return verify_invalid_expression(
            "newline not trimmed",
            "pi_plus_pi_plus,\nk_plus_p"
        );
    }

}  // namespace

/**
 * @brief Run the HBT enabled-channels parser unit tests.
 *
 * @return EXIT_SUCCESS when every test passes, otherwise EXIT_FAILURE.
 */
int main() {
    bool success = true;

    if (!verify_single_product()) {
        success = false;
    }

    if (!verify_multiple_products()) {
        success = false;
    }

    if (!verify_composed_product()) {
        success = false;
    }

    if (!verify_spaces_around_products()) {
        success = false;
    }

    if (!verify_product_expression_metadata()) {
        success = false;
    }

    if (!verify_tabs_around_products()) {
        success = false;
    }

    if (!verify_product_order()) {
        success = false;
    }

    if (!verify_duplicate_channel_in_product()) {
        success = false;
    }

    if (!verify_duplicate_products()) {
        success = false;
    }

    if (!verify_order_equivalent_duplicate_products()) {
        success = false;
    }

    if (!verify_empty_expression()) {
        success = false;
    }

    if (!verify_whitespace_only_expression()) {
        success = false;
    }

    if (!verify_comma_only_expression()) {
        success = false;
    }

    if (!verify_leading_comma()) {
        success = false;
    }

    if (!verify_trailing_comma()) {
        success = false;
    }

    if (!verify_double_comma()) {
        success = false;
    }

    if (!verify_whitespace_only_product()) {
        success = false;
    }

    if (!verify_unknown_token()) {
        success = false;
    }

    if (!verify_invalid_composed_product()) {
        success = false;
    }

    if (!verify_newline_not_trimmed()) {
        success = false;
    }

    return success ? EXIT_SUCCESS : EXIT_FAILURE;
}
