/**
 * @file analysis_product_parser_test.cpp
 * @brief Unit tests for parsing one HBT analysis-product expression.
 *
 * This test verifies conversion of one textual analysis-product expression
 * into an AnalysisProduct. It also verifies the defined rejection cases for
 * invalid single-product expressions.
 */

#include "hbt/config/analysis_product_parser.h"

#include <cstdlib>
#include <iostream>
#include <stdexcept>
#include <string_view>
#include <vector>

namespace {

    /**
     * @brief Verify the primitive-channel sequence of one parsed product.
     *
     * @param test_name Human-readable test-case name.
     * @param expression Product expression to parse.
     * @param expected Expected primitive channels in exact order.
     *
     * @return true when parsing succeeds and the resulting channel sequence
     *         matches @p expected, otherwise false.
     */
    bool verify_product(
        const char* test_name,
        std::string_view expression,
        const std::vector<hbt::PrimitiveChannelId>& expected
    ) {
        try {
            const hbt::AnalysisProduct product =
                hbt::parse_analysis_product(expression);

            if (product.primitive_channels.size() != expected.size()) {
                std::cerr
                    << "analysis_product_parser_test: "
                    << test_name
                    << " produced "
                    << product.primitive_channels.size()
                    << " channels; expected "
                    << expected.size()
                    << ".\n";
                return false;
            }

            for (std::size_t index = 0; index < expected.size(); ++index) {
                if (product.primitive_channels[index] != expected[index]) {
                    std::cerr
                        << "analysis_product_parser_test: "
                        << test_name
                        << " produced the wrong channel at index "
                        << index
                        << ".\n";
                    return false;
                }
            }

            return true;
        } catch (const std::exception& exception) {
            std::cerr
                << "analysis_product_parser_test: "
                << test_name
                << " unexpectedly threw: "
                << exception.what()
                << ".\n";
            return false;
        }
    }

    /**
     * @brief Verify that one invalid expression throws std::invalid_argument.
     *
     * @param test_name Human-readable test-case name.
     * @param expression Invalid product expression.
     *
     * @return true when std::invalid_argument is thrown, otherwise false.
     */
    bool verify_invalid_expression(
        const char* test_name,
        std::string_view expression
    ) {
        try {
            static_cast<void>(
                hbt::parse_analysis_product(expression)
            );
        } catch (const std::invalid_argument&) {
            return true;
        } catch (const std::exception& exception) {
            std::cerr
                << "analysis_product_parser_test: "
                << test_name
                << " threw an unexpected exception type: "
                << exception.what()
                << ".\n";
            return false;
        }

        std::cerr
            << "analysis_product_parser_test: "
            << test_name
            << " did not throw std::invalid_argument.\n";
        return false;
    }

    /**
     * @brief Verify parsing of one primitive channel.
     *
     * @return true when the product contains exactly the expected channel.
     */
    bool verify_single_channel() {
        return verify_product(
            "single channel",
            "k_plus_p",
            {
                hbt::PrimitiveChannelId::KPlusProton
            }
        );
    }

    /**
     * @brief Verify parsing of multiple channels separated by '+'.
     *
     * @return true when both channels are preserved in expression order.
     */
    bool verify_multiple_channels() {
        return verify_product(
            "multiple channels",
            "k_plus_p+k_minus_p_bar",
            {
                hbt::PrimitiveChannelId::KPlusProton,
                hbt::PrimitiveChannelId::KMinusProtonBar
            }
        );
    }

    /**
     * @brief Verify surrounding ASCII spaces around primitive-channel tokens.
     *
     * @return true when spaces around tokens are ignored.
     */
    bool verify_spaces_around_tokens() {
        return verify_product(
            "spaces around tokens",
            "  k_plus_p  +  k_minus_p_bar  ",
            {
                hbt::PrimitiveChannelId::KPlusProton,
                hbt::PrimitiveChannelId::KMinusProtonBar
            }
        );
    }

    /**
     * @brief Verify preservation of the configured ASCII product expression.
     * @return true when only surrounding horizontal whitespace is removed.
     */
    bool verify_configured_expression_metadata() {
        const hbt::AnalysisProduct product = hbt::parse_analysis_product(
            "  k_plus_p  +  k_minus_p_bar  "
        );
        if (product.configured_expression !=
            "k_plus_p  +  k_minus_p_bar") {
            std::cerr
                << "analysis_product_parser_test: configured expression "
                << "metadata was not preserved.\n";
            return false;
        }
        return true;
    }

    /**
     * @brief Verify surrounding horizontal tabs around primitive-channel
     * tokens.
     *
     * @return true when horizontal tabs around tokens are ignored.
     */
    bool verify_tabs_around_tokens() {
        return verify_product(
            "tabs around tokens",
            "\tk_plus_p\t+\tk_minus_p_bar\t",
            {
                hbt::PrimitiveChannelId::KPlusProton,
                hbt::PrimitiveChannelId::KMinusProtonBar
            }
        );
    }

    /**
     * @brief Verify that primitive-channel order is preserved.
     *
     * @return true when the parsed order matches the expression order.
     */
    bool verify_channel_order() {
        return verify_product(
            "channel order",
            "pi_plus_p+k_minus_p+lambda_lambda",
            {
                hbt::PrimitiveChannelId::PiPlusProton,
                hbt::PrimitiveChannelId::KMinusProton,
                hbt::PrimitiveChannelId::LambdaLambda
            }
        );
    }

    /**
     * @brief Verify rejection of a duplicate primitive channel.
     *
     * @return true when repeated channels are rejected explicitly.
     */
    bool verify_duplicate_channels_rejected() {
        return verify_invalid_expression(
            "duplicate channels",
            "k_plus_p+k_plus_p"
        );
    }

    /**
     * @brief Verify rejection of an empty expression.
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
     * @brief Verify rejection of a whitespace-only expression.
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
     * @brief Verify rejection of a leading empty primitive-channel token.
     *
     * @return true when the expression is rejected.
     */
    bool verify_leading_separator() {
        return verify_invalid_expression(
            "leading separator",
            "+k_plus_p"
        );
    }

    /**
     * @brief Verify rejection of a trailing empty primitive-channel token.
     *
     * @return true when the expression is rejected.
     */
    bool verify_trailing_separator() {
        return verify_invalid_expression(
            "trailing separator",
            "k_plus_p+"
        );
    }

    /**
     * @brief Verify rejection of an empty token between two separators.
     *
     * @return true when the expression is rejected.
     */
    bool verify_double_separator() {
        return verify_invalid_expression(
            "double separator",
            "k_plus_p++p_p"
        );
    }

    /**
     * @brief Verify rejection of multiple-product syntax.
     *
     * @return true when an expression containing ',' is rejected.
     */
    bool verify_comma_rejected() {
        return verify_invalid_expression(
            "comma rejected",
            "k_plus_p,p_p"
        );
    }

    /**
     * @brief Verify rejection of an unknown primitive-channel token.
     *
     * @return true when the unknown token is rejected.
     */
    bool verify_unknown_token() {
        return verify_invalid_expression(
            "unknown token",
            "unknown_channel"
        );
    }

    /**
     * @brief Verify that non-horizontal whitespace is not trimmed.
     *
     * @return true when a token containing a newline is rejected.
     */
    bool verify_newline_not_trimmed() {
        return verify_invalid_expression(
            "newline not trimmed",
            "\nk_plus_p"
        );
    }

}  // namespace

/**
 * @brief Run the analysis-product parser unit tests.
 *
 * @return EXIT_SUCCESS when every test passes, otherwise EXIT_FAILURE.
 */
int main() {
    bool success = true;

    if (!verify_single_channel()) {
        success = false;
    }

    if (!verify_multiple_channels()) {
        success = false;
    }

    if (!verify_spaces_around_tokens()) {
        success = false;
    }

    if (!verify_configured_expression_metadata()) {
        success = false;
    }

    if (!verify_tabs_around_tokens()) {
        success = false;
    }

    if (!verify_channel_order()) {
        success = false;
    }

    if (!verify_duplicate_channels_rejected()) {
        success = false;
    }

    if (!verify_empty_expression()) {
        success = false;
    }

    if (!verify_whitespace_only_expression()) {
        success = false;
    }

    if (!verify_leading_separator()) {
        success = false;
    }

    if (!verify_trailing_separator()) {
        success = false;
    }

    if (!verify_double_separator()) {
        success = false;
    }

    if (!verify_comma_rejected()) {
        success = false;
    }

    if (!verify_unknown_token()) {
        success = false;
    }

    if (!verify_newline_not_trimmed()) {
        success = false;
    }

    return success ? EXIT_SUCCESS : EXIT_FAILURE;
}
