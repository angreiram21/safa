/**
 * @file primitive_channel_token_test.cpp
 * @brief Unit tests for primitive HBT channel configuration tokens.
 *
 * This test verifies the exact mapping between the canonical ASCII
 * configuration tokens and all primitive HBT channel identifiers. It also
 * verifies that non-canonical token forms are rejected.
 */

#include "hbt/config/primitive_channel_token.h"

#include <array>
#include <cstdlib>
#include <iostream>
#include <string_view>

namespace {

    /**
     * @brief Expected mapping between one configuration token and channel ID.
     */
    struct TokenExpectation {
        std::string_view token;             ///< Exact canonical ASCII token.
        hbt::PrimitiveChannelId channel;     ///< Expected primitive channel ID.
    };

    /**
     * @brief Verify all canonical primitive-channel tokens.
     *
     * @return true when every canonical token resolves to the expected channel
     *         identifier, otherwise false.
     */
    bool verify_canonical_tokens() {
        constexpr std::array<TokenExpectation, 21> expectations{{
            {"pi_plus_pi_plus", hbt::PrimitiveChannelId::PiPlusPiPlus},
            {"pi_minus_pi_minus", hbt::PrimitiveChannelId::PiMinusPiMinus},
            {"pi_zero_pi_zero", hbt::PrimitiveChannelId::PiZeroPiZero},

            {"k_plus_k_plus", hbt::PrimitiveChannelId::KPlusKPlus},
            {"k_minus_k_minus", hbt::PrimitiveChannelId::KMinusKMinus},
            {"k_zero_k_zero", hbt::PrimitiveChannelId::KZeroKZero},
            {"k_zero_bar_k_zero_bar",
             hbt::PrimitiveChannelId::KZeroBarKZeroBar},

            {"p_p", hbt::PrimitiveChannelId::ProtonProton},
            {"p_bar_p_bar", hbt::PrimitiveChannelId::ProtonBarProtonBar},

            {"lambda_lambda", hbt::PrimitiveChannelId::LambdaLambda},
            {"lambda_bar_lambda_bar",
             hbt::PrimitiveChannelId::LambdaBarLambdaBar},

            {"k_minus_p", hbt::PrimitiveChannelId::KMinusProton},
            {"k_plus_p", hbt::PrimitiveChannelId::KPlusProton},
            {"k_zero_bar_n", hbt::PrimitiveChannelId::KZeroBarNeutron},
            {"k_plus_p_bar", hbt::PrimitiveChannelId::KPlusProtonBar},
            {"k_minus_p_bar", hbt::PrimitiveChannelId::KMinusProtonBar},

            {"pi_plus_p", hbt::PrimitiveChannelId::PiPlusProton},
            {"pi_minus_p_bar", hbt::PrimitiveChannelId::PiMinusProtonBar},

            {"pi_minus_sigma_plus",
             hbt::PrimitiveChannelId::PiMinusSigmaPlus},
            {"pi_plus_sigma_bar_minus",
             hbt::PrimitiveChannelId::PiPlusSigmaBarMinus},
            {"pi_zero_sigma_zero",
             hbt::PrimitiveChannelId::PiZeroSigmaZero}
        }};

        for (const TokenExpectation& expectation : expectations) {
            const auto result =
                hbt::primitive_channel_from_token(expectation.token);

            if (!result.has_value()) {
                std::cerr
                    << "primitive_channel_token_test: canonical token '"
                    << expectation.token
                    << "' was not recognized.\n";
                return false;
            }

            if (*result != expectation.channel) {
                std::cerr
                    << "primitive_channel_token_test: canonical token '"
                    << expectation.token
                    << "' resolved to the wrong channel.\n";
                return false;
            }
        }

        return true;
    }

    /**
     * @brief Verify rejection of one non-canonical configuration token.
     *
     * @param token Token expected not to resolve to a primitive channel.
     *
     * @return true when @p token is rejected, otherwise false.
     */
    bool verify_rejected_token(std::string_view token) {
        const auto result = hbt::primitive_channel_from_token(token);

        if (result.has_value()) {
            std::cerr
                << "primitive_channel_token_test: non-canonical token '"
                << token
                << "' was unexpectedly recognized.\n";
            return false;
        }

        return true;
    }

    /**
     * @brief Verify that an empty token is rejected.
     *
     * @return true when the empty token is rejected, otherwise false.
     */
    bool verify_empty_token() {
        return verify_rejected_token("");
    }

    /**
     * @brief Verify that token matching is case-sensitive.
     *
     * @return true when a token containing uppercase characters is rejected,
     *         otherwise false.
     */
    bool verify_uppercase_token() {
        return verify_rejected_token("K_plus_p");
    }

    /**
     * @brief Verify that leading whitespace is not removed.
     *
     * @return true when a token with leading whitespace is rejected, otherwise
     *         false.
     */
    bool verify_leading_whitespace() {
        return verify_rejected_token(" k_plus_p");
    }

    /**
     * @brief Verify that trailing whitespace is not removed.
     *
     * @return true when a token with trailing whitespace is rejected, otherwise
     *         false.
     */
    bool verify_trailing_whitespace() {
        return verify_rejected_token("k_plus_p ");
    }

    /**
     * @brief Verify that a complete product expression is not accepted as
     * a token.
     *
     * @return true when an expression containing the product-composition
     *         separator is rejected, otherwise false.
     */
    bool verify_product_expression() {
        return verify_rejected_token("k_plus_p+k_minus_p_bar");
    }

    /**
     * @brief Verify that scientific-display notation is not a
     * configuration token.
     *
     * @return true when a Unicode scientific-display form is rejected,
     * otherwise
     *         false.
     */
    bool verify_scientific_display_token() {
        return verify_rejected_token("K⁺p");
    }

}  // namespace

/**
 * @brief Run the primitive-channel configuration-token unit tests.
 *
 * @return EXIT_SUCCESS when every test passes, otherwise EXIT_FAILURE.
 */
int main() {
    bool success = true;

    if (!verify_canonical_tokens()) {
        success = false;
    }

    if (!verify_empty_token()) {
        success = false;
    }

    if (!verify_uppercase_token()) {
        success = false;
    }

    if (!verify_leading_whitespace()) {
        success = false;
    }

    if (!verify_trailing_whitespace()) {
        success = false;
    }

    if (!verify_product_expression()) {
        success = false;
    }

    if (!verify_scientific_display_token()) {
        success = false;
    }

    return success ? EXIT_SUCCESS : EXIT_FAILURE;
}
