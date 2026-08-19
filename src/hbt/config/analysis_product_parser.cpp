/**
 * @file analysis_product_parser.cpp
 * @brief Parsing of one HBT analysis-product configuration expression.
 */

#include "hbt/config/analysis_product_parser.h"

#include "hbt/config/primitive_channel_token.h"

#include <algorithm>
#include <stdexcept>
#include <string>
#include <string_view>

namespace hbt {

namespace {

    /**
     * @brief Remove surrounding ASCII horizontal whitespace from a token.
     *
     * Only the ASCII space character and horizontal tab are removed. Other
     * whitespace characters are preserved and therefore remain part of the
     * token presented to primitive-channel token resolution.
     *
     * @param token Token whose surrounding horizontal whitespace is removed.
     *
     * @return View of the trimmed portion of @p token.
     */
    std::string_view trim_ascii_horizontal_whitespace(
        std::string_view token
    ) noexcept {
        std::size_t first = 0;

        while (
            first < token.size()
            && (token[first] == ' ' || token[first] == '\t')
        ) {
            ++first;
        }

        std::size_t last = token.size();

        while (
            last > first
            && (token[last - 1] == ' ' || token[last - 1] == '\t')
        ) {
            --last;
        }

        return token.substr(first, last - first);
    }

}  // namespace

AnalysisProduct parse_analysis_product(
    std::string_view expression
) {
    const std::string_view configured_expression =
        trim_ascii_horizontal_whitespace(expression);
    if (configured_expression.empty()) {
        throw std::invalid_argument(
            "analysis-product expression must not be empty"
        );
    }

    if (configured_expression.find(',') != std::string_view::npos) {
        throw std::invalid_argument(
            "analysis-product expression must not contain ','"
        );
    }

    AnalysisProduct product{};

    std::size_t token_begin = 0;

    while (token_begin <= configured_expression.size()) {
        const std::size_t separator =
            configured_expression.find('+', token_begin);

        const std::size_t token_end =
            separator == std::string_view::npos
                ? configured_expression.size()
                : separator;

        const std::string_view raw_token =
            configured_expression.substr(
                token_begin,
                token_end - token_begin
            );

        const std::string_view token =
            trim_ascii_horizontal_whitespace(raw_token);

        if (token.empty()) {
            throw std::invalid_argument(
                "analysis-product expression contains an empty "
                "primitive-channel token"
            );
        }

        const auto channel =
            primitive_channel_from_token(token);

        if (!channel.has_value()) {
            throw std::invalid_argument(
                "analysis-product expression contains an unrecognized "
                "primitive-channel token"
            );
        }

        if (std::find(
                product.primitive_channels.begin(),
                product.primitive_channels.end(),
                *channel
            ) != product.primitive_channels.end()) {
            throw std::invalid_argument(
                "analysis-product expression contains a duplicate "
                "primitive channel"
            );
        }

        product.primitive_channels.push_back(*channel);

        if (separator == std::string_view::npos) {
            break;
        }

        token_begin = separator + 1;
    }

    product.configured_expression = std::string(configured_expression);
    return product;
}

}  // namespace hbt
