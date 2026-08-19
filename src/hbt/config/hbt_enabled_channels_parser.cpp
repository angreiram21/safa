/**
 * @file hbt_enabled_channels_parser.cpp
 * @brief Parsing of the complete HBT enabled-channels configuration expression.
 */

#include "hbt/config/hbt_enabled_channels_parser.h"

#include "hbt/config/analysis_product_parser.h"

#include <algorithm>
#include <set>
#include <stdexcept>
#include <string_view>
#include <utility>
#include <vector>

namespace hbt {

HBTSelection parse_hbt_enabled_channels(
    std::string_view expression
) {
    if (expression.empty()) {
        throw std::invalid_argument(
            "hbt_enabled_channels expression must not be empty"
        );
    }

    HBTSelection selection{};
    std::set<std::vector<PrimitiveChannelId>> canonical_products;

    std::size_t product_begin = 0;

    while (product_begin <= expression.size()) {
        const std::size_t separator =
            expression.find(',', product_begin);

        const std::size_t product_end =
            separator == std::string_view::npos
                ? expression.size()
                : separator;

        const std::string_view product_expression =
            expression.substr(
                product_begin,
                product_end - product_begin
            );

        AnalysisProduct product =
            parse_analysis_product(product_expression);
        std::vector<PrimitiveChannelId> canonical =
            product.primitive_channels;
        std::sort(canonical.begin(), canonical.end());
        if (!canonical_products.insert(canonical).second) {
            throw std::invalid_argument(
                "hbt_enabled_channels contains a duplicate final product"
            );
        }
        selection.products.push_back(std::move(product));

        if (separator == std::string_view::npos) {
            break;
        }

        product_begin = separator + 1;
    }

    return selection;
}

}  // namespace hbt
