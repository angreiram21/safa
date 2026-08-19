/**
 * @file hbt_startup_builder.cpp
 * @brief Construction of resolved HBT startup state.
 */

#include "hbt/startup/hbt_startup_builder.h"

#include "hbt/channels/channel_catalog.h"
#include "hbt/selection/required_channels.h"
#include "hbt/selection/required_species.h"

#include <algorithm>
#include <stdexcept>
#include <utility>
#include <vector>

namespace hbt {
namespace {

/**
 * @brief Return an order-independent channel signature for one product.
 * @param product Product whose primitive channels are already validated.
 * @return Sorted primitive-channel identifiers preserving multiplicity.
 */
std::vector<PrimitiveChannelId> canonical_product_signature(
    const AnalysisProduct& product
) {
    std::vector<PrimitiveChannelId> signature = product.primitive_channels;
    std::sort(signature.begin(), signature.end());
    return signature;
}

/**
 * @brief Validate final-product structure before deriving startup requirements.
 * @param selection Parsed HBT final-product selection.
 * @throws std::invalid_argument If the selection is empty, one product is
 *         empty, one product repeats a primitive channel, a channel identifier
 *         is invalid, or two products are semantically identical.
 *
 * Duplicate channels and products are rejected rather than deduplicated. A
 * repeated channel would otherwise assign fictitious statistical weight to
 * the same physical pairs.
 */
void validate_selection(const HBTSelection& selection) {
    if (selection.products.empty()) {
        throw std::invalid_argument(
            "HBT startup: selection must contain at least one final product"
        );
    }

    std::vector<std::vector<PrimitiveChannelId>> signatures;
    signatures.reserve(selection.products.size());

    for (const AnalysisProduct& product : selection.products) {
        if (product.primitive_channels.empty()) {
            throw std::invalid_argument(
                "HBT startup: final product must contain at least one channel"
            );
        }

        for (std::size_t index = 0U;
             index < product.primitive_channels.size();
             ++index) {
            const PrimitiveChannelId channel =
                product.primitive_channels[index];
            static_cast<void>(primitive_channel_definition(channel));

            for (std::size_t previous = 0U;
                 previous < index;
                 ++previous) {
                if (product.primitive_channels[previous] == channel) {
                    throw std::invalid_argument(
                        "HBT startup: duplicate primitive channel in product"
                    );
                }
            }
        }

        std::vector<PrimitiveChannelId> signature =
            canonical_product_signature(product);
        if (std::find(signatures.begin(), signatures.end(), signature) !=
            signatures.end()) {
            throw std::invalid_argument(
                "HBT startup: duplicate final analysis product"
            );
        }
        signatures.push_back(std::move(signature));
    }
}

}  // namespace

HBTStartupState build_hbt_startup_state(
    const HBTConfig& config
) {
    validate_selection(config.selection);

    std::vector<PrimitiveChannelId> required_channels =
        required_primitive_channels(config.selection);
    std::vector<SpeciesId> required_particle_species =
        required_species(required_channels);

    return HBTStartupState{
        config.selection,
        std::move(required_channels),
        std::move(required_particle_species)
    };
}

}  // namespace hbt
