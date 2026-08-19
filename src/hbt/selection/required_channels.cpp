/**
 * @file required_channels.cpp
 * @brief Derivation of primitive HBT channels required by a selection.
 */

#include "hbt/selection/required_channels.h"

#include <algorithm>

namespace hbt {

    std::vector<PrimitiveChannelId> required_primitive_channels(
        const HBTSelection& selection
    ) {
        std::vector<PrimitiveChannelId> required_channels;

        for (const AnalysisProduct& product : selection.products) {
            for (const PrimitiveChannelId channel :
                 product.primitive_channels) {
                const auto existing_channel = std::find(
                    required_channels.begin(),
                    required_channels.end(),
                    channel
                );

                if (existing_channel == required_channels.end()) {
                    required_channels.push_back(channel);
                }
            }
        }

        return required_channels;
    }

}  // namespace hbt
