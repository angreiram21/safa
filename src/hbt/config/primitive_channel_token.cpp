/**
 * @file primitive_channel_token.cpp
 * @brief Conversion of HBT configuration tokens to primitive channel IDs.
 */

#include "hbt/config/primitive_channel_token.h"

#include "hbt/channels/channel_catalog.h"

namespace hbt {

std::optional<PrimitiveChannelId> primitive_channel_from_token(
    std::string_view token
) noexcept {
    return primitive_channel_id_from_name(token);
}

}  // namespace hbt
