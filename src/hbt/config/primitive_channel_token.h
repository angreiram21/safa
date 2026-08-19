/**
 * @file primitive_channel_token.h
 * @brief Conversion of HBT configuration tokens to primitive channel IDs.
 *
 * This file declares the configuration-facing operation used to resolve an
 * exact primitive-channel token to its canonical PrimitiveChannelId.
 *
 * Tokenization of complete configuration expressions, product composition,
 * whitespace handling, YAML parsing, and HBT selection construction are
 * outside this operation's responsibility.
 */

#ifndef HBT_CONFIG_PRIMITIVE_CHANNEL_TOKEN_H
#define HBT_CONFIG_PRIMITIVE_CHANNEL_TOKEN_H

#include "hbt/channels/primitive_channel.h"

#include <optional>
#include <string_view>

namespace hbt {

    /**
     * @brief Resolve an exact primitive-channel configuration token.
     *
     * The function performs an exact match against the canonical textual tokens
     * accepted for primitive HBT channels.
     *
     * The token is not trimmed, normalized, or modified. Handling separators,
     * surrounding whitespace, and complete configuration expressions belongs
     * to the higher-level configuration parser.
     *
     * @param token Exact primitive-channel token to resolve.
     *
     * @return The corresponding PrimitiveChannelId when @p token is recognized;
     *         std::nullopt otherwise.
     */
    std::optional<PrimitiveChannelId> primitive_channel_from_token(
        std::string_view token
    ) noexcept;

}  // namespace hbt

#endif  // HBT_CONFIG_PRIMITIVE_CHANNEL_TOKEN_H
