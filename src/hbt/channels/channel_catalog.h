/**
 * @file channel_catalog.h
 * @brief Access to the canonical definitions of primitive HBT channels.
 *
 * The catalogue is the single source of truth for the relationship between a
 * PrimitiveChannelId, its canonical ordered particle species, and its stable
 * ASCII name used by configuration and output layers.
 */

#ifndef HBT_CHANNELS_CHANNEL_CATALOG_H
#define HBT_CHANNELS_CHANNEL_CATALOG_H

#include "hbt/channels/primitive_channel.h"

#include <optional>
#include <string_view>

namespace hbt {

/**
 * @brief Obtain the canonical definition of a primitive HBT channel.
 * @param channel Canonical primitive-channel identifier.
 * @return Immutable canonical definition associated with @p channel.
 * @throws std::invalid_argument If @p channel is not a valid identifier.
 */
const PrimitiveChannel& primitive_channel_definition(
    PrimitiveChannelId channel
);

/**
 * @brief Resolve one exact canonical primitive-channel ASCII name.
 * @param canonical_name Exact stable channel name to resolve.
 * @return Matching channel identifier, or std::nullopt when unknown.
 *
 * The name is matched exactly. No trimming, normalization, aliasing, or case
 * conversion is performed.
 */
std::optional<PrimitiveChannelId> primitive_channel_id_from_name(
    std::string_view canonical_name
) noexcept;

}  // namespace hbt

#endif  // HBT_CHANNELS_CHANNEL_CATALOG_H
