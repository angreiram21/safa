/**
 * @file hbt_startup_builder.h
 * @brief Construction of resolved HBT startup state.
 */

#ifndef HBT_STARTUP_HBT_STARTUP_BUILDER_H
#define HBT_STARTUP_HBT_STARTUP_BUILDER_H

#include "hbt/config/hbt_config.h"
#include "hbt/startup/hbt_startup_state.h"

namespace hbt {

/**
 * @brief Build the resolved HBT state required before event processing.
 *
 * The returned state contains the configured HBT selection together with the
 * unique primitive channels and particle species required by that selection.
 *
 * This operation performs no file access, YAML parsing, event reading,
 * particle identification, pair construction, observable calculation, or
 * accumulator allocation.
 *
 * @param config Resolved scientific HBT configuration.
 *
 * @return HBT startup state derived from the configured selection.
 *
 * @throws std::invalid_argument If the selection is empty, contains an empty
 *         product, repeats one primitive channel within a product, contains
 *         semantically duplicate products, or contains an invalid channel.
 */
HBTStartupState build_hbt_startup_state(
    const HBTConfig& config
);

}  // namespace hbt

#endif
