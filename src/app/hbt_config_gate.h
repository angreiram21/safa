/**
 * @file hbt_config_gate.h
 * @brief Application-level gating of HBT scientific configuration.
 *
 * This file declares the application-layer operation that coordinates global
 * run controls with conditional loading of HBT scientific configuration.
 */

#ifndef SMASH_AFTERBURNER_ANALYSIS_APP_HBT_CONFIG_GATE_H
#define SMASH_AFTERBURNER_ANALYSIS_APP_HBT_CONFIG_GATE_H

#include "config/run_config.h"
#include "hbt/config/hbt_config.h"

#include <optional>

namespace app {

/**
 * @brief Load the HBT scientific configuration when HBT is enabled.
 *
 * This operation coordinates the resolved global run controls with the
 * HBT-specific scientific configuration loader.
 *
 * When HBT is disabled, the function returns std::nullopt and does not access
 * the configured HBT configuration file, even if config::RunConfig contains an
 * HBT configuration path.
 *
 * When HBT is enabled, the function loads the HBT scientific configuration
 * from config::RunConfig::hbt_config_path and returns the resolved
 * hbt::HBTConfig.
 *
 * @param run_config Resolved global run configuration.
 *
 * @return The resolved HBT configuration when HBT is enabled, otherwise
 *         std::nullopt.
 *
 * @throws std::invalid_argument if HBT is enabled but hbt_config_path is
 *         absent from run_config.
 * @throws std::runtime_error if the enabled HBT configuration file cannot be
 *         read or its YAML document violates the HBT configuration structure
 *         contract.
 * @throws std::invalid_argument if the enabled HBT configuration contains an
 *         invalid scientific configuration expression.
 */
std::optional<hbt::HBTConfig> load_hbt_config_if_enabled(
    const config::RunConfig& run_config
);

}  // namespace app

#endif  // SMASH_AFTERBURNER_ANALYSIS_APP_HBT_CONFIG_GATE_H
