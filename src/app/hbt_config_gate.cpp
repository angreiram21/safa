/**
 * @file hbt_config_gate.cpp
 * @brief Application-level gating of HBT scientific configuration.
 *
 * This file implements the application-layer operation that coordinates global
 * run controls with conditional loading of HBT scientific configuration.
 */

#include "app/hbt_config_gate.h"

#include "hbt/config/hbt_config_loader.h"

#include <stdexcept>

namespace app {

/*
 * @brief Load the HBT scientific configuration when HBT is enabled.
 *
 * When HBT is disabled, no HBT configuration file is accessed.
 *
 * When HBT is enabled, the resolved HBT configuration path is required and
 * loading is delegated to the HBT scientific configuration loader.
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
) {
    if (!run_config.hbt_enabled) {
        return std::nullopt;
    }

    if (!run_config.hbt_config_path.has_value()) {
        throw std::invalid_argument(
            "HBT is enabled but RunConfig does not contain hbt_config_path"
        );
    }

    return hbt::load_hbt_config(
        run_config.hbt_config_path.value()
    );
}

}  // namespace app
