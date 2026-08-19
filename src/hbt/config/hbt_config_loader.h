/**
 * @file hbt_config_loader.h
 * @brief Loading of the scientific HBT module configuration from a YAML file.
 *
 * This file declares the configuration boundary that converts an HBT module
 * configuration file into the resolved HBTConfig used by the scientific
 * module.
 *
 * The public interface is independent of the YAML implementation library.
 */

#ifndef HBT_CONFIG_HBT_CONFIG_LOADER_H
#define HBT_CONFIG_HBT_CONFIG_LOADER_H

#include "hbt/config/hbt_config.h"

#include <filesystem>

namespace hbt {

    /**
     * @brief Load and resolve the scientific HBT configuration from YAML.
     *
     * The configuration file must provide the HBT-specific scientific settings
     * required to construct a complete HBTConfig:
     *
     *     hbt_enabled_channels
     *     hbt_particle_acceptance
     *     hbt_pair_slicing
     *     hbt_histograms
     *     hbt_origin_mode
     *
     * hbt_enabled_channels is resolved into the HBTSelection stored in
     * HBTConfig.
     *
     * hbt_particle_acceptance provides the longitudinal-variable choice and the
     * particle-level kinematic acceptance cuts for every configured particle
     * group.
     *
     * hbt_pair_slicing provides the explicit optional kT/mT routing axes.
     *
     * hbt_histograms provides explicit OSL, radial, and relative-time raw
     * histogram binning. No scientific histogram defaults are supplied.
     *
     * hbt_origin_mode provides the requested nested HBT origin-selection mode.
     *
     * Global run controls, module activation, input/output paths, resource
     * settings, and module-configuration file locations are not interpreted by
     * this operation.
     *
     * The loader interface does not expose YAML-library-specific types.
     *
     * @param path Path to the HBT module YAML configuration file.
     *
     * @return Fully resolved scientific HBT configuration.
     *
     * @throws std::runtime_error if the configuration file cannot be read, the
     *         YAML document cannot be interpreted, or a required configuration
     *         entry is missing or has an invalid YAML representation.
     * @throws std::invalid_argument if a textual or numeric scientific value
     *         violates the HBT configuration contract.
     */
    HBTConfig load_hbt_config(
        const std::filesystem::path& path
    );

}  // namespace hbt

#endif  // HBT_CONFIG_HBT_CONFIG_LOADER_H
