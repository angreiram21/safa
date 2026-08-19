/**
 * @file run_config_loader.h
 * @brief Loading of the resolved global run configuration from a YAML file.
 *
 * This file declares the configuration boundary that converts the global
 * run-configuration YAML file into the resolved RunConfig used to assemble
 * one analysis run.
 *
 * The public interface is independent of the YAML implementation library.
 */

#ifndef CONFIG_RUN_CONFIG_LOADER_H
#define CONFIG_RUN_CONFIG_LOADER_H

#include "config/run_config.h"

#include <filesystem>

namespace config {

    /**
     * @brief Load and resolve the global configuration for one analysis run.
     *
     * The configuration file provides run-level controls such as the event
     * source location, the number of events, the number of independent
     * subevents associated with each event, analysis-module activation, and
     * paths to module-specific configuration files.
     *
     * At the current stage, the recognized entries are events_path,
     * output_path, number_of_events, number_of_subevents, hbt_enabled, and
     * hbt_config_path.
     *
     * events_path is required and must be a non-empty scalar. It identifies the
     * directory containing the outer event directories. A relative events_path
     * is resolved relative to the directory containing the run-configuration
     * file. An absolute path is preserved. This operation does not require the
     * resolved directory to exist.
     *
     * output_path is required and must be a non-empty scalar. It identifies
     * the production-output root. A relative output_path is resolved relative
     * to the directory containing the run-configuration file. An absolute
     * path is preserved. This operation does not create or require the
     * resolved directory to exist.
     *
     * number_of_events is required and must be a scalar representing a strictly
     * positive integer that can be represented by std::size_t.
     *
     * number_of_subevents is required and must be a scalar representing a
     * strictly positive integer that can be represented by std::size_t. It
     * specifies the number of independent subevents associated with each event,
     * not the total number of subevents in the complete run.
     *
     * hbt_enabled is required and determines whether the HBT module is active.
     *
     * When hbt_enabled is true, hbt_config_path is required. A relative
     * hbt_config_path is resolved relative to the directory containing the
     * run-configuration file. An absolute path is preserved.
     *
     * When hbt_enabled is false, hbt_config_path is optional. If present, it is
     * resolved in the same way, but the HBT configuration file is not loaded or
     * parsed by this operation.
     *
     * This operation resolves only the global run configuration. It does not
     * enumerate event directories, open event input files, load or interpret
     * module-specific scientific configuration files, or perform event or
     * subevent processing.
     *
     * The loader interface does not expose YAML-library-specific types.
     *
     * @param path Path to the global run-configuration YAML file.
     *
     * @return Fully resolved global run configuration.
     *
     * @throws std::runtime_error if the configuration file cannot be read, the
     *         YAML document cannot be interpreted, the document violates the
     *         global configuration structure, or a required configuration entry
     *         is missing or has an invalid YAML representation or value.
     */
    RunConfig load_run_config(
        const std::filesystem::path& path
    );

}  // namespace config

#endif  // CONFIG_RUN_CONFIG_LOADER_H
