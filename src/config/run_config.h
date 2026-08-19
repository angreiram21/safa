/**
 * @file run_config.h
 * @brief Resolved global configuration for one analysis run.
 *
 * This file defines the global configuration data required to assemble an
 * analysis run after the main configuration file has been loaded and resolved.
 *
 * Scientific configuration belonging to individual analysis modules is not
 * stored here.
 */

#ifndef CONFIG_RUN_CONFIG_H
#define CONFIG_RUN_CONFIG_H

#include <cstddef>
#include <filesystem>
#include <optional>

namespace config {

    /**
     * @brief Global configuration required to assemble one analysis run.
     *
     * This type contains run-level controls rather than module-specific
     * scientific settings.
     *
     * events_path is the resolved path to the directory containing the outer
     * event directories used by the analysis run.
     *
     * output_path is the resolved root directory for production output. It is
     * a run-level filesystem setting, not a scientific configuration value.
     *
     * number_of_events is the number of events processed by the analysis run.
     * It must be strictly positive.
     *
     * number_of_subevents is the number of subevents processed for each event.
     * It is a per-event quantity, not a total number of subevents for the
     * complete run, and it must be strictly positive.
     *
     * threads controls Phase-8 outer-event parallelism. A value of zero
     * requests automatic hardware concurrency, one preserves serial event
     * processing, and values greater than one request that many workers at
     * most. The effective worker count never exceeds number_of_events. This
     * control changes execution only; it does not alter scientific defaults or
     * post-sample analysis.
     *
     * The HBT configuration path, when present, is the resolved path to the
     * HBT module configuration file. The run-configuration loader must
     * guarantee that hbt_config_path contains a value whenever hbt_enabled is
     * true.
     *
     * When HBT is disabled, the HBT configuration file must not be loaded or
     * parsed.
     */
    struct RunConfig {
        /// Resolved directory containing the outer-event directories.
        std::filesystem::path events_path;
        /// Resolved root directory for production output.
        std::filesystem::path output_path;
        /// Number of outer events processed by the run.
        std::size_t number_of_events;
        /// Number of independent subevents processed per outer event.
        std::size_t number_of_subevents;
        /// Whether the HBT module is enabled for this run.
        bool hbt_enabled;
        /// Resolved HBT configuration path when HBT is enabled.
        std::optional<std::filesystem::path> hbt_config_path;
        /// Requested Phase-8 outer-event worker count; zero means automatic.
        std::size_t threads{1U};
    };

}  // namespace config

#endif  // CONFIG_RUN_CONFIG_H
