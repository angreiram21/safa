/**
 * @file run_config_loader.cpp
 * @brief Loading of the resolved global run configuration from a YAML file.
 */

#include "config/run_config_loader.h"

#include <yaml-cpp/yaml.h>

#include <charconv>
#include <cstddef>
#include <filesystem>
#include <optional>
#include <stdexcept>
#include <string>
#include <string_view>
#include <system_error>

namespace config {

namespace {

    /**
     * @brief Intermediate values extracted from the global YAML mapping.
     *
     * Optional members distinguish an absent configuration entry from an entry
     * whose parsed value is otherwise valid.
     */
    struct RunConfigEntries {
        /// Raw configured events path when present.
        std::optional<std::string> events_path;
        /// Raw configured production-output path when present.
        std::optional<std::string> output_path;
        /// Raw configured outer-event count when present.
        std::optional<std::size_t> number_of_events;
        /// Raw configured subevent count when present.
        std::optional<std::size_t> number_of_subevents;
        /// Raw configured Phase-8 worker count when present.
        std::optional<std::size_t> threads;
        /// Raw configured HBT activation flag when present.
        std::optional<bool> hbt_enabled;
        /// Raw configured HBT configuration path when present.
        std::optional<std::string> hbt_config_path;
    };

    /**
     * @brief Load one YAML document from the global run-configuration file.
     *
     * YAML-library exceptions produced while opening or parsing the file are
     * translated into std::runtime_error so that YAML-specific exception types
     * do not escape through the public configuration interface.
     *
     * @param path Path to the global run-configuration YAML file.
     *
     * @return Parsed YAML root node.
     *
     * @throws std::runtime_error if the file cannot be read or the YAML
     *         document cannot be parsed.
     */
    YAML::Node load_yaml_document(
        const std::filesystem::path& path
    ) {
        try {
            return YAML::LoadFile(path.string());
        } catch (const YAML::Exception& exception) {
            throw std::runtime_error(
                "failed to load run configuration file '" +
                path.string() +
                "': " +
                exception.what()
            );
        }
    }

    /**
     * @brief Parse one strictly positive std::size_t configuration value.
     *
     * The YAML node must be a scalar whose complete textual representation is
     * a base-10 unsigned integer greater than zero. Leading signs, negative
     * values, fractional values, surrounding whitespace, trailing characters,
     * and values outside the range of std::size_t are rejected.
     *
     * @param node YAML node containing the configured value.
     * @param entry_name Name of the global configuration entry.
     *
     * @return Parsed strictly positive integer.
     *
     * @throws std::runtime_error if the value is not a valid strictly positive
     *         integer representable by std::size_t.
     */
    std::size_t parse_positive_size_t(
        const YAML::Node& node,
        std::string_view entry_name
    ) {
        const std::string name{entry_name};

        if (!node.IsScalar()) {
            throw std::runtime_error(
                "run configuration entry '" +
                name +
                "' must be a YAML scalar containing a strictly positive integer"
            );
        }

        const std::string text = node.Scalar();

        if (text.empty()) {
            throw std::runtime_error(
                "run configuration entry '" +
                name +
                "' must contain a strictly positive integer"
            );
        }

        std::size_t value = 0;

        const char* const begin = text.data();
        const char* const end = begin + text.size();

        const auto result =
            std::from_chars(begin, end, value, 10);

        if (
            result.ec != std::errc{} ||
            result.ptr != end ||
            value == 0
        ) {
            throw std::runtime_error(
                "run configuration entry '" +
                name +
                "' must contain a strictly positive integer "
                "representable by std::size_t"
            );
        }

        return value;
    }

    /**
     * @brief Parse one non-negative std::size_t configuration value.
     *
     * The complete scalar must be a base-10 unsigned integer representable by
     * std::size_t. Zero is accepted for controls whose contract defines it.
     *
     * @param node YAML node containing the configured value.
     * @param entry_name Name of the global configuration entry.
     * @return Parsed non-negative integer.
     * @throws std::runtime_error If the value is not a valid std::size_t.
     */
    std::size_t parse_nonnegative_size_t(
        const YAML::Node& node,
        std::string_view entry_name
    ) {
        const std::string name{entry_name};

        if (!node.IsScalar()) {
            throw std::runtime_error(
                "run configuration entry '" + name +
                "' must be a YAML scalar containing a non-negative integer"
            );
        }

        const std::string text = node.Scalar();
        std::size_t value = 0U;
        const char* const begin = text.data();
        const char* const end = begin + text.size();
        const auto result = std::from_chars(begin, end, value, 10);

        if (text.empty() || result.ec != std::errc{} || result.ptr != end) {
            throw std::runtime_error(
                "run configuration entry '" + name +
                "' must contain a non-negative integer representable by "
                "std::size_t"
            );
        }

        return value;
    }

    /**
     * @brief Parse the hbt_enabled value from one YAML node.
     *
     * @param node YAML node containing the configured activation value.
     *
     * @return Parsed HBT activation state.
     *
     * @throws std::runtime_error if the value is not a scalar that can be
     *         interpreted as a boolean.
     */
    bool parse_hbt_enabled(
        const YAML::Node& node
    ) {
        if (!node.IsScalar()) {
            throw std::runtime_error(
                "run configuration entry 'hbt_enabled' "
                "must be a YAML boolean scalar"
            );
        }

        try {
            return node.as<bool>();
        } catch (const YAML::Exception&) {
            throw std::runtime_error(
                "run configuration entry 'hbt_enabled' "
                "must be a YAML boolean scalar"
            );
        }
    }

    /**
     * @brief Parse one required or optional non-empty configured path.
     *
     * @param node YAML node containing the configured path.
     * @param entry_name Name of the global configuration entry.
     *
     * @return Non-empty textual path exactly as represented by the scalar.
     *
     * @throws std::runtime_error if the value is not a scalar or the resulting
     *         path string is empty.
     */
    std::string parse_nonempty_path(
        const YAML::Node& node,
        std::string_view entry_name
    ) {
        const std::string name{entry_name};

        if (!node.IsScalar()) {
            throw std::runtime_error(
                "run configuration entry '" +
                name +
                "' must be a YAML scalar"
            );
        }

        const std::string configured_path = node.Scalar();

        if (configured_path.empty()) {
            throw std::runtime_error(
                "run configuration entry '" +
                name +
                "' must not be empty"
            );
        }

        return configured_path;
    }

    /**
     * @brief Extract the recognized entries from the global YAML root mapping.
     *
     * The current global configuration contract recognizes events_path,
     * output_path, number_of_events, number_of_subevents, threads,
     * hbt_enabled, and hbt_config_path.
     * Unknown keys, duplicate occurrences, and non-scalar keys are rejected.
     *
     * Required-entry and conditional-entry relationships are validated after
     * extraction when the resolved RunConfig is constructed.
     *
     * @param root Parsed root node of the global YAML document.
     *
     * @return Parsed entries while preserving whether each key was present.
     *
     * @throws std::runtime_error if the document does not satisfy the current
     *         global YAML structure contract.
     */
    RunConfigEntries parse_run_config_entries(
        const YAML::Node& root
    ) {
        if (!root.IsMap()) {
            throw std::runtime_error(
                "run configuration root must be a YAML mapping"
            );
        }

        RunConfigEntries entries{};

        for (const auto& entry : root) {
            if (!entry.first.IsScalar()) {
                throw std::runtime_error(
                    "run configuration keys must be YAML scalars"
                );
            }

            const std::string key = entry.first.Scalar();

            if (key == "events_path") {
                if (entries.events_path.has_value()) {
                    throw std::runtime_error(
                        "run configuration contains duplicate key "
                        "'events_path'"
                    );
                }

                entries.events_path =
                    parse_nonempty_path(
                        entry.second,
                        "events_path"
                    );

                continue;
            }

            if (key == "output_path") {
                if (entries.output_path.has_value()) {
                    throw std::runtime_error(
                        "run configuration contains duplicate key "
                        "'output_path'"
                    );
                }

                entries.output_path =
                    parse_nonempty_path(
                        entry.second,
                        "output_path"
                    );

                continue;
            }

            if (key == "number_of_events") {
                if (entries.number_of_events.has_value()) {
                    throw std::runtime_error(
                        "run configuration contains duplicate key "
                        "'number_of_events'"
                    );
                }

                entries.number_of_events =
                    parse_positive_size_t(
                        entry.second,
                        "number_of_events"
                    );

                continue;
            }

            if (key == "number_of_subevents") {
                if (entries.number_of_subevents.has_value()) {
                    throw std::runtime_error(
                        "run configuration contains duplicate key "
                        "'number_of_subevents'"
                    );
                }

                entries.number_of_subevents =
                    parse_positive_size_t(
                        entry.second,
                        "number_of_subevents"
                    );

                continue;
            }

            if (key == "threads") {
                if (entries.threads.has_value()) {
                    throw std::runtime_error(
                        "run configuration contains duplicate key 'threads'"
                    );
                }

                entries.threads = parse_nonnegative_size_t(
                    entry.second,
                    "threads"
                );
                continue;
            }

            if (key == "hbt_enabled") {
                if (entries.hbt_enabled.has_value()) {
                    throw std::runtime_error(
                        "run configuration contains duplicate key "
                        "'hbt_enabled'"
                    );
                }

                entries.hbt_enabled =
                    parse_hbt_enabled(entry.second);

                continue;
            }

            if (key == "hbt_config_path") {
                if (entries.hbt_config_path.has_value()) {
                    throw std::runtime_error(
                        "run configuration contains duplicate key "
                        "'hbt_config_path'"
                    );
                }

                entries.hbt_config_path =
                    parse_nonempty_path(
                        entry.second,
                        "hbt_config_path"
                    );

                continue;
            }

            throw std::runtime_error(
                "run configuration contains unknown key '" +
                key +
                "'"
            );
        }

        return entries;
    }

    /**
     * @brief Resolve one configured path relative to the run configuration.
     *
     * Absolute configured paths are returned unchanged. Relative configured
     * paths are interpreted relative to the directory containing the global
     * run-configuration file.
     *
     * This operation performs lexical path composition only. It does not open
     * the target path or require the target to exist.
     *
     * @param run_config_path Path to the global run-configuration file.
     * @param configured_path Path obtained from the YAML document.
     *
     * @return Resolved filesystem path.
     */
    std::filesystem::path resolve_configured_path(
        const std::filesystem::path& run_config_path,
        const std::string& configured_path
    ) {
        const std::filesystem::path path{configured_path};

        if (path.is_absolute()) {
            return path;
        }

        return run_config_path.parent_path() / path;
    }

    /**
     * @brief Construct the resolved RunConfig from parsed global entries.
     *
     * events_path, output_path, number_of_events, number_of_subevents, and
     * hbt_enabled are required. threads is optional and defaults to one.
     * hbt_config_path is additionally required when
     * HBT is enabled.
     * If an HBT path is provided while HBT is disabled, it remains available
     * in RunConfig but no HBT configuration file is opened here.
     *
     * events_path and output_path are resolved but are not opened, enumerated,
     * created, or required to exist by this operation.
     *
     * number_of_subevents is stored as the number of independent subevents
     * associated with each event. No total subevent count is calculated by this
     * operation.
     *
     * @param entries Parsed global configuration entries.
     * @param run_config_path Path to the global run-configuration file.
     *
     * @return Resolved global run configuration.
     *
     * @throws std::runtime_error if a required global entry is absent or if HBT
     *         is enabled without hbt_config_path.
     */
    RunConfig resolve_run_config(
        const RunConfigEntries& entries,
        const std::filesystem::path& run_config_path
    ) {
        if (!entries.events_path.has_value()) {
            throw std::runtime_error(
                "run configuration is missing required entry "
                "'events_path'"
            );
        }

        if (!entries.output_path.has_value()) {
            throw std::runtime_error(
                "run configuration is missing required entry "
                "'output_path'"
            );
        }

        if (!entries.number_of_events.has_value()) {
            throw std::runtime_error(
                "run configuration is missing required entry "
                "'number_of_events'"
            );
        }

        if (!entries.number_of_subevents.has_value()) {
            throw std::runtime_error(
                "run configuration is missing required entry "
                "'number_of_subevents'"
            );
        }

        if (!entries.hbt_enabled.has_value()) {
            throw std::runtime_error(
                "run configuration is missing required entry "
                "'hbt_enabled'"
            );
        }

        if (
            entries.hbt_enabled.value() &&
            !entries.hbt_config_path.has_value()
        ) {
            throw std::runtime_error(
                "run configuration is missing required entry "
                "'hbt_config_path' while HBT is enabled"
            );
        }

        const std::filesystem::path resolved_events_path =
            resolve_configured_path(
                run_config_path,
                entries.events_path.value()
            );
        const std::filesystem::path resolved_output_path =
            resolve_configured_path(
                run_config_path,
                entries.output_path.value()
            );

        std::optional<std::filesystem::path> resolved_hbt_config_path;

        if (entries.hbt_config_path.has_value()) {
            resolved_hbt_config_path =
                resolve_configured_path(
                    run_config_path,
                    entries.hbt_config_path.value()
                );
        }

        return RunConfig{
            resolved_events_path,
            resolved_output_path,
            entries.number_of_events.value(),
            entries.number_of_subevents.value(),
            entries.hbt_enabled.value(),
            resolved_hbt_config_path,
            entries.threads.value_or(1U)
        };
    }

}  // namespace

RunConfig load_run_config(
    const std::filesystem::path& path
) {
    const YAML::Node root = load_yaml_document(path);

    const RunConfigEntries entries =
        parse_run_config_entries(root);

    return resolve_run_config(entries, path);
}

}  // namespace config
