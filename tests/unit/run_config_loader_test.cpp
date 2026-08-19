/**
 * @file run_config_loader_test.cpp
 * @brief Unit tests for loading the resolved global run configuration.
 *
 * This test verifies conversion of the current global YAML configuration
 * contract into RunConfig. It also verifies run-size controls, event-source,
 * production-output, and module path resolution, conditional HBT
 * configuration requirements,
 * structural YAML errors, and file-loading errors.
 */

#include "config/run_config_loader.h"

#include <cstddef>
#include <cstdlib>
#include <exception>
#include <filesystem>
#include <fstream>
#include <iostream>
#include <optional>
#include <stdexcept>
#include <string_view>
#include <system_error>

namespace {

    /**
     * @brief Return the directory used by the run-configuration loader tests.
     *
     * @return Relative path to the test working directory.
     */
    std::filesystem::path test_directory() {
        return "run_config_loader_test_data";
    }

    /**
     * @brief Prepare an empty directory for temporary test configuration files.
     *
     * @return true when the directory is ready for use, otherwise false.
     */
    bool prepare_test_directory() {
        const std::filesystem::path directory = test_directory();

        std::error_code error;

        static_cast<void>(
            std::filesystem::remove_all(directory, error)
        );

        error.clear();

        static_cast<void>(
            std::filesystem::create_directories(directory, error)
        );

        if (error) {
            std::cerr
                << "run_config_loader_test: failed to prepare test directory '"
                << directory.string()
                << "'.\n";
            return false;
        }

        return true;
    }

    /**
     * @brief Remove the temporary test directory without throwing.
     */
    void cleanup_test_directory() noexcept {
        std::error_code error;

        static_cast<void>(
            std::filesystem::remove_all(
                test_directory(),
                error
            )
        );
    }

    /**
     * @brief Write textual YAML content to one temporary test file.
     *
     * @param path Path of the file to create.
     * @param contents Complete YAML document to write.
     *
     * @return true when the file is created successfully, otherwise false.
     */
    bool write_test_file(
        const std::filesystem::path& path,
        std::string_view contents
    ) {
        std::ofstream file(path);

        if (!file) {
            std::cerr
                << "run_config_loader_test: failed to create test file '"
                << path.string()
                << "'.\n";
            return false;
        }

        file.write(
            contents.data(),
            static_cast<std::streamsize>(contents.size())
        );

        if (!file) {
            std::cerr
                << "run_config_loader_test: failed to write test file '"
                << path.string()
                << "'.\n";
            return false;
        }

        return true;
    }

    /**
     * @brief Verify one successfully loaded global run configuration.
     *
     * @param test_name Human-readable test-case name.
     * @param path Path used for the temporary main YAML document.
     * @param contents YAML document to write.
     * @param expected_events_path Expected resolved outer-event directory.
     * @param expected_output_path Expected resolved production-output root.
     * @param expected_number_of_events Expected number of events.
     * @param expected_number_of_subevents Expected subevents per event.
     * @param expected_hbt_enabled Expected HBT activation state.
     * @param expected_hbt_config_path Expected resolved HBT configuration path.
     * @param expected_threads Expected Phase-8 thread control.
     *
     * @return true when the loaded RunConfig matches all expected values,
     *         otherwise false.
     */
    bool verify_valid_config(
        const char* test_name,
        const std::filesystem::path& path,
        std::string_view contents,
        const std::filesystem::path& expected_events_path,
        const std::filesystem::path& expected_output_path,
        std::size_t expected_number_of_events,
        std::size_t expected_number_of_subevents,
        bool expected_hbt_enabled,
        const std::optional<std::filesystem::path>&
            expected_hbt_config_path,
        std::size_t expected_threads = 1U
    ) {
        if (!write_test_file(path, contents)) {
            return false;
        }

        try {
            const config::RunConfig run_config =
                config::load_run_config(path);

            if (run_config.events_path != expected_events_path) {
                std::cerr
                    << "run_config_loader_test: "
                    << test_name
                    << " produced the wrong events_path.\n";
                return false;
            }

            if (run_config.output_path != expected_output_path) {
                std::cerr
                    << "run_config_loader_test: "
                    << test_name
                    << " produced the wrong output_path.\n";
                return false;
            }

            if (
                run_config.number_of_events !=
                expected_number_of_events
            ) {
                std::cerr
                    << "run_config_loader_test: "
                    << test_name
                    << " produced the wrong number_of_events value.\n";
                return false;
            }

            if (
                run_config.number_of_subevents !=
                expected_number_of_subevents
            ) {
                std::cerr
                    << "run_config_loader_test: "
                    << test_name
                    << " produced the wrong number_of_subevents value.\n";
                return false;
            }

            if (run_config.threads != expected_threads) {
                std::cerr
                    << "run_config_loader_test: "
                    << test_name
                    << " produced the wrong threads value.\n";
                return false;
            }

            if (run_config.hbt_enabled != expected_hbt_enabled) {
                std::cerr
                    << "run_config_loader_test: "
                    << test_name
                    << " produced the wrong hbt_enabled value.\n";
                return false;
            }

            if (
                run_config.hbt_config_path !=
                expected_hbt_config_path
            ) {
                std::cerr
                    << "run_config_loader_test: "
                    << test_name
                    << " produced the wrong hbt_config_path.\n";
                return false;
            }

            return true;
        } catch (const std::exception& exception) {
            std::cerr
                << "run_config_loader_test: "
                << test_name
                << " unexpectedly threw: "
                << exception.what()
                << ".\n";
            return false;
        }
    }

    /**
     * @brief Verify that one YAML document produces std::runtime_error.
     *
     * @param test_name Human-readable test-case name.
     * @param path Path used for the temporary main YAML document.
     * @param contents YAML document to write.
     *
     * @return true when std::runtime_error is thrown, otherwise false.
     */
    bool verify_runtime_error(
        const char* test_name,
        const std::filesystem::path& path,
        std::string_view contents
    ) {
        if (!write_test_file(path, contents)) {
            return false;
        }

        try {
            static_cast<void>(
                config::load_run_config(path)
            );
        } catch (const std::runtime_error&) {
            return true;
        } catch (const std::exception& exception) {
            std::cerr
                << "run_config_loader_test: "
                << test_name
                << " threw an unexpected exception type: "
                << exception.what()
                << ".\n";
            return false;
        }

        std::cerr
            << "run_config_loader_test: "
            << test_name
            << " did not throw std::runtime_error.\n";
        return false;
    }

    /**
     * @brief Verify resolution of a relative outer-event directory path.
     *
     * The target directory is deliberately not created. Successful loading
     * therefore also verifies that load_run_config() resolves events_path
     * without accessing the configured event source.
     *
     * @return true when the expected resolved events_path is produced.
     */
    bool verify_relative_events_path() {
        const std::filesystem::path main_path =
            test_directory() / "main_relative_events_path.yaml";

        const std::filesystem::path expected_events_path =
            test_directory() / "relative_events";

        return verify_valid_config(
            "relative events_path",
            main_path,
            "events_path: \"relative_events\"\n"
            "output_path: \"output\"\n"
            "number_of_events: 1\n"
            "number_of_subevents: 2\n"
            "hbt_enabled: false\n",
            expected_events_path,
            test_directory() / "output",
            1,
            2,
            false,
            std::nullopt
        );
    }

    /**
     * @brief Verify preservation of an absolute outer-event directory path.
     *
     * @return true when the absolute events_path is returned unchanged.
     */
    bool verify_absolute_events_path() {
        const std::filesystem::path main_path =
            test_directory() / "main_absolute_events_path.yaml";

        const std::filesystem::path expected_events_path =
            "/run_config_loader_test_absolute_events";

        return verify_valid_config(
            "absolute events_path",
            main_path,
            "events_path: \"/run_config_loader_test_absolute_events\"\n"
            "output_path: \"output\"\n"
            "number_of_events: 1\n"
            "number_of_subevents: 2\n"
            "hbt_enabled: false\n",
            expected_events_path,
            test_directory() / "output",
            1,
            2,
            false,
            std::nullopt
        );
    }

    /**
     * @brief Verify preservation of an absolute production-output path.
     *
     * @return true when the absolute output_path is returned unchanged.
     */
    bool verify_absolute_output_path() {
        const std::filesystem::path expected_output_path =
            "/run_config_loader_test_absolute_output";

        return verify_valid_config(
            "absolute output_path",
            test_directory() / "main_absolute_output_path.yaml",
            "events_path: \"events\"\n"
            "output_path: \"/run_config_loader_test_absolute_output\"\n"
            "number_of_events: 1\n"
            "number_of_subevents: 2\n"
            "hbt_enabled: false\n",
            test_directory() / "events",
            expected_output_path,
            1,
            2,
            false,
            std::nullopt
        );
    }

    /**
     * @brief Verify enabled HBT with a relative module-configuration path.
     *
     * The target hbt.yaml file is deliberately not created. Successful loading
     * therefore also verifies that load_run_config() resolves the path without
     * opening the module-specific configuration file.
     *
     * @return true when the expected RunConfig is produced.
     */
    bool verify_enabled_relative_path() {
        const std::filesystem::path main_path =
            test_directory() / "main_enabled_relative.yaml";

        const std::filesystem::path expected_hbt_path =
            test_directory() / "hbt.yaml";

        return verify_valid_config(
            "enabled HBT with relative path",
            main_path,
            "events_path: \"events\"\n"
            "output_path: \"output\"\n"
            "number_of_events: 500\n"
            "number_of_subevents: 2000\n"
            "hbt_enabled: true\n"
            "hbt_config_path: \"hbt.yaml\"\n",
            test_directory() / "events",
            test_directory() / "output",
            500,
            2000,
            true,
            std::optional<std::filesystem::path>{
                expected_hbt_path
            }
        );
    }

    /**
     * @brief Verify preservation of an absolute HBT configuration path.
     *
     * @return true when the absolute path is returned unchanged.
     */
    bool verify_enabled_absolute_path() {
        const std::filesystem::path main_path =
            test_directory() / "main_enabled_absolute.yaml";

        const std::filesystem::path expected_hbt_path =
            "/run_config_loader_test_absolute_hbt.yaml";

        return verify_valid_config(
            "enabled HBT with absolute path",
            main_path,
            "events_path: \"events\"\n"
            "output_path: \"output\"\n"
            "number_of_events: 3\n"
            "number_of_subevents: 4\n"
            "hbt_enabled: true\n"
            "hbt_config_path: "
            "\"/run_config_loader_test_absolute_hbt.yaml\"\n",
            test_directory() / "events",
            test_directory() / "output",
            3,
            4,
            true,
            std::optional<std::filesystem::path>{
                expected_hbt_path
            }
        );
    }

    /**
     * @brief Verify that disabled HBT does not require a configuration path.
     *
     * @return true when HBT is disabled and no path is stored.
     */
    bool verify_disabled_without_path() {
        return verify_valid_config(
            "disabled HBT without path",
            test_directory() / "main_disabled_without_path.yaml",
            "events_path: \"events\"\n"
            "output_path: \"output\"\n"
            "number_of_events: 2\n"
            "number_of_subevents: 5\n"
            "hbt_enabled: false\n",
            test_directory() / "events",
            test_directory() / "output",
            2,
            5,
            false,
            std::nullopt
        );
    }

    /**
     * @brief Verify that a path may remain present while HBT is disabled.
     *
     * The target file is deliberately absent. Successful loading verifies that
     * the global loader does not attempt to open module-specific configuration
     * files even when a path is present.
     *
     * @return true when the path is resolved and retained without being loaded.
     */
    bool verify_disabled_with_path() {
        const std::filesystem::path expected_hbt_path =
            test_directory() / "disabled_hbt.yaml";

        return verify_valid_config(
            "disabled HBT with path",
            test_directory() / "main_disabled_with_path.yaml",
            "events_path: \"events\"\n"
            "output_path: \"output\"\n"
            "number_of_events: 7\n"
            "number_of_subevents: 11\n"
            "hbt_enabled: false\n"
            "hbt_config_path: \"disabled_hbt.yaml\"\n",
            test_directory() / "events",
            test_directory() / "output",
            7,
            11,
            false,
            std::optional<std::filesystem::path>{
                expected_hbt_path
            }
        );
    }

    /**
     * @brief Verify rejection of a configuration file that does not exist.
     *
     * @return true when std::runtime_error is thrown.
     */
    bool verify_missing_main_file() {
        const std::filesystem::path path =
            test_directory() / "missing_main.yaml";

        std::error_code error;

        static_cast<void>(
            std::filesystem::remove(path, error)
        );

        try {
            static_cast<void>(
                config::load_run_config(path)
            );
        } catch (const std::runtime_error&) {
            return true;
        } catch (const std::exception& exception) {
            std::cerr
                << "run_config_loader_test: missing main file threw an "
                << "unexpected exception type: "
                << exception.what()
                << ".\n";
            return false;
        }

        std::cerr
            << "run_config_loader_test: missing main file did not throw "
            << "std::runtime_error.\n";
        return false;
    }

    /**
     * @brief Verify rejection of syntactically invalid YAML.
     *
     * @return true when std::runtime_error is thrown.
     */
    bool verify_invalid_yaml() {
        return verify_runtime_error(
            "invalid YAML",
            test_directory() / "main_invalid_yaml.yaml",
            "events_path: \"events\"\n"
            "output_path: \"output\"\n"
            "number_of_events: 1\n"
            "number_of_subevents: 1\n"
            "hbt_enabled: [\n"
        );
    }

    /**
     * @brief Verify rejection of a non-mapping YAML root.
     *
     * @return true when std::runtime_error is thrown.
     */
    bool verify_non_mapping_root() {
        return verify_runtime_error(
            "non-mapping root",
            test_directory() / "main_non_mapping_root.yaml",
            "- number_of_events\n"
            "- 1\n"
            "- number_of_subevents\n"
            "- 1\n"
            "- hbt_enabled\n"
            "- true\n"
        );
    }

    /**
     * @brief Verify that events_path is required.
     *
     * @return true when std::runtime_error is thrown.
     */
    bool verify_missing_events_path() {
        return verify_runtime_error(
            "missing events_path",
            test_directory() / "main_missing_events_path.yaml",
            "output_path: \"output\"\n"
            "number_of_events: 1\n"
            "number_of_subevents: 1\n"
            "hbt_enabled: false\n"
        );
    }

    /**
     * @brief Verify rejection of a non-scalar events_path value.
     *
     * @return true when std::runtime_error is thrown.
     */
    bool verify_non_scalar_events_path() {
        return verify_runtime_error(
            "non-scalar events_path",
            test_directory() / "main_non_scalar_events_path.yaml",
            "events_path:\n"
            "  - events\n"
            "output_path: \"output\"\n"
            "number_of_events: 1\n"
            "number_of_subevents: 1\n"
            "hbt_enabled: false\n"
        );
    }

    /**
     * @brief Verify rejection of an empty events_path scalar.
     *
     * @return true when std::runtime_error is thrown.
     */
    bool verify_empty_events_path() {
        return verify_runtime_error(
            "empty events_path",
            test_directory() / "main_empty_events_path.yaml",
            "events_path: \"\"\n"
            "output_path: \"output\"\n"
            "number_of_events: 1\n"
            "number_of_subevents: 1\n"
            "hbt_enabled: false\n"
        );
    }

    /**
     * @brief Verify rejection of a duplicate events_path key.
     *
     * @return true when std::runtime_error is thrown.
     */
    bool verify_duplicate_events_path() {
        return verify_runtime_error(
            "duplicate events_path",
            test_directory() / "main_duplicate_events_path.yaml",
            "events_path: \"first\"\n"
            "output_path: \"output\"\n"
            "events_path: \"second\"\n"
            "number_of_events: 1\n"
            "number_of_subevents: 1\n"
            "hbt_enabled: false\n"
        );
    }

    /**
     * @brief Verify that output_path is required.
     *
     * @return true when std::runtime_error is thrown.
     */
    bool verify_missing_output_path() {
        return verify_runtime_error(
            "missing output_path",
            test_directory() / "main_missing_output_path.yaml",
            "events_path: \"events\"\n"
            "number_of_events: 1\n"
            "number_of_subevents: 1\n"
            "hbt_enabled: false\n"
        );
    }

    /**
     * @brief Verify rejection of a non-scalar output_path value.
     *
     * @return true when std::runtime_error is thrown.
     */
    bool verify_non_scalar_output_path() {
        return verify_runtime_error(
            "non-scalar output_path",
            test_directory() / "main_non_scalar_output_path.yaml",
            "events_path: \"events\"\n"
            "output_path:\n"
            "  - output\n"
            "number_of_events: 1\n"
            "number_of_subevents: 1\n"
            "hbt_enabled: false\n"
        );
    }

    /**
     * @brief Verify rejection of an empty output_path scalar.
     *
     * @return true when std::runtime_error is thrown.
     */
    bool verify_empty_output_path() {
        return verify_runtime_error(
            "empty output_path",
            test_directory() / "main_empty_output_path.yaml",
            "events_path: \"events\"\n"
            "output_path: \"\"\n"
            "number_of_events: 1\n"
            "number_of_subevents: 1\n"
            "hbt_enabled: false\n"
        );
    }

    /**
     * @brief Verify rejection of a duplicate output_path key.
     *
     * @return true when std::runtime_error is thrown.
     */
    bool verify_duplicate_output_path() {
        return verify_runtime_error(
            "duplicate output_path",
            test_directory() / "main_duplicate_output_path.yaml",
            "events_path: \"events\"\n"
            "output_path: \"first\"\n"
            "output_path: \"second\"\n"
            "number_of_events: 1\n"
            "number_of_subevents: 1\n"
            "hbt_enabled: false\n"
        );
    }

    /**
     * @brief Verify that number_of_events is required.
     *
     * @return true when std::runtime_error is thrown.
     */
    bool verify_missing_number_of_events() {
        return verify_runtime_error(
            "missing number_of_events",
            test_directory() / "main_missing_number_of_events.yaml",
            "events_path: \"events\"\n"
            "output_path: \"output\"\n"
            "number_of_subevents: 2000\n"
            "hbt_enabled: false\n"
        );
    }

    /**
     * @brief Verify that number_of_subevents is required.
     *
     * @return true when std::runtime_error is thrown.
     */
    bool verify_missing_number_of_subevents() {
        return verify_runtime_error(
            "missing number_of_subevents",
            test_directory() / "main_missing_number_of_subevents.yaml",
            "events_path: \"events\"\n"
            "output_path: \"output\"\n"
            "number_of_events: 500\n"
            "hbt_enabled: false\n"
        );
    }

    /**
     * @brief Verify rejection of zero number_of_events.
     *
     * @return true when std::runtime_error is thrown.
     */
    bool verify_zero_number_of_events() {
        return verify_runtime_error(
            "zero number_of_events",
            test_directory() / "main_zero_number_of_events.yaml",
            "events_path: \"events\"\n"
            "output_path: \"output\"\n"
            "number_of_events: 0\n"
            "number_of_subevents: 2000\n"
            "hbt_enabled: false\n"
        );
    }

    /**
     * @brief Verify rejection of zero number_of_subevents.
     *
     * @return true when std::runtime_error is thrown.
     */
    bool verify_zero_number_of_subevents() {
        return verify_runtime_error(
            "zero number_of_subevents",
            test_directory() / "main_zero_number_of_subevents.yaml",
            "events_path: \"events\"\n"
            "output_path: \"output\"\n"
            "number_of_events: 500\n"
            "number_of_subevents: 0\n"
            "hbt_enabled: false\n"
        );
    }

    /**
     * @brief Verify rejection of negative number_of_events.
     *
     * @return true when std::runtime_error is thrown.
     */
    bool verify_negative_number_of_events() {
        return verify_runtime_error(
            "negative number_of_events",
            test_directory() / "main_negative_number_of_events.yaml",
            "events_path: \"events\"\n"
            "output_path: \"output\"\n"
            "number_of_events: -1\n"
            "number_of_subevents: 2000\n"
            "hbt_enabled: false\n"
        );
    }

    /**
     * @brief Verify rejection of negative number_of_subevents.
     *
     * @return true when std::runtime_error is thrown.
     */
    bool verify_negative_number_of_subevents() {
        return verify_runtime_error(
            "negative number_of_subevents",
            test_directory() / "main_negative_number_of_subevents.yaml",
            "events_path: \"events\"\n"
            "output_path: \"output\"\n"
            "number_of_events: 500\n"
            "number_of_subevents: -1\n"
            "hbt_enabled: false\n"
        );
    }

    /**
     * @brief Verify rejection of a fractional event count.
     *
     * @return true when std::runtime_error is thrown.
     */
    bool verify_fractional_number_of_events() {
        return verify_runtime_error(
            "fractional number_of_events",
            test_directory() / "main_fractional_number_of_events.yaml",
            "events_path: \"events\"\n"
            "output_path: \"output\"\n"
            "number_of_events: 1.5\n"
            "number_of_subevents: 2000\n"
            "hbt_enabled: false\n"
        );
    }

    /**
     * @brief Verify rejection of textual number_of_subevents.
     *
     * @return true when std::runtime_error is thrown.
     */
    bool verify_textual_number_of_subevents() {
        return verify_runtime_error(
            "textual number_of_subevents",
            test_directory() / "main_textual_number_of_subevents.yaml",
            "events_path: \"events\"\n"
            "output_path: \"output\"\n"
            "number_of_events: 500\n"
            "number_of_subevents: many\n"
            "hbt_enabled: false\n"
        );
    }

    /**
     * @brief Verify rejection of a leading plus sign in a run-count scalar.
     *
     * @return true when std::runtime_error is thrown.
     */
    bool verify_signed_number_of_events() {
        return verify_runtime_error(
            "signed number_of_events",
            test_directory() / "main_signed_number_of_events.yaml",
            "events_path: \"events\"\n"
            "output_path: \"output\"\n"
            "number_of_events: +500\n"
            "number_of_subevents: 2000\n"
            "hbt_enabled: false\n"
        );
    }

    /**
     * @brief Verify rejection of preserved surrounding whitespace in a count.
     *
     * Quotation marks preserve the surrounding spaces in the YAML scalar so
     * that the strict textual count parser receives them.
     *
     * @return true when std::runtime_error is thrown.
     */
    bool verify_spaced_number_of_subevents() {
        return verify_runtime_error(
            "spaced number_of_subevents",
            test_directory() / "main_spaced_number_of_subevents.yaml",
            "events_path: \"events\"\n"
            "output_path: \"output\"\n"
            "number_of_events: 500\n"
            "number_of_subevents: \" 2000 \"\n"
            "hbt_enabled: false\n"
        );
    }

    /**
     * @brief Verify rejection of trailing characters in a run-count scalar.
     *
     * @return true when std::runtime_error is thrown.
     */
    bool verify_trailing_characters_number_of_events() {
        return verify_runtime_error(
            "trailing characters in number_of_events",
            test_directory() / "main_trailing_number_of_events.yaml",
            "events_path: \"events\"\n"
            "output_path: \"output\"\n"
            "number_of_events: 500events\n"
            "number_of_subevents: 2000\n"
            "hbt_enabled: false\n"
        );
    }

    /**
     * @brief Verify rejection of a run-count value outside std::size_t range.
     *
     * The decimal scalar is deliberately much larger than the maximum value of
     * std::size_t on every supported target.
     *
     * @return true when std::runtime_error is thrown.
     */
    bool verify_overflow_number_of_subevents() {
        return verify_runtime_error(
            "overflow number_of_subevents",
            test_directory() / "main_overflow_number_of_subevents.yaml",
            "events_path: \"events\"\n"
            "output_path: \"output\"\n"
            "number_of_events: 500\n"
            "number_of_subevents: "
            "999999999999999999999999999999999999999999999999999999999999\n"
            "hbt_enabled: false\n"
        );
    }

    /**
     * @brief Verify rejection of a non-scalar number_of_events value.
     *
     * @return true when std::runtime_error is thrown.
     */
    bool verify_non_scalar_number_of_events() {
        return verify_runtime_error(
            "non-scalar number_of_events",
            test_directory() / "main_non_scalar_number_of_events.yaml",
            "events_path: \"events\"\n"
            "output_path: \"output\"\n"
            "number_of_events:\n"
            "  - 500\n"
            "number_of_subevents: 2000\n"
            "hbt_enabled: false\n"
        );
    }

    /**
     * @brief Verify rejection of a non-scalar number_of_subevents value.
     *
     * @return true when std::runtime_error is thrown.
     */
    bool verify_non_scalar_number_of_subevents() {
        return verify_runtime_error(
            "non-scalar number_of_subevents",
            test_directory() / "main_non_scalar_number_of_subevents.yaml",
            "events_path: \"events\"\n"
            "output_path: \"output\"\n"
            "number_of_events: 500\n"
            "number_of_subevents:\n"
            "  - 2000\n"
            "hbt_enabled: false\n"
        );
    }

    /**
     * @brief Verify rejection of a duplicate number_of_events key.
     *
     * @return true when std::runtime_error is thrown.
     */
    bool verify_duplicate_number_of_events() {
        return verify_runtime_error(
            "duplicate number_of_events",
            test_directory() / "main_duplicate_number_of_events.yaml",
            "events_path: \"events\"\n"
            "output_path: \"output\"\n"
            "number_of_events: 500\n"
            "number_of_events: 600\n"
            "number_of_subevents: 2000\n"
            "hbt_enabled: false\n"
        );
    }

    /**
     * @brief Verify rejection of a duplicate number_of_subevents key.
     *
     * @return true when std::runtime_error is thrown.
     */
    bool verify_duplicate_number_of_subevents() {
        return verify_runtime_error(
            "duplicate number_of_subevents",
            test_directory() / "main_duplicate_number_of_subevents.yaml",
            "events_path: \"events\"\n"
            "output_path: \"output\"\n"
            "number_of_events: 500\n"
            "number_of_subevents: 2000\n"
            "number_of_subevents: 3000\n"
            "hbt_enabled: false\n"
        );
    }

    /**
     * @brief Verify rejection of a missing hbt_enabled entry.
     *
     * @return true when std::runtime_error is thrown.
     */
    bool verify_missing_hbt_enabled() {
        return verify_runtime_error(
            "missing hbt_enabled",
            test_directory() / "main_missing_hbt_enabled.yaml",
            "events_path: \"events\"\n"
            "output_path: \"output\"\n"
            "number_of_events: 1\n"
            "number_of_subevents: 1\n"
            "hbt_config_path: \"hbt.yaml\"\n"
        );
    }

    /**
     * @brief Verify rejection of an invalid boolean scalar.
     *
     * @return true when std::runtime_error is thrown.
     */
    bool verify_invalid_hbt_enabled_scalar() {
        return verify_runtime_error(
            "invalid hbt_enabled scalar",
            test_directory() / "main_invalid_boolean.yaml",
            "events_path: \"events\"\n"
            "output_path: \"output\"\n"
            "number_of_events: 1\n"
            "number_of_subevents: 1\n"
            "hbt_enabled: \"not_a_boolean\"\n"
        );
    }

    /**
     * @brief Verify rejection of a non-scalar hbt_enabled value.
     *
     * @return true when std::runtime_error is thrown.
     */
    bool verify_non_scalar_hbt_enabled() {
        return verify_runtime_error(
            "non-scalar hbt_enabled",
            test_directory() / "main_non_scalar_boolean.yaml",
            "events_path: \"events\"\n"
            "output_path: \"output\"\n"
            "number_of_events: 1\n"
            "number_of_subevents: 1\n"
            "hbt_enabled:\n"
            "  - true\n"
        );
    }

    /**
     * @brief Verify that enabled HBT requires hbt_config_path.
     *
     * @return true when std::runtime_error is thrown.
     */
    bool verify_enabled_without_path() {
        return verify_runtime_error(
            "enabled HBT without path",
            test_directory() / "main_enabled_without_path.yaml",
            "events_path: \"events\"\n"
            "output_path: \"output\"\n"
            "number_of_events: 1\n"
            "number_of_subevents: 1\n"
            "hbt_enabled: true\n"
        );
    }

    /**
     * @brief Verify rejection of an unknown global configuration key.
     *
     * @return true when std::runtime_error is thrown.
     */
    bool verify_unknown_key() {
        return verify_runtime_error(
            "unknown key",
            test_directory() / "main_unknown_key.yaml",
            "events_path: \"events\"\n"
            "output_path: \"output\"\n"
            "number_of_events: 1\n"
            "number_of_subevents: 1\n"
            "hbt_enabled: true\n"
            "hbt_config_path: \"hbt.yaml\"\n"
            "unknown_setting: 123\n"
        );
    }

    /**
     * @brief Verify rejection of a duplicate hbt_enabled key.
     *
     * @return true when std::runtime_error is thrown.
     */
    bool verify_duplicate_hbt_enabled() {
        return verify_runtime_error(
            "duplicate hbt_enabled",
            test_directory() / "main_duplicate_hbt_enabled.yaml",
            "events_path: \"events\"\n"
            "output_path: \"output\"\n"
            "number_of_events: 1\n"
            "number_of_subevents: 1\n"
            "hbt_enabled: false\n"
            "hbt_enabled: true\n"
            "hbt_config_path: \"hbt.yaml\"\n"
        );
    }

    /**
     * @brief Verify rejection of a duplicate hbt_config_path key.
     *
     * @return true when std::runtime_error is thrown.
     */
    bool verify_duplicate_hbt_config_path() {
        return verify_runtime_error(
            "duplicate hbt_config_path",
            test_directory() / "main_duplicate_hbt_path.yaml",
            "events_path: \"events\"\n"
            "output_path: \"output\"\n"
            "number_of_events: 1\n"
            "number_of_subevents: 1\n"
            "hbt_enabled: true\n"
            "hbt_config_path: \"first.yaml\"\n"
            "hbt_config_path: \"second.yaml\"\n"
        );
    }

    /**
     * @brief Verify rejection of a non-scalar YAML key.
     *
     * @return true when std::runtime_error is thrown.
     */
    bool verify_non_scalar_key() {
        return verify_runtime_error(
            "non-scalar key",
            test_directory() / "main_non_scalar_key.yaml",
            "events_path: \"events\"\n"
            "output_path: \"output\"\n"
            "number_of_events: 1\n"
            "number_of_subevents: 1\n"
            "? [hbt_enabled]\n"
            ": true\n"
        );
    }

    /**
     * @brief Verify rejection of a non-scalar hbt_config_path value.
     *
     * @return true when std::runtime_error is thrown.
     */
    bool verify_non_scalar_hbt_config_path() {
        return verify_runtime_error(
            "non-scalar hbt_config_path",
            test_directory() / "main_non_scalar_path.yaml",
            "events_path: \"events\"\n"
            "output_path: \"output\"\n"
            "number_of_events: 1\n"
            "number_of_subevents: 1\n"
            "hbt_enabled: true\n"
            "hbt_config_path:\n"
            "  - hbt.yaml\n"
        );
    }

    /**
     * @brief Verify threads defaults to one when omitted.
     * @return true when the backward-compatible serial default is resolved.
     */
    bool verify_threads_default_is_serial() {
        return verify_valid_config(
            "threads default",
            test_directory() / "main_threads_default.yaml",
            "events_path: \"events\"\n"
            "output_path: \"output\"\n"
            "number_of_events: 4\n"
            "number_of_subevents: 2\n"
            "hbt_enabled: false\n",
            test_directory() / "events",
            test_directory() / "output",
            4U,
            2U,
            false,
            std::nullopt,
            1U
        );
    }

    /**
     * @brief Verify zero and positive thread controls are accepted exactly.
     * @return true when automatic and explicit values are preserved.
     */
    bool verify_threads_values() {
        const bool automatic = verify_valid_config(
            "threads automatic",
            test_directory() / "main_threads_auto.yaml",
            "events_path: \"events\"\n"
            "output_path: \"output\"\n"
            "number_of_events: 4\n"
            "number_of_subevents: 2\n"
            "threads: 0\n"
            "hbt_enabled: false\n",
            test_directory() / "events",
            test_directory() / "output",
            4U,
            2U,
            false,
            std::nullopt,
            0U
        );
        const bool explicit_workers = verify_valid_config(
            "threads explicit",
            test_directory() / "main_threads_explicit.yaml",
            "events_path: \"events\"\n"
            "output_path: \"output\"\n"
            "number_of_events: 4\n"
            "number_of_subevents: 2\n"
            "threads: 8\n"
            "hbt_enabled: false\n",
            test_directory() / "events",
            test_directory() / "output",
            4U,
            2U,
            false,
            std::nullopt,
            8U
        );
        return automatic && explicit_workers;
    }

    /**
     * @brief Verify malformed thread controls remain configuration errors.
     * @return true when negative, non-scalar, and duplicate values are
     *         rejected.
     */
    bool verify_invalid_threads_values() {
        const bool negative = verify_runtime_error(
            "negative threads",
            test_directory() / "main_threads_negative.yaml",
            "events_path: \"events\"\n"
            "output_path: \"output\"\n"
            "number_of_events: 1\n"
            "number_of_subevents: 1\n"
            "threads: -1\n"
            "hbt_enabled: false\n"
        );
        const bool non_scalar = verify_runtime_error(
            "non-scalar threads",
            test_directory() / "main_threads_non_scalar.yaml",
            "events_path: \"events\"\n"
            "output_path: \"output\"\n"
            "number_of_events: 1\n"
            "number_of_subevents: 1\n"
            "threads: [1, 2]\n"
            "hbt_enabled: false\n"
        );
        const bool duplicate = verify_runtime_error(
            "duplicate threads",
            test_directory() / "main_threads_duplicate.yaml",
            "events_path: \"events\"\n"
            "output_path: \"output\"\n"
            "number_of_events: 1\n"
            "number_of_subevents: 1\n"
            "threads: 1\n"
            "threads: 2\n"
            "hbt_enabled: false\n"
        );
        return negative && non_scalar && duplicate;
    }

    /**
     * @brief Verify rejection of an empty hbt_config_path scalar.
     *
     * @return true when std::runtime_error is thrown.
     */
    bool verify_empty_hbt_config_path() {
        return verify_runtime_error(
            "empty hbt_config_path",
            test_directory() / "main_empty_path.yaml",
            "events_path: \"events\"\n"
            "output_path: \"output\"\n"
            "number_of_events: 1\n"
            "number_of_subevents: 1\n"
            "hbt_enabled: true\n"
            "hbt_config_path: \"\"\n"
        );
    }

}  // namespace

/**
 * @brief Run the global run-configuration loader unit tests.
 *
 * @return EXIT_SUCCESS when every test passes, otherwise EXIT_FAILURE.
 */
int main() {
    if (!prepare_test_directory()) {
        return EXIT_FAILURE;
    }

    bool success = true;

    if (!verify_relative_events_path()) {
        success = false;
    }

    if (!verify_absolute_events_path()) {
        success = false;
    }

    if (!verify_absolute_output_path()) {
        success = false;
    }

    if (!verify_enabled_relative_path()) {
        success = false;
    }

    if (!verify_enabled_absolute_path()) {
        success = false;
    }

    if (!verify_disabled_without_path()) {
        success = false;
    }

    if (!verify_disabled_with_path()) {
        success = false;
    }

    if (!verify_threads_default_is_serial()) {
        success = false;
    }

    if (!verify_threads_values()) {
        success = false;
    }

    if (!verify_invalid_threads_values()) {
        success = false;
    }

    if (!verify_missing_main_file()) {
        success = false;
    }

    if (!verify_invalid_yaml()) {
        success = false;
    }

    if (!verify_non_mapping_root()) {
        success = false;
    }

    if (!verify_missing_events_path()) {
        success = false;
    }

    if (!verify_non_scalar_events_path()) {
        success = false;
    }

    if (!verify_empty_events_path()) {
        success = false;
    }

    if (!verify_duplicate_events_path()) {
        success = false;
    }

    if (!verify_missing_output_path()) {
        success = false;
    }

    if (!verify_non_scalar_output_path()) {
        success = false;
    }

    if (!verify_empty_output_path()) {
        success = false;
    }

    if (!verify_duplicate_output_path()) {
        success = false;
    }

    if (!verify_missing_number_of_events()) {
        success = false;
    }

    if (!verify_missing_number_of_subevents()) {
        success = false;
    }

    if (!verify_zero_number_of_events()) {
        success = false;
    }

    if (!verify_zero_number_of_subevents()) {
        success = false;
    }

    if (!verify_negative_number_of_events()) {
        success = false;
    }

    if (!verify_negative_number_of_subevents()) {
        success = false;
    }

    if (!verify_fractional_number_of_events()) {
        success = false;
    }

    if (!verify_textual_number_of_subevents()) {
        success = false;
    }

    if (!verify_signed_number_of_events()) {
        success = false;
    }

    if (!verify_spaced_number_of_subevents()) {
        success = false;
    }

    if (!verify_trailing_characters_number_of_events()) {
        success = false;
    }

    if (!verify_overflow_number_of_subevents()) {
        success = false;
    }

    if (!verify_non_scalar_number_of_events()) {
        success = false;
    }

    if (!verify_non_scalar_number_of_subevents()) {
        success = false;
    }

    if (!verify_duplicate_number_of_events()) {
        success = false;
    }

    if (!verify_duplicate_number_of_subevents()) {
        success = false;
    }

    if (!verify_missing_hbt_enabled()) {
        success = false;
    }

    if (!verify_invalid_hbt_enabled_scalar()) {
        success = false;
    }

    if (!verify_non_scalar_hbt_enabled()) {
        success = false;
    }

    if (!verify_enabled_without_path()) {
        success = false;
    }

    if (!verify_unknown_key()) {
        success = false;
    }

    if (!verify_duplicate_hbt_enabled()) {
        success = false;
    }

    if (!verify_duplicate_hbt_config_path()) {
        success = false;
    }

    if (!verify_non_scalar_key()) {
        success = false;
    }

    if (!verify_non_scalar_hbt_config_path()) {
        success = false;
    }

    if (!verify_empty_hbt_config_path()) {
        success = false;
    }

    cleanup_test_directory();

    return success ? EXIT_SUCCESS : EXIT_FAILURE;
}
