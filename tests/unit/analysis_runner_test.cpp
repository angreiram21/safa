/**
 * @file analysis_runner_test.cpp
 * @brief Unit tests for ordered application startup orchestration.
 */

#include "app/analysis_runner.h"
#include "app/console_progress_bar.h"

#include "hbt/channels/primitive_channel.h"
#include "hbt/species/species.h"

#include <cstdlib>
#include <exception>
#include <filesystem>
#include <fstream>
#include <iostream>
#include <sstream>
#include <stdexcept>
#include <string>
#include <system_error>
#include <vector>

namespace {

/**
 * @brief Return the directory used by the analysis-runner tests.
 *
 * @return Relative path to the temporary test directory.
 */
std::filesystem::path test_directory() {
    return "analysis_runner_test_data";
}

/**
 * @brief Prepare an empty directory for temporary test files.
 *
 * @return true when the directory is ready, otherwise false.
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
            << "analysis_runner_test: failed to prepare test directory '"
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
 * @brief Write text to one temporary test file.
 *
 * @param path Path of the file to create.
 * @param content Complete text to write.
 *
 * @return true when the file is written successfully, otherwise false.
 */
bool write_text_file(
    const std::filesystem::path& path,
    const std::string& content
) {
    std::ofstream file(path);

    if (!file) {
        std::cerr
            << "analysis_runner_test: failed to create test file '"
            << path.string()
            << "'.\n";
        return false;
    }

    file << content;

    if (!file) {
        std::cerr
            << "analysis_runner_test: failed to write test file '"
            << path.string()
            << "'.\n";
        return false;
    }

    return true;
}

/**
 * @brief Return a complete valid HBT YAML configuration.
 *
 * @return YAML text selecting the PiPlusPiPlus primitive channel.
 */
std::string valid_hbt_config_text() {
    return
        "hbt_enabled_channels: \"pi_plus_pi_plus\"\n"
        "\n"
        "hbt_particle_acceptance:\n"
        "  longitudinal_variable: \"pseudorapidity\"\n"
        "\n"
        "  groups:\n"
        "    pions:\n"
        "      longitudinal_abs_max: 0.8\n"
        "      pt_min_gev: 0.14\n"
        "      pt_max_gev: 4.0\n"
        "\n"
        "    kaons:\n"
        "      longitudinal_abs_max: 0.8\n"
        "      pt_min_gev: 0.4\n"
        "      pt_max_gev: 1.4\n"
        "\n"
        "    nucleons:\n"
        "      longitudinal_abs_max: 0.8\n"
        "      pt_min_gev: 0.5\n"
        "      pt_max_gev: 4.05\n"
        "\n"
        "    sigmas:\n"
        "      longitudinal_abs_max: 0.8\n"
        "      pt_min_gev: 1.0\n"
        "      pt_max_gev: 10000.0\n"
        "\n"
        "    lambdas:\n"
        "      longitudinal_abs_max: 0.8\n"
        "      pt_min_gev: 0.3\n"
        "      pt_max_gev: 10000.0\n"
        "\n"
        "hbt_pair_slicing:\n"
        "  kt:\n"
        "    enabled: false\n"
        "  mt:\n"
        "    enabled: false\n"
        "\n"
        "hbt_histograms:\n"
        "  osl:\n"
        "    nbins: 30\n"
        "    min_fm: 0.0\n"
        "    max_fm: 3.0\n"
        "  radial:\n"
        "    nbins: 20\n"
        "    min_fm: 0.0\n"
        "    max_fm: 4.0\n"
        "  delta_t:\n"
        "    nbins: 40\n"
        "    min_fm_c: -2.0\n"
        "    max_fm_c: 2.0\n"
        "\n"
        "hbt_origin_mode: \"all\"\n"
        "hbt_fit_estimator: \"all\"\n";
}

/**
 * @brief Verify disabled HBT produces only global startup state.
 *
 * events_path, output_path, and hbt_config_path deliberately identify
 * nonexistent locations. Successful startup therefore also verifies that this
 * stage does not access input or output paths and does not load disabled HBT
 * configuration.
 *
 * @return true when only RunConfig is populated as expected.
 */
bool verify_disabled_hbt_builds_global_startup_only() {
    const std::filesystem::path main_path =
        test_directory() / "disabled_main.yaml";

    const std::string main_config =
        "events_path: \"missing_events\"\n"
        "output_path: \"output\"\n"
        "number_of_events: 2\n"
        "number_of_subevents: 3\n"
        "hbt_enabled: false\n"
        "hbt_config_path: \"missing_hbt.yaml\"\n";

    if (!write_text_file(main_path, main_config)) {
        return false;
    }

    try {
        const app::AnalysisStartupState state =
            app::AnalysisRunner(main_path).prepare_startup();

        const std::filesystem::path expected_events_path =
            test_directory() / "missing_events";

        if (state.run_config.events_path != expected_events_path) {
            std::cerr
                << "analysis_runner_test: disabled startup resolved the "
                << "wrong events_path.\n";
            return false;
        }

        const std::filesystem::path expected_output_path =
            test_directory() / "output";

        if (state.run_config.output_path != expected_output_path) {
            std::cerr
                << "analysis_runner_test: disabled startup resolved the "
                << "wrong output_path.\n";
            return false;
        }

        if (
            state.run_config.number_of_events != 2U ||
            state.run_config.number_of_subevents != 3U
        ) {
            std::cerr
                << "analysis_runner_test: disabled startup changed run "
                << "cardinality.\n";
            return false;
        }

        if (state.hbt_config.has_value()) {
            std::cerr
                << "analysis_runner_test: disabled HBT produced HBTConfig.\n";
            return false;
        }

        if (state.hbt_startup_state.has_value()) {
            std::cerr
                << "analysis_runner_test: disabled HBT produced "
                << "HBTStartupState.\n";
            return false;
        }

        return true;
    } catch (const std::exception& exception) {
        std::cerr
            << "analysis_runner_test: disabled startup unexpectedly threw: "
            << exception.what()
            << ".\n";
        return false;
    }
}

/**
 * @brief Verify enabled HBT is loaded before startup requirements are derived.
 *
 * @return true when HBTConfig and its derived HBTStartupState are both
 *         available and contain the expected selection requirements.
 */
bool verify_enabled_hbt_builds_complete_startup() {
    const std::filesystem::path hbt_path =
        test_directory() / "hbt.yaml";

    if (!write_text_file(hbt_path, valid_hbt_config_text())) {
        return false;
    }

    const std::filesystem::path main_path =
        test_directory() / "enabled_main.yaml";

    const std::string main_config =
        "events_path: \"missing_events\"\n"
        "output_path: \"output\"\n"
        "number_of_events: 4\n"
        "number_of_subevents: 5\n"
        "hbt_enabled: true\n"
        "hbt_config_path: \"hbt.yaml\"\n";

    if (!write_text_file(main_path, main_config)) {
        return false;
    }

    try {
        const app::AnalysisStartupState state =
            app::AnalysisRunner(main_path).prepare_startup();

        if (!state.hbt_config.has_value()) {
            std::cerr
                << "analysis_runner_test: enabled HBT did not produce "
                << "HBTConfig.\n";
            return false;
        }

        if (!state.hbt_startup_state.has_value()) {
            std::cerr
                << "analysis_runner_test: enabled HBT did not produce "
                << "HBTStartupState.\n";
            return false;
        }

        const hbt::HBTStartupState& hbt_startup =
            state.hbt_startup_state.value();

        const std::vector<hbt::PrimitiveChannelId> expected_channels{
            hbt::PrimitiveChannelId::PiPlusPiPlus
        };

        if (
            hbt_startup.required_primitive_channels !=
            expected_channels
        ) {
            std::cerr
                << "analysis_runner_test: enabled startup derived the "
                << "wrong primitive-channel requirements.\n";
            return false;
        }

        const std::vector<hbt::SpeciesId> expected_species{
            hbt::SpeciesId::PiPlus
        };

        if (hbt_startup.required_species != expected_species) {
            std::cerr
                << "analysis_runner_test: enabled startup derived the "
                << "wrong species requirements.\n";
            return false;
        }

        return true;
    } catch (const std::exception& exception) {
        std::cerr
            << "analysis_runner_test: enabled startup unexpectedly threw: "
            << exception.what()
            << ".\n";
        return false;
    }
}

/**
 * @brief Verify global configuration failure occurs before HBT file loading.
 *
 * The global configuration deliberately omits events_path while naming an HBT
 * file that also does not exist. The reported failure must therefore be the
 * missing global entry, proving that RunConfig resolution is the first startup
 * stage.
 *
 * @return true when the failure identifies the missing events_path entry.
 */
bool verify_global_config_failure_precedes_hbt_loading() {
    const std::filesystem::path main_path =
        test_directory() / "invalid_main.yaml";

    const std::string main_config =
        "output_path: \"output\"\n"
        "number_of_events: 1\n"
        "number_of_subevents: 1\n"
        "hbt_enabled: true\n"
        "hbt_config_path: \"does_not_exist.yaml\"\n";

    if (!write_text_file(main_path, main_config)) {
        return false;
    }

    try {
        static_cast<void>(
            app::AnalysisRunner(main_path).prepare_startup()
        );
    } catch (const std::runtime_error& exception) {
        const std::string message = exception.what();

        if (message.find("'events_path'") != std::string::npos) {
            return true;
        }

        std::cerr
            << "analysis_runner_test: invalid global configuration failed "
            << "for the wrong reason: "
            << message
            << ".\n";
        return false;
    } catch (const std::exception& exception) {
        std::cerr
            << "analysis_runner_test: invalid global configuration threw "
            << "an unexpected exception type: "
            << exception.what()
            << ".\n";
        return false;
    }

    std::cerr
        << "analysis_runner_test: invalid global configuration did not "
        << "throw.\n";
    return false;
}


/**
 * @brief Remove ANSI CSI sequences from one captured terminal rendering.
 * @param text Captured terminal text.
 * @return Text with terminal control sequences removed.
 */
std::string strip_ansi_sequences(const std::string& text) {
    std::string plain;
    plain.reserve(text.size());

    for (std::size_t index = 0U; index < text.size();) {
        if (
            text[index] == '\033' &&
            index + 1U < text.size() &&
            text[index + 1U] == '['
        ) {
            index += 2U;

            while (index < text.size()) {
                const unsigned char character =
                    static_cast<unsigned char>(text[index]);
                ++index;

                if (character >= 0x40U && character <= 0x7eU) {
                    break;
                }
            }

            continue;
        }

        plain.push_back(text[index]);
        ++index;
    }

    return plain;
}

/**
 * @brief Verify interactive progress is green and carries global context.
 * @return true when the expected colored three-line block is present.
 */
bool verify_console_progress_interactive_rendering() {
    std::ostringstream output;
    app::ConsoleProgressBar progress(output, true, 40U);

    progress.begin(2U, 10U);
    const std::string rendered = output.str();

    if (rendered.find("\033[92m") == std::string::npos) {
        std::cerr
            << "analysis_runner_test: progress is not bright green.\n";
        return false;
    }

    const std::string plain = strip_ansi_sequences(rendered);

    if (
        plain.find("HBT [") == std::string::npos ||
        plain.find("  0.0%") == std::string::npos ||
        plain.find("0 / 10 subevents") == std::string::npos ||
        plain.find("event 1 / 2") == std::string::npos
    ) {
        std::cerr
            << "analysis_runner_test: interactive progress content is "
            << "incomplete.\n";
        return false;
    }

    return true;
}

/**
 * @brief Verify adaptive progress never exceeds a narrow terminal width.
 * @return true when every visible line fits the configured width.
 */
bool verify_console_progress_adapts_to_terminal_width() {
    constexpr std::size_t width = 20U;
    std::ostringstream output;
    app::ConsoleProgressBar progress(output, true, width);

    progress.begin(12U, 12000U);
    const std::string plain = strip_ansi_sequences(output.str());
    std::istringstream lines(plain);
    std::string line;

    while (std::getline(lines, line)) {
        if (line.size() > width) {
            std::cerr
                << "analysis_runner_test: progress line exceeds terminal "
                << "width.\n";
            return false;
        }
    }

    return true;
}

/**
 * @brief Verify redirected progress is plain and throttled to milestones.
 * @return true when log output has no ANSI escapes and includes final stages.
 */
bool verify_console_progress_redirected_rendering() {
    std::ostringstream output;
    app::ConsoleProgressBar progress(output, false, 80U);

    progress.begin(1U, 20U);

    for (std::size_t completed = 1U; completed <= 20U; ++completed) {
        progress.subevent_completed(completed, 1U, completed);
    }

    progress.begin_postprocessing();
    progress.analysis_complete();
    progress.begin_output();
    progress.finish();

    const std::string rendered = output.str();

    if (rendered.find("\033[") != std::string::npos) {
        std::cerr
            << "analysis_runner_test: redirected progress contains ANSI "
            << "escapes.\n";
        return false;
    }

    if (
        rendered.find("HBT progress: 100% (20/20 subevents)") ==
            std::string::npos ||
        rendered.find("HBT stage: post-sample statistical analysis") ==
            std::string::npos ||
        rendered.find("HBT stage: writing production output") ==
            std::string::npos ||
        rendered.find("HBT stage: complete") == std::string::npos
    ) {
        std::cerr
            << "analysis_runner_test: redirected progress is incomplete.\n";
        return false;
    }

    return true;
}

}  // namespace

/**
 * @brief Run the analysis-runner unit tests.
 *
 * @return EXIT_SUCCESS when every test passes, otherwise EXIT_FAILURE.
 */
int main() {
    if (!prepare_test_directory()) {
        return EXIT_FAILURE;
    }

    bool success = true;

    if (!verify_disabled_hbt_builds_global_startup_only()) {
        success = false;
    }

    if (!verify_enabled_hbt_builds_complete_startup()) {
        success = false;
    }

    if (!verify_global_config_failure_precedes_hbt_loading()) {
        success = false;
    }

    if (!verify_console_progress_interactive_rendering()) {
        success = false;
    }

    if (!verify_console_progress_adapts_to_terminal_width()) {
        success = false;
    }

    if (!verify_console_progress_redirected_rendering()) {
        success = false;
    }

    cleanup_test_directory();

    return success ? EXIT_SUCCESS : EXIT_FAILURE;
}
