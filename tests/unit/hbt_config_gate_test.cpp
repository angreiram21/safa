/**
 * @file hbt_config_gate_test.cpp
 * @brief Unit tests for conditional loading of HBT scientific configuration.
 */

#include "app/hbt_config_gate.h"

#include "hbt/channels/primitive_channel.h"

#include <cstdlib>
#include <exception>
#include <filesystem>
#include <fstream>
#include <iostream>
#include <optional>
#include <stdexcept>
#include <system_error>

namespace {

    /**
     * @brief Return the directory used by the HBT configuration-gate tests.
     *
     * @return Relative path to the temporary test directory.
     */
    std::filesystem::path test_directory() {
        return "hbt_config_gate_test_data";
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
                << "hbt_config_gate_test: failed to prepare test directory '"
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
     * @brief Write a complete valid HBT YAML configuration file.
     *
     * @param path Path of the file to create.
     *
     * @return true when the file is written successfully, otherwise false.
     */
    bool write_valid_hbt_config(
        const std::filesystem::path& path
    ) {
        std::ofstream file(path);

        if (!file) {
            std::cerr
                << "hbt_config_gate_test: failed to create test file '"
                << path.string()
                << "'.\n";
            return false;
        }

        file
            << "hbt_enabled_channels: \"pi_plus_pi_plus\"\n"
            << "\n"
            << "hbt_particle_acceptance:\n"
            << "  longitudinal_variable: \"pseudorapidity\"\n"
            << "\n"
            << "  groups:\n"
            << "    pions:\n"
            << "      longitudinal_abs_max: 0.8\n"
            << "      pt_min_gev: 0.14\n"
            << "      pt_max_gev: 4.0\n"
            << "\n"
            << "    kaons:\n"
            << "      longitudinal_abs_max: 0.8\n"
            << "      pt_min_gev: 0.4\n"
            << "      pt_max_gev: 1.4\n"
            << "\n"
            << "    nucleons:\n"
            << "      longitudinal_abs_max: 0.8\n"
            << "      pt_min_gev: 0.5\n"
            << "      pt_max_gev: 4.05\n"
            << "\n"
            << "    sigmas:\n"
            << "      longitudinal_abs_max: 0.8\n"
            << "      pt_min_gev: 1.0\n"
            << "      pt_max_gev: 10000.0\n"
            << "\n"
            << "    lambdas:\n"
            << "      longitudinal_abs_max: 0.8\n"
            << "      pt_min_gev: 0.3\n"
            << "      pt_max_gev: 10000.0\n"
            << "\n"
            << "hbt_pair_slicing:\n"
            << "  kt:\n"
            << "    enabled: false\n"
            << "  mt:\n"
            << "    enabled: false\n"
            << "\n"
            << "hbt_histograms:\n"
            << "  osl:\n"
            << "    nbins: 30\n"
            << "    min_fm: 0.0\n"
            << "    max_fm: 3.0\n"
            << "  radial:\n"
            << "    nbins: 20\n"
            << "    min_fm: 0.0\n"
            << "    max_fm: 4.0\n"
            << "  delta_t:\n"
            << "    nbins: 40\n"
            << "    min_fm_c: -2.0\n"
            << "    max_fm_c: 2.0\n"
            << "\n"
            << "hbt_origin_mode: \"all\"\n"
            << "hbt_fit_estimator: \"all\"\n";

        if (!file) {
            std::cerr
                << "hbt_config_gate_test: failed to write test file '"
                << path.string()
                << "'.\n";
            return false;
        }

        return true;
    }

    /**
     * @brief Verify that disabled HBT does not access its configuration file.
     *
     * The configured path deliberately names a file that does not exist.
     * Returning std::nullopt without an exception therefore verifies that the
     * gate does not attempt to load HBT scientific configuration when HBT is
     * disabled.
     *
     * @return true when std::nullopt is returned without file access.
     */
    bool verify_disabled_hbt_does_not_load_config() {
        const config::RunConfig run_config{
            test_directory() / "events",
            test_directory() / "output",
            1U,
            1U,
            false,
            test_directory() / "does_not_exist.yaml"
        };

        try {
            const std::optional<hbt::HBTConfig> hbt_config =
                app::load_hbt_config_if_enabled(run_config);

            if (hbt_config.has_value()) {
                std::cerr
                    << "hbt_config_gate_test: disabled HBT unexpectedly "
                    << "returned an HBTConfig.\n";
                return false;
            }

            return true;
        } catch (const std::exception& exception) {
            std::cerr
                << "hbt_config_gate_test: disabled HBT unexpectedly threw: "
                << exception.what()
                << ".\n";
            return false;
        }
    }

    /**
     * @brief Verify rejection of enabled HBT without a configuration path.
     *
     * @return true when std::invalid_argument is thrown.
     */
    bool verify_enabled_hbt_requires_config_path() {
        const config::RunConfig run_config{
            test_directory() / "events",
            test_directory() / "output",
            1U,
            1U,
            true,
            std::nullopt
        };

        try {
            static_cast<void>(
                app::load_hbt_config_if_enabled(run_config)
            );
        } catch (const std::invalid_argument&) {
            return true;
        } catch (const std::exception& exception) {
            std::cerr
                << "hbt_config_gate_test: enabled HBT without a path threw "
                << "an unexpected exception type: "
                << exception.what()
                << ".\n";
            return false;
        }

        std::cerr
            << "hbt_config_gate_test: enabled HBT without a path did not "
            << "throw std::invalid_argument.\n";
        return false;
    }

    /**
     * @brief Verify complete scientific configuration loading for enabled HBT.
     *
     * @return true when the expected resolved HBTConfig is returned.
     */
    bool verify_enabled_hbt_loads_config() {
        const std::filesystem::path hbt_config_path =
            test_directory() / "hbt.yaml";

        if (!write_valid_hbt_config(hbt_config_path)) {
            return false;
        }

        const config::RunConfig run_config{
            test_directory() / "events",
            test_directory() / "output",
            1U,
            1U,
            true,
            hbt_config_path
        };

        try {
            const std::optional<hbt::HBTConfig> hbt_config =
                app::load_hbt_config_if_enabled(run_config);

            if (!hbt_config.has_value()) {
                std::cerr
                    << "hbt_config_gate_test: enabled HBT did not return "
                    << "an HBTConfig.\n";
                return false;
            }

            if (hbt_config->selection.products.size() != 1U) {
                std::cerr
                    << "hbt_config_gate_test: loaded HBTConfig contains the "
                    << "wrong number of analysis products.\n";
                return false;
            }

            const hbt::AnalysisProduct& product =
                hbt_config->selection.products.front();

            if (product.primitive_channels.size() != 1U) {
                std::cerr
                    << "hbt_config_gate_test: loaded analysis product "
                    << "contains the wrong number of primitive channels.\n";
                return false;
            }

            if (
                product.primitive_channels.front() !=
                hbt::PrimitiveChannelId::PiPlusPiPlus
            ) {
                std::cerr
                    << "hbt_config_gate_test: loaded analysis product "
                    << "contains the wrong primitive channel.\n";
                return false;
            }

            const hbt::ParticleAcceptanceConfig& acceptance =
                hbt_config->particle_acceptance;

            if (
                acceptance.longitudinal_variable !=
                hbt::LongitudinalVariable::Pseudorapidity
            ) {
                std::cerr
                    << "hbt_config_gate_test: loaded HBTConfig contains the "
                    << "wrong longitudinal variable.\n";
                return false;
            }

            if (
                acceptance.pions.longitudinal_abs_max != 0.8 ||
                acceptance.pions.pt_min_gev != 0.14 ||
                acceptance.pions.pt_max_gev != 4.0 ||
                acceptance.kaons.longitudinal_abs_max != 0.8 ||
                acceptance.kaons.pt_min_gev != 0.4 ||
                acceptance.kaons.pt_max_gev != 1.4 ||
                acceptance.nucleons.longitudinal_abs_max != 0.8 ||
                acceptance.nucleons.pt_min_gev != 0.5 ||
                acceptance.nucleons.pt_max_gev != 4.05 ||
                acceptance.sigmas.longitudinal_abs_max != 0.8 ||
                acceptance.sigmas.pt_min_gev != 1.0 ||
                acceptance.sigmas.pt_max_gev != 10000.0 ||
                acceptance.lambdas.longitudinal_abs_max != 0.8 ||
                acceptance.lambdas.pt_min_gev != 0.3 ||
                acceptance.lambdas.pt_max_gev != 10000.0
            ) {
                std::cerr
                    << "hbt_config_gate_test: loaded HBTConfig contains "
                    << "incorrect particle-acceptance cuts.\n";
                return false;
            }

            const hbt::HBTHistogramConfig& histograms =
                hbt_config->histogram_config;
            if (histograms.osl.nbins != 30U ||
                histograms.osl.minimum != 0.0 ||
                histograms.osl.maximum != 3.0 ||
                histograms.osl.inverse_bin_width != 10.0 ||
                histograms.radial.nbins != 20U ||
                histograms.radial.minimum != 0.0 ||
                histograms.radial.maximum != 4.0 ||
                histograms.radial.inverse_bin_width != 5.0 ||
                histograms.delta_t.nbins != 40U ||
                histograms.delta_t.minimum != -2.0 ||
                histograms.delta_t.maximum != 2.0 ||
                histograms.delta_t.inverse_bin_width != 10.0) {
                std::cerr
                    << "hbt_config_gate_test: loaded HBTConfig contains "
                    << "incorrect histogram binning.\n";
                return false;
            }

            if (hbt_config->origin_mode != hbt::OriginMode::All) {
                std::cerr
                    << "hbt_config_gate_test: loaded HBTConfig contains the "
                    << "wrong origin mode.\n";
                return false;
            }

            if (hbt_config->fit_estimator_mode !=
                hbt::FitEstimatorMode::All) {
                std::cerr
                    << "hbt_config_gate_test: loaded HBTConfig has the "
                    << "wrong fit estimator mode.\n";
                return false;
            }

            return true;
        } catch (const std::exception& exception) {
            std::cerr
                << "hbt_config_gate_test: enabled HBT with a valid "
                << "configuration unexpectedly threw: "
                << exception.what()
                << ".\n";
            return false;
        }
    }

}  // namespace

/**
 * @brief Run the HBT configuration-gate unit tests.
 *
 * @return EXIT_SUCCESS when every test passes, otherwise EXIT_FAILURE.
 */
int main() {
    if (!prepare_test_directory()) {
        return EXIT_FAILURE;
    }

    bool success = true;

    if (!verify_disabled_hbt_does_not_load_config()) {
        success = false;
    }

    if (!verify_enabled_hbt_requires_config_path()) {
        success = false;
    }

    if (!verify_enabled_hbt_loads_config()) {
        success = false;
    }

    cleanup_test_directory();

    return success ? EXIT_SUCCESS : EXIT_FAILURE;
}
