/**
 * @file hbt_config_loader_test.cpp
 * @brief Unit tests for loading the scientific HBT configuration from YAML.
 *
 * This test verifies conversion of the complete HBT YAML configuration
 * contract into HBTConfig. It also verifies strict YAML structure and
 * scientific-value validation.
 */

#include "hbt/config/hbt_config_loader.h"

#include <cstdlib>
#include <exception>
#include <filesystem>
#include <fstream>
#include <iostream>
#include <stdexcept>
#include <string>
#include <string_view>
#include <system_error>
#include <vector>

namespace {

    /**
     * @brief Description of one YAML replacement-based test case.
     *
     * original identifies one exact substring in the valid base document.
     * replacement provides the modified text used by the test.
     */
    struct ReplacementCase {
        /// Stable description of the test case.
        const char* test_name;
        /// Exact substring replaced in the valid base document.
        const char* original;
        /// Replacement text injected by the test.
        const char* replacement;
        /// Temporary YAML filename used by the test.
        const char* filename;
    };

    /**
     * @brief Write textual YAML content to one temporary test file.
     *
     * @param path Path of the file to create.
     * @param contents Complete YAML document to write.
     *
     * @return true when the file is written successfully, otherwise false.
     */
    bool write_test_file(
        const std::filesystem::path& path,
        std::string_view contents
    ) {
        std::ofstream file(path);

        if (!file) {
            std::cerr
                << "hbt_config_loader_test: failed to create test file '"
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
                << "hbt_config_loader_test: failed to write test file '"
                << path.string()
                << "'.\n";
            return false;
        }

        return true;
    }

    /**
     * @brief Remove one temporary test file without throwing.
     *
     * @param path Path of the file to remove.
     */
    void remove_test_file(
        const std::filesystem::path& path
    ) noexcept {
        std::error_code error;

        static_cast<void>(
            std::filesystem::remove(path, error)
        );
    }

    /**
     * @brief Construct one complete valid HBT YAML document.
     *
     * @param longitudinal_variable Canonical longitudinal-variable token.
     * @param origin_mode Canonical origin-selection token.
     *
     * @return Complete valid HBT YAML document.
     */
    std::string make_valid_config(
        std::string_view longitudinal_variable = "pseudorapidity",
        std::string_view origin_mode = "all"
    ) {
        std::string text;

        text +=
            "hbt_enabled_channels: "
            "\"pi_plus_pi_plus,k_plus_p+k_minus_p_bar\"\n"
            "\n"
            "hbt_particle_acceptance:\n"
            "  longitudinal_variable: \"";
        text += longitudinal_variable;
        text +=
            "\"\n"
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
            "hbt_origin_mode: \"";
        text += origin_mode;
        text += "\"\n";

        return text;
    }

    /**
     * @brief Replace one exact substring in a YAML test document.
     *
     * @param text YAML document to modify.
     * @param original Exact substring that must be present.
     * @param replacement Replacement text.
     *
     * @return true when the substring is found and replaced, otherwise false.
     */
    bool replace_once(
        std::string& text,
        std::string_view original,
        std::string_view replacement
    ) {
        const std::size_t position =
            text.find(original);

        if (position == std::string::npos) {
            return false;
        }

        text.replace(
            position,
            original.size(),
            replacement
        );

        return true;
    }

    /**
     * @brief Verify the complete scientific state of one valid configuration.
     *
     * @param test_name Human-readable test-case name.
     * @param path Path used for the temporary YAML file.
     * @param contents Complete YAML document.
     * @param expected_variable Expected longitudinal-variable identifier.
     * @param expected_origin_mode Expected origin-selection mode.
     *
     * @return true when every resolved HBTConfig member matches expectation.
     */
    bool verify_valid_config(
        const char* test_name,
        const std::filesystem::path& path,
        std::string_view contents,
        hbt::LongitudinalVariable expected_variable,
        hbt::OriginMode expected_origin_mode
    ) {
        remove_test_file(path);

        if (!write_test_file(path, contents)) {
            return false;
        }

        try {
            const hbt::HBTConfig config =
                hbt::load_hbt_config(path);

            remove_test_file(path);

            if (config.selection.products.size() != 2U) {
                std::cerr
                    << "hbt_config_loader_test: "
                    << test_name
                    << " produced the wrong product count.\n";
                return false;
            }

            const auto& first_channels =
                config.selection.products[0].primitive_channels;

            if (
                first_channels.size() != 1U ||
                first_channels[0] !=
                    hbt::PrimitiveChannelId::PiPlusPiPlus
            ) {
                std::cerr
                    << "hbt_config_loader_test: "
                    << test_name
                    << " produced an incorrect first product.\n";
                return false;
            }

            const auto& second_channels =
                config.selection.products[1].primitive_channels;

            if (
                second_channels.size() != 2U ||
                second_channels[0] !=
                    hbt::PrimitiveChannelId::KPlusProton ||
                second_channels[1] !=
                    hbt::PrimitiveChannelId::KMinusProtonBar
            ) {
                std::cerr
                    << "hbt_config_loader_test: "
                    << test_name
                    << " produced an incorrect second product.\n";
                return false;
            }

            const hbt::ParticleAcceptanceConfig& acceptance =
                config.particle_acceptance;

            if (
                acceptance.longitudinal_variable !=
                expected_variable
            ) {
                std::cerr
                    << "hbt_config_loader_test: "
                    << test_name
                    << " resolved the wrong longitudinal variable.\n";
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
                    << "hbt_config_loader_test: "
                    << test_name
                    << " produced incorrect particle-acceptance cuts.\n";
                return false;
            }

            const hbt::PairSlicingConfig& slicing =
                config.pair_slicing;

            if (
                slicing.kt.enabled ||
                !slicing.kt.bin_edges_gev.empty() ||
                slicing.mt.enabled ||
                !slicing.mt.bin_edges_gev.empty()
            ) {
                std::cerr
                    << "hbt_config_loader_test: "
                    << test_name
                    << " produced incorrect disabled pair slicing.\n";
                return false;
            }

            const hbt::HBTHistogramConfig& histograms =
                config.histogram_config;
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
                    << "hbt_config_loader_test: "
                    << test_name
                    << " produced incorrect histogram binning.\n";
                return false;
            }

            if (config.origin_mode != expected_origin_mode) {
                std::cerr
                    << "hbt_config_loader_test: "
                    << test_name
                    << " resolved the wrong origin mode.\n";
                return false;
            }

            return true;
        } catch (const std::exception& exception) {
            remove_test_file(path);

            std::cerr
                << "hbt_config_loader_test: "
                << test_name
                << " unexpectedly threw: "
                << exception.what()
                << ".\n";
            return false;
        }
    }

    /**
     * @brief Verify one complete pair-slicing configuration.
     *
     * @param test_name Human-readable test-case name.
     * @param filename Temporary YAML filename.
     * @param contents Complete YAML document.
     * @param expected_kt_enabled Expected kT slicing activation.
     * @param expected_kt_edges Expected kT bin edges in GeV.
     * @param expected_mt_enabled Expected mT slicing activation.
     * @param expected_mt_edges Expected mT bin edges in GeV.
     *
     * @return true when both resolved slicing axes match expectation.
     */
    bool verify_pair_slicing_config(
        const char* test_name,
        const char* filename,
        std::string_view contents,
        bool expected_kt_enabled,
        const std::vector<double>& expected_kt_edges,
        bool expected_mt_enabled,
        const std::vector<double>& expected_mt_edges
    ) {
        const std::filesystem::path path = filename;
        remove_test_file(path);

        if (!write_test_file(path, contents)) {
            return false;
        }

        try {
            const hbt::HBTConfig config =
                hbt::load_hbt_config(path);

            remove_test_file(path);

            if (
                config.pair_slicing.kt.enabled != expected_kt_enabled ||
                config.pair_slicing.kt.bin_edges_gev != expected_kt_edges ||
                config.pair_slicing.mt.enabled != expected_mt_enabled ||
                config.pair_slicing.mt.bin_edges_gev != expected_mt_edges
            ) {
                std::cerr
                    << "hbt_config_loader_test: "
                    << test_name
                    << " produced incorrect pair slicing.\n";
                return false;
            }

            return true;
        } catch (const std::exception& exception) {
            remove_test_file(path);

            std::cerr
                << "hbt_config_loader_test: "
                << test_name
                << " unexpectedly threw: "
                << exception.what()
                << ".\n";
            return false;
        }
    }

    /**
     * @brief Verify that one YAML document throws std::runtime_error.
     *
     * @param test_name Human-readable test-case name.
     * @param path Path used for the temporary YAML file.
     * @param contents Complete YAML document.
     *
     * @return true when std::runtime_error is thrown, otherwise false.
     */
    bool verify_runtime_error(
        const char* test_name,
        const std::filesystem::path& path,
        std::string_view contents
    ) {
        remove_test_file(path);

        if (!write_test_file(path, contents)) {
            return false;
        }

        try {
            static_cast<void>(
                hbt::load_hbt_config(path)
            );
        } catch (const std::runtime_error&) {
            remove_test_file(path);
            return true;
        } catch (const std::exception& exception) {
            remove_test_file(path);

            std::cerr
                << "hbt_config_loader_test: "
                << test_name
                << " threw an unexpected exception type: "
                << exception.what()
                << ".\n";
            return false;
        }

        remove_test_file(path);

        std::cerr
            << "hbt_config_loader_test: "
            << test_name
            << " did not throw std::runtime_error.\n";
        return false;
    }

    /**
     * @brief Verify that one YAML document throws std::invalid_argument.
     *
     * @param test_name Human-readable test-case name.
     * @param path Path used for the temporary YAML file.
     * @param contents Complete YAML document.
     *
     * @return true when std::invalid_argument is thrown, otherwise false.
     */
    bool verify_invalid_argument(
        const char* test_name,
        const std::filesystem::path& path,
        std::string_view contents
    ) {
        remove_test_file(path);

        if (!write_test_file(path, contents)) {
            return false;
        }

        try {
            static_cast<void>(
                hbt::load_hbt_config(path)
            );
        } catch (const std::invalid_argument&) {
            remove_test_file(path);
            return true;
        } catch (const std::exception& exception) {
            remove_test_file(path);

            std::cerr
                << "hbt_config_loader_test: "
                << test_name
                << " threw an unexpected exception type: "
                << exception.what()
                << ".\n";
            return false;
        }

        remove_test_file(path);

        std::cerr
            << "hbt_config_loader_test: "
            << test_name
            << " did not throw std::invalid_argument.\n";
        return false;
    }

    /**
     * @brief Verify all four independent kT/mT slicing activation states.
     *
     * @return true when disabled, kT-only, mT-only, and combined slicing all
     *         resolve exactly as configured.
     */
    bool verify_pair_slicing_combinations() {
        const std::string original =
            "hbt_pair_slicing:\n"
            "  kt:\n"
            "    enabled: false\n"
            "  mt:\n"
            "    enabled: false\n";

        {
            const std::string contents = make_valid_config();

            if (!verify_pair_slicing_config(
                    "both slicing axes disabled",
                    "hbt_config_loader_test_slicing_disabled.yaml",
                    contents,
                    false,
                    {},
                    false,
                    {}
                )) {
                return false;
            }
        }

        {
            std::string contents = make_valid_config();
            const std::string replacement =
                "hbt_pair_slicing:\n"
                "  kt:\n"
                "    enabled: false\n"
                "  mt:\n"
                "    enabled: false\n"
                "    bin_edges_gev: [0.5, 0.75, 1.0]\n";

            if (!replace_once(contents, original, replacement)) {
                return false;
            }

            if (!verify_pair_slicing_config(
                    "disabled mT slicing with retained binning",
                    "hbt_config_loader_test_slicing_retained.yaml",
                    contents,
                    false,
                    {},
                    false,
                    {0.5, 0.75, 1.0}
                )) {
                return false;
            }
        }

        {
            std::string contents = make_valid_config();
            const std::string replacement =
                "hbt_pair_slicing:\n"
                "  kt:\n"
                "    enabled: true\n"
                "    bin_edges_gev: [0.2, 0.4, 0.6]\n"
                "  mt:\n"
                "    enabled: false\n";

            if (!replace_once(contents, original, replacement)) {
                return false;
            }

            if (!verify_pair_slicing_config(
                    "kT-only slicing",
                    "hbt_config_loader_test_slicing_kt.yaml",
                    contents,
                    true,
                    {0.2, 0.4, 0.6},
                    false,
                    {}
                )) {
                return false;
            }
        }

        {
            std::string contents = make_valid_config();
            const std::string replacement =
                "hbt_pair_slicing:\n"
                "  kt:\n"
                "    enabled: false\n"
                "  mt:\n"
                "    enabled: true\n"
                "    bin_edges_gev: [0.5, 0.75, 1.0]\n";

            if (!replace_once(contents, original, replacement)) {
                return false;
            }

            if (!verify_pair_slicing_config(
                    "mT-only slicing",
                    "hbt_config_loader_test_slicing_mt.yaml",
                    contents,
                    false,
                    {},
                    true,
                    {0.5, 0.75, 1.0}
                )) {
                return false;
            }
        }

        {
            std::string contents = make_valid_config();
            const std::string replacement =
                "hbt_pair_slicing:\n"
                "  kt:\n"
                "    enabled: true\n"
                "    bin_edges_gev: [0.2, 0.4]\n"
                "  mt:\n"
                "    enabled: true\n"
                "    bin_edges_gev: [0.5, 0.75, 1.0, 1.25]\n";

            if (!replace_once(contents, original, replacement)) {
                return false;
            }

            if (!verify_pair_slicing_config(
                    "combined kT and mT slicing",
                    "hbt_config_loader_test_slicing_both.yaml",
                    contents,
                    true,
                    {0.2, 0.4},
                    true,
                    {0.5, 0.75, 1.0, 1.25}
                )) {
                return false;
            }
        }

        return true;
    }

    /**
     * @brief Verify all canonical longitudinal and origin tokens.
     *
     * @return true when every accepted token resolves correctly.
     */
    bool verify_valid_tokens() {
        /**
         * @brief One accepted longitudinal/origin token combination.
         */
        struct ValidCase {
            /// Configured longitudinal-variable token.
            const char* longitudinal_variable;
            /// Configured origin-mode token.
            const char* origin_mode;
            /// Expected parsed longitudinal-variable value.
            hbt::LongitudinalVariable expected_variable;
            /// Expected parsed origin-mode value.
            hbt::OriginMode expected_origin_mode;
            /// Temporary YAML filename used by the test.
            const char* filename;
        };

        const ValidCase cases[] = {
            {
                "rapidity",
                "primordial",
                hbt::LongitudinalVariable::Rapidity,
                hbt::OriginMode::Primordial,
                "hbt_config_loader_test_rapidity_primordial.yaml"
            },
            {
                "pseudorapidity",
                "primordial_rescattering",
                hbt::LongitudinalVariable::Pseudorapidity,
                hbt::OriginMode::PrimordialRescattering,
                "hbt_config_loader_test_origin_rescattering.yaml"
            },
            {
                "pseudorapidity",
                "primordial_rescattering_decay",
                hbt::LongitudinalVariable::Pseudorapidity,
                hbt::OriginMode::PrimordialRescatteringDecay,
                "hbt_config_loader_test_origin_decay.yaml"
            },
            {
                "pseudorapidity",
                "all",
                hbt::LongitudinalVariable::Pseudorapidity,
                hbt::OriginMode::All,
                "hbt_config_loader_test_origin_all.yaml"
            }
        };

        for (const ValidCase& test_case : cases) {
            if (
                !verify_valid_config(
                    test_case.origin_mode,
                    test_case.filename,
                    make_valid_config(
                        test_case.longitudinal_variable,
                        test_case.origin_mode
                    ),
                    test_case.expected_variable,
                    test_case.expected_origin_mode
                )
            ) {
                return false;
            }
        }

        return true;
    }

    /**
     * @brief Verify that unquoted enabled-channel scalars remain valid YAML.
     *
     * @return true when the unquoted scalar resolves successfully.
     */
    bool verify_unquoted_enabled_channels() {
        std::string contents =
            make_valid_config();

        if (
            !replace_once(
                contents,
                "\"pi_plus_pi_plus,k_plus_p+k_minus_p_bar\"",
                "pi_plus_pi_plus,k_plus_p+k_minus_p_bar"
            )
        ) {
            return false;
        }

        return verify_valid_config(
            "unquoted enabled channels",
            "hbt_config_loader_test_unquoted_channels.yaml",
            contents,
            hbt::LongitudinalVariable::Pseudorapidity,
            hbt::OriginMode::All
        );
    }

    /**
     * @brief Verify folded enabled-channel syntax used in production.
     *
     * YAML >- folds physical line breaks into spaces and strips the final line
     * break. The enabled-channel parser already accepts surrounding ASCII
     * spaces around products, so the resolved selection must remain unchanged.
     *
     * @return true when the folded block scalar resolves successfully.
     */
    bool verify_folded_enabled_channels() {
        std::string contents =
            make_valid_config();

        if (
            !replace_once(
                contents,
                "hbt_enabled_channels: "
                "\"pi_plus_pi_plus,k_plus_p+k_minus_p_bar\"\n",
                "hbt_enabled_channels: >-\n"
                "  pi_plus_pi_plus,\n"
                "  k_plus_p+k_minus_p_bar\n"
            )
        ) {
            return false;
        }

        return verify_valid_config(
            "folded enabled channels",
            "hbt_config_loader_test_folded_channels.yaml",
            contents,
            hbt::LongitudinalVariable::Pseudorapidity,
            hbt::OriginMode::All
        );
    }

    /**
     * @brief Verify direct file and YAML structure failures.
     *
     * @return true when every direct structural failure is rejected.
     */
    bool verify_direct_runtime_errors() {
        const std::filesystem::path missing_path =
            "hbt_config_loader_test_missing.yaml";

        remove_test_file(missing_path);

        try {
            static_cast<void>(
                hbt::load_hbt_config(missing_path)
            );

            std::cerr
                << "hbt_config_loader_test: missing file did not throw.\n";
            return false;
        } catch (const std::runtime_error&) {
        } catch (const std::exception& exception) {
            std::cerr
                << "hbt_config_loader_test: missing file threw: "
                << exception.what()
                << ".\n";
            return false;
        }

        /**
         * @brief One direct YAML structure failure case.
         */
        struct DirectCase {
            /// Stable description of the direct failure case.
            const char* test_name;
            /// Complete YAML document supplied to the loader.
            const char* contents;
            /// Temporary YAML filename used by the test.
            const char* filename;
        };

        const DirectCase cases[] = {
            {
                "invalid YAML",
                "hbt_enabled_channels: [\n",
                "hbt_config_loader_test_invalid_yaml.yaml"
            },
            {
                "non-mapping root",
                "- pi_plus_pi_plus\n",
                "hbt_config_loader_test_non_mapping_root.yaml"
            },
            {
                "non-scalar root key",
                "? [one, two]\n"
                ": value\n",
                "hbt_config_loader_test_non_scalar_root_key.yaml"
            },
            {
                "missing hbt_particle_acceptance",
                "hbt_enabled_channels: \"pi_plus_pi_plus\"\n"
                "hbt_origin_mode: \"all\"\n",
                "hbt_config_loader_test_missing_acceptance.yaml"
            },
            {
                "non-mapping particle acceptance",
                "hbt_enabled_channels: \"pi_plus_pi_plus\"\n"
                "hbt_particle_acceptance: 1\n"
                "hbt_origin_mode: \"all\"\n",
                "hbt_config_loader_test_non_mapping_acceptance.yaml"
            },
            {
                "missing groups",
                "hbt_enabled_channels: \"pi_plus_pi_plus\"\n"
                "hbt_particle_acceptance:\n"
                "  longitudinal_variable: \"pseudorapidity\"\n"
                "hbt_origin_mode: \"all\"\n",
                "hbt_config_loader_test_groups_missing.yaml"
            },
            {
                "non-mapping groups",
                "hbt_enabled_channels: \"pi_plus_pi_plus\"\n"
                "hbt_particle_acceptance:\n"
                "  longitudinal_variable: \"pseudorapidity\"\n"
                "  groups: 1\n"
                "hbt_origin_mode: \"all\"\n",
                "hbt_config_loader_test_groups_non_mapping.yaml"
            }
        };

        for (const DirectCase& test_case : cases) {
            if (
                !verify_runtime_error(
                    test_case.test_name,
                    test_case.filename,
                    test_case.contents
                )
            ) {
                return false;
            }
        }

        return true;
    }

    /**
     * @brief Verify strict mapping and scalar structure by YAML replacement.
     *
     * @return true when every malformed structure is rejected.
     */
    bool verify_replacement_runtime_errors() {
        const ReplacementCase cases[] = {
            {
                "missing hbt_enabled_channels",
                "hbt_enabled_channels: "
                "\"pi_plus_pi_plus,k_plus_p+k_minus_p_bar\"\n"
                "\n",
                "",
                "hbt_config_loader_test_missing_channels.yaml"
            },
            {
                "missing hbt_origin_mode",
                "hbt_origin_mode: \"all\"\n",
                "",
                "hbt_config_loader_test_missing_origin.yaml"
            },
            {
                "unknown root key",
                "hbt_origin_mode: \"all\"\n",
                "hbt_origin_mode: \"all\"\n"
                "unknown_setting: 1\n",
                "hbt_config_loader_test_unknown_root.yaml"
            },
            {
                "duplicate root key",
                "hbt_origin_mode: \"all\"\n",
                "hbt_origin_mode: \"all\"\n"
                "hbt_origin_mode: \"all\"\n",
                "hbt_config_loader_test_duplicate_root.yaml"
            },
            {
                "missing hbt_histograms",
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
                "\n",
                "",
                "hbt_config_loader_test_missing_histograms.yaml"
            },
            {
                "unknown histogram family",
                "  delta_t:\n"
                "    nbins: 40\n"
                "    min_fm_c: -2.0\n"
                "    max_fm_c: 2.0\n",
                "  delta_t:\n"
                "    nbins: 40\n"
                "    min_fm_c: -2.0\n"
                "    max_fm_c: 2.0\n"
                "  extra: {}\n",
                "hbt_config_loader_test_unknown_histogram_family.yaml"
            },
            {
                "missing histogram nbins",
                "  osl:\n"
                "    nbins: 30\n"
                "    min_fm: 0.0\n"
                "    max_fm: 3.0\n",
                "  osl:\n"
                "    min_fm: 0.0\n"
                "    max_fm: 3.0\n",
                "hbt_config_loader_test_missing_histogram_nbins.yaml"
            },
            {
                "non-scalar hbt_enabled_channels",
                "hbt_enabled_channels: "
                "\"pi_plus_pi_plus,k_plus_p+k_minus_p_bar\"\n",
                "hbt_enabled_channels:\n"
                "  - pi_plus_pi_plus\n",
                "hbt_config_loader_test_non_scalar_channels.yaml"
            },
            {
                "unknown particle-acceptance key",
                "  groups:\n",
                "  unknown: 1\n"
                "  groups:\n",
                "hbt_config_loader_test_acceptance_unknown.yaml"
            },
            {
                "duplicate longitudinal variable",
                "  longitudinal_variable: \"pseudorapidity\"\n",
                "  longitudinal_variable: \"pseudorapidity\"\n"
                "  longitudinal_variable: \"pseudorapidity\"\n",
                "hbt_config_loader_test_acceptance_duplicate.yaml"
            },
            {
                "missing longitudinal variable",
                "  longitudinal_variable: \"pseudorapidity\"\n"
                "\n",
                "",
                "hbt_config_loader_test_acceptance_missing.yaml"
            },
            {
                "unknown particle group",
                "    pions:\n",
                "    unknown: {}\n"
                "    pions:\n",
                "hbt_config_loader_test_group_unknown.yaml"
            },
            {
                "duplicate particle group",
                "    pions:\n",
                "    pions: {}\n"
                "    pions:\n",
                "hbt_config_loader_test_group_duplicate.yaml"
            },
            {
                "missing particle group",
                "    lambdas:\n"
                "      longitudinal_abs_max: 0.8\n"
                "      pt_min_gev: 0.3\n"
                "      pt_max_gev: 10000.0\n"
                "\n",
                "",
                "hbt_config_loader_test_group_missing.yaml"
            },
            {
                "non-mapping particle group",
                "    pions:\n"
                "      longitudinal_abs_max: 0.8\n"
                "      pt_min_gev: 0.14\n"
                "      pt_max_gev: 4.0\n",
                "    pions: 1\n",
                "hbt_config_loader_test_group_non_mapping.yaml"
            },
            {
                "unknown group cut",
                "      pt_min_gev: 0.14\n",
                "      unknown: 1\n"
                "      pt_min_gev: 0.14\n",
                "hbt_config_loader_test_cut_unknown.yaml"
            },
            {
                "duplicate group cut",
                "      pt_min_gev: 0.14\n",
                "      pt_min_gev: 0.14\n"
                "      pt_min_gev: 0.14\n",
                "hbt_config_loader_test_cut_duplicate.yaml"
            },
            {
                "missing group cut",
                "      pt_max_gev: 4.0\n",
                "",
                "hbt_config_loader_test_cut_missing.yaml"
            },
            {
                "non-scalar group cut",
                "      pt_min_gev: 0.14\n",
                "      pt_min_gev:\n"
                "        - 0.14\n",
                "hbt_config_loader_test_cut_non_scalar.yaml"
            },
            {
                "non-scalar longitudinal variable",
                "  longitudinal_variable: \"pseudorapidity\"\n",
                "  longitudinal_variable:\n"
                "    - pseudorapidity\n",
                "hbt_config_loader_test_longitudinal_non_scalar.yaml"
            },
            {
                "missing pair slicing",
                "hbt_pair_slicing:\n"
                "  kt:\n"
                "    enabled: false\n"
                "  mt:\n"
                "    enabled: false\n"
                "\n",
                "",
                "hbt_config_loader_test_slicing_missing.yaml"
            },
            {
                "non-mapping pair slicing",
                "hbt_pair_slicing:\n"
                "  kt:\n"
                "    enabled: false\n"
                "  mt:\n"
                "    enabled: false\n",
                "hbt_pair_slicing: 1\n",
                "hbt_config_loader_test_slicing_non_mapping.yaml"
            },
            {
                "unknown pair-slicing axis",
                "hbt_pair_slicing:\n",
                "hbt_pair_slicing:\n"
                "  unknown: {}\n",
                "hbt_config_loader_test_slicing_unknown_axis.yaml"
            },
            {
                "missing mT slicing axis",
                "  mt:\n"
                "    enabled: false\n",
                "",
                "hbt_config_loader_test_slicing_missing_mt.yaml"
            },
            {
                "non-mapping kT slicing axis",
                "  kt:\n"
                "    enabled: false\n",
                "  kt: 1\n",
                "hbt_config_loader_test_slicing_kt_non_mapping.yaml"
            },
            {
                "missing slicing enabled flag",
                "  kt:\n"
                "    enabled: false\n",
                "  kt: {}\n",
                "hbt_config_loader_test_slicing_enabled_missing.yaml"
            },
            {
                "duplicate slicing enabled flag",
                "  kt:\n"
                "    enabled: false\n",
                "  kt:\n"
                "    enabled: false\n"
                "    enabled: false\n",
                "hbt_config_loader_test_slicing_enabled_duplicate.yaml"
            },
            {
                "unknown slicing-axis key",
                "  kt:\n"
                "    enabled: false\n",
                "  kt:\n"
                "    enabled: false\n"
                "    unknown: 1\n",
                "hbt_config_loader_test_slicing_axis_unknown.yaml"
            },
            {
                "non-boolean slicing enabled flag",
                "  kt:\n"
                "    enabled: false\n",
                "  kt:\n"
                "    enabled: not_a_boolean\n",
                "hbt_config_loader_test_slicing_enabled_invalid.yaml"
            },
            {
                "enabled slicing without bin edges",
                "  kt:\n"
                "    enabled: false\n",
                "  kt:\n"
                "    enabled: true\n",
                "hbt_config_loader_test_slicing_edges_missing.yaml"
            },
            {
                "non-sequence slicing bin edges",
                "  kt:\n"
                "    enabled: false\n",
                "  kt:\n"
                "    enabled: true\n"
                "    bin_edges_gev: 0.2\n",
                "hbt_config_loader_test_slicing_edges_non_sequence.yaml"
            },
            {
                "non-scalar origin mode",
                "hbt_origin_mode: \"all\"\n",
                "hbt_origin_mode:\n"
                "  - all\n",
                "hbt_config_loader_test_origin_non_scalar.yaml"
            }
        };

        for (const ReplacementCase& test_case : cases) {
            std::string contents =
                make_valid_config();

            if (
                !replace_once(
                    contents,
                    test_case.original,
                    test_case.replacement
                )
            ) {
                std::cerr
                    << "hbt_config_loader_test: replacement failed for "
                    << test_case.test_name
                    << ".\n";
                return false;
            }

            if (
                !verify_runtime_error(
                    test_case.test_name,
                    test_case.filename,
                    contents
                )
            ) {
                return false;
            }
        }

        return true;
    }

    /**
     * @brief Verify scientific-value validation by YAML replacement.
     *
     * @return true when every invalid scientific value is rejected.
     */
    bool verify_invalid_argument_cases() {
        const ReplacementCase cases[] = {
            {
                "invalid enabled channels",
                "\"pi_plus_pi_plus,k_plus_p+k_minus_p_bar\"",
                "\"pi_plus_pi_plus,,k_plus_p\"",
                "hbt_config_loader_test_invalid_channels.yaml"
            },
            {
                "invalid longitudinal variable",
                "\"pseudorapidity\"",
                "\"Pseudorapidity\"",
                "hbt_config_loader_test_invalid_longitudinal.yaml"
            },
            {
                "invalid origin mode",
                "hbt_origin_mode: \"all\"\n",
                "hbt_origin_mode: \"Primordial\"\n",
                "hbt_config_loader_test_invalid_origin.yaml"
            },
            {
                "zero longitudinal cut",
                "      longitudinal_abs_max: 0.8\n",
                "      longitudinal_abs_max: 0.0\n",
                "hbt_config_loader_test_zero_longitudinal.yaml"
            },
            {
                "negative pt minimum",
                "      pt_min_gev: 0.14\n",
                "      pt_min_gev: -0.1\n",
                "hbt_config_loader_test_negative_pt_min.yaml"
            },
            {
                "invalid pt interval",
                "      pt_max_gev: 4.0\n",
                "      pt_max_gev: 0.14\n",
                "hbt_config_loader_test_invalid_pt_interval.yaml"
            },
            {
                "non-numeric cut",
                "      pt_min_gev: 0.14\n",
                "      pt_min_gev: not_a_number\n",
                "hbt_config_loader_test_non_numeric_cut.yaml"
            },
            {
                "non-finite cut",
                "      pt_max_gev: 4.0\n",
                "      pt_max_gev: .nan\n",
                "hbt_config_loader_test_non_finite_cut.yaml"
            },
            {
                "zero OSL histogram bin count",
                "    nbins: 30\n",
                "    nbins: 0\n",
                "hbt_config_loader_test_zero_osl_nbins.yaml"
            },
            {
                "fractional OSL histogram bin count",
                "    nbins: 30\n",
                "    nbins: 30.5\n",
                "hbt_config_loader_test_fractional_osl_nbins.yaml"
            },
            {
                "negative OSL histogram minimum",
                "    min_fm: 0.0\n"
                "    max_fm: 3.0\n",
                "    min_fm: -1.0\n"
                "    max_fm: 3.0\n",
                "hbt_config_loader_test_negative_osl_min.yaml"
            },
            {
                "positive radial histogram minimum",
                "  radial:\n"
                "    nbins: 20\n"
                "    min_fm: 0.0\n",
                "  radial:\n"
                "    nbins: 20\n"
                "    min_fm: 0.1\n",
                "hbt_config_loader_test_positive_radial_min.yaml"
            },
            {
                "delta-t histogram does not span negative values",
                "    min_fm_c: -2.0\n"
                "    max_fm_c: 2.0\n",
                "    min_fm_c: 0.0\n"
                "    max_fm_c: 2.0\n",
                "hbt_config_loader_test_delta_t_no_negative.yaml"
            },
            {
                "delta-t histogram does not span positive values",
                "    min_fm_c: -2.0\n"
                "    max_fm_c: 2.0\n",
                "    min_fm_c: -2.0\n"
                "    max_fm_c: 0.0\n",
                "hbt_config_loader_test_delta_t_no_positive.yaml"
            },
            {
                "invalid retained bin edges on disabled slicing",
                "  kt:\n"
                "    enabled: false\n",
                "  kt:\n"
                "    enabled: false\n"
                "    bin_edges_gev: [0.2]\n",
                "hbt_config_loader_test_slicing_disabled_edges_invalid.yaml"
            },
            {
                "too few slicing bin edges",
                "  kt:\n"
                "    enabled: false\n",
                "  kt:\n"
                "    enabled: true\n"
                "    bin_edges_gev: [0.2]\n",
                "hbt_config_loader_test_slicing_edges_too_few.yaml"
            },
            {
                "negative slicing bin edge",
                "  kt:\n"
                "    enabled: false\n",
                "  kt:\n"
                "    enabled: true\n"
                "    bin_edges_gev: [-0.1, 0.2]\n",
                "hbt_config_loader_test_slicing_edge_negative.yaml"
            },
            {
                "duplicate slicing bin edge",
                "  kt:\n"
                "    enabled: false\n",
                "  kt:\n"
                "    enabled: true\n"
                "    bin_edges_gev: [0.2, 0.2]\n",
                "hbt_config_loader_test_slicing_edge_duplicate.yaml"
            },
            {
                "decreasing slicing bin edges",
                "  kt:\n"
                "    enabled: false\n",
                "  kt:\n"
                "    enabled: true\n"
                "    bin_edges_gev: [0.4, 0.2]\n",
                "hbt_config_loader_test_slicing_edges_decreasing.yaml"
            },
            {
                "non-numeric slicing bin edge",
                "  kt:\n"
                "    enabled: false\n",
                "  kt:\n"
                "    enabled: true\n"
                "    bin_edges_gev: [0.2, bad]\n",
                "hbt_config_loader_test_slicing_edge_non_numeric.yaml"
            },
            {
                "non-finite slicing bin edge",
                "  kt:\n"
                "    enabled: false\n",
                "  kt:\n"
                "    enabled: true\n"
                "    bin_edges_gev: [0.2, .nan]\n",
                "hbt_config_loader_test_slicing_edge_non_finite.yaml"
            }
        };

        for (const ReplacementCase& test_case : cases) {
            std::string contents =
                make_valid_config();

            if (
                !replace_once(
                    contents,
                    test_case.original,
                    test_case.replacement
                )
            ) {
                std::cerr
                    << "hbt_config_loader_test: replacement failed for "
                    << test_case.test_name
                    << ".\n";
                return false;
            }

            if (
                !verify_invalid_argument(
                    test_case.test_name,
                    test_case.filename,
                    contents
                )
            ) {
                return false;
            }
        }

        return true;
    }

}  // namespace

/**
 * @brief Run the HBT configuration-loader unit tests.
 *
 * @return EXIT_SUCCESS when every test passes, otherwise EXIT_FAILURE.
 */
int main() {
    bool success = true;

    if (!verify_pair_slicing_combinations()) {
        success = false;
    }

    if (!verify_valid_tokens()) {
        success = false;
    }

    if (!verify_unquoted_enabled_channels()) {
        success = false;
    }

    if (!verify_folded_enabled_channels()) {
        success = false;
    }

    if (!verify_direct_runtime_errors()) {
        success = false;
    }

    if (!verify_replacement_runtime_errors()) {
        success = false;
    }

    if (!verify_invalid_argument_cases()) {
        success = false;
    }

    return success ? EXIT_SUCCESS : EXIT_FAILURE;
}
