/**
 * @file event_preparation_integration_test.cpp
 * @brief Integration tests for HBT preparation and pair-count orchestration.
 */

#include "app/analysis_runner.h"
#include "output/analysis_output_writer.h"

#include "hbt/channels/primitive_channel.h"
#include "hbt/species/species.h"

#include <algorithm>
#include <array>
#include <cmath>
#include <cstdint>
#include <cstdlib>
#include <exception>
#include <filesystem>
#include <fstream>
#include <iostream>
#include <map>
#include <sstream>
#include <stdexcept>
#include <string>
#include <system_error>
#include <thread>

namespace {

/**
 * @brief Return the temporary root used by event-preparation integration tests.
 * @return Relative temporary-directory path.
 */
std::filesystem::path test_root() {
    return "event_preparation_integration_test_data";
}

/**
 * @brief Prepare an empty temporary root for integration fixtures.
 * @return `true` when the directory is ready, otherwise `false`.
 */
bool prepare_test_root() {
    std::error_code error;
    static_cast<void>(std::filesystem::remove_all(test_root(), error));
    error.clear();
    static_cast<void>(std::filesystem::create_directories(test_root(), error));

    if (error) {
        std::cerr
            << "event_preparation_integration_test: failed to prepare root.\n";
        return false;
    }

    return true;
}

/**
 * @brief Remove the temporary integration-test root without throwing.
 */
void cleanup_test_root() noexcept {
    std::error_code error;
    static_cast<void>(std::filesystem::remove_all(test_root(), error));
}

/**
 * @brief Write one complete text fixture, creating parent directories.
 * @param path Path of the file to create.
 * @param content Complete text content.
 * @return `true` when the file is written successfully.
 */
bool write_text_file(
    const std::filesystem::path& path,
    const std::string& content
) {
    std::error_code error;
    static_cast<void>(
        std::filesystem::create_directories(path.parent_path(), error)
    );

    if (error) {
        std::cerr
            << "event_preparation_integration_test: failed to create parent.\n";
        return false;
    }

    std::ofstream output(path);

    if (!output) {
        std::cerr
            << "event_preparation_integration_test: failed to create file.\n";
        return false;
    }

    output << content;

    if (!output) {
        std::cerr
            << "event_preparation_integration_test: failed to write file.\n";
        return false;
    }

    return true;
}

/**
 * @brief Return the supported Afterburner three-line file header.
 * @return Complete header with trailing newlines.
 */
std::string afterburner_header() {
    return
        "#!OSCAR2013Extended particle_lists "
        "t x y z mass p0 px py pz pdg ID charge ncoll form_time xsecfac "
        "proc_id_origin proc_type_origin time_last_coll pdg_mother1 "
        "pdg_mother2 baryon_number strangeness\n"
        "# Units: fm fm fm fm GeV GeV GeV GeV GeV none none e none fm "
        "none none none fm none none none none\n"
        "# SMASH-3.2.1\n";
}

/**
 * @brief Return the supported Sampler three-line file header.
 * @return Complete header with trailing newlines.
 */
std::string sampler_header() {
    return
        "#!OSCAR2013 particle_lists t x y z mass p0 px py pz pdg ID charge\n"
        "# Units: fm fm fm fm GeV GeV GeV GeV GeV none none e\n"
        "# SMASH-3.2.1\n";
}

/**
 * @brief Build one subevent opening marker.
 * @param subevent_id Subevent identifier.
 * @param particle_count Declared particle-row count.
 * @return Complete opening marker.
 */
std::string opening_marker(int subevent_id, std::size_t particle_count) {
    return "# event " + std::to_string(subevent_id) +
           " ensemble 0 out " + std::to_string(particle_count) + "\n";
}

/**
 * @brief Build one valid subevent closing marker.
 * @param subevent_id Subevent identifier.
 * @return Complete closing marker.
 */
std::string closing_marker(int subevent_id) {
    return "# event " + std::to_string(subevent_id) +
           " ensemble 0 end 0 impact 0.000 "
           "scattering_projectile_target yes\n";
}

/**
 * @brief Build one complete Afterburner particle row for integration tests.
 * @param pdg Raw signed PDG code.
 * @param id Raw particle identifier.
 * @param charge Raw electric charge.
 * @param ncoll Raw collision count.
 * @param px First transverse momentum component in GeV.
 * @param py Second transverse momentum component in GeV.
 * @param time_last_coll Last-interaction time in fm.
 * @param mother1 First raw mother-PDG field.
 * @param mother2 Second raw mother-PDG field.
 * @return Complete 22-column row with trailing newline.
 */
std::string afterburner_row(
    int pdg,
    int id,
    int charge,
    int ncoll,
    double px,
    double py,
    double time_last_coll,
    int mother1,
    int mother2
) {
    return
        "10 5 7 9 0.138 2 " + std::to_string(px) + " " +
        std::to_string(py) + " 0 " + std::to_string(pdg) + " " +
        std::to_string(id) + " " + std::to_string(charge) + " " +
        std::to_string(ncoll) + " 0 1 0 0 " +
        std::to_string(time_last_coll) + " " +
        std::to_string(mother1) + " " + std::to_string(mother2) +
        " 0 0\n";
}

/**
 * @brief Build one Afterburner row with controlled energy for mass tests.
 * @param energy Raw p0 component in GeV.
 * @param px First transverse momentum component in GeV.
 * @param id Raw particle identifier.
 * @return Primordial positive-pion row with zero longitudinal momentum.
 */
std::string afterburner_mass_test_row(
    double energy,
    double px,
    int id
) {
    return
        "10 5 7 9 0.138 " + std::to_string(energy) + " " +
        std::to_string(px) + " 0 0 211 " + std::to_string(id) +
        " 1 0 0 1 0 0 4 0 0 0 0\n";
}

/**
 * @brief Build one row whose finite propagation inputs overflow in position.
 * @param id Raw particle identifier.
 * @return Non-primordial positive-pion row with finite input components.
 */
std::string afterburner_propagation_overflow_row(int id) {
    constexpr const char* huge = "1.7976931348623157e308";

    return
        std::string(huge) + " " + huge +
        " 0 0 0.138 2 -1 0 0 211 " + std::to_string(id) +
        " 1 1 0 1 0 0 0 1 2 0 0\n";
}

/**
 * @brief Build one complete Sampler row for a lookup key.
 * @param t Sampler time coordinate in fm.
 * @param id Raw particle identifier.
 * @param pdg Raw signed PDG code.
 * @param charge Raw electric charge.
 * @return Complete 12-column Sampler row.
 */
std::string sampler_row(double t, int id, int pdg, int charge) {
    return std::to_string(t) +
           " 1 2 3 0.138 9 0.1 0.2 0.3 " +
           std::to_string(pdg) + " " + std::to_string(id) + " " +
           std::to_string(charge) + "\n";
}

/**
 * @brief Return the default configuration with both pair-slicing axes off.
 * @return Complete hbt_pair_slicing YAML block.
 */
std::string disabled_pair_slicing_text() {
    return
        "hbt_pair_slicing:\n"
        "  kt:\n"
        "    enabled: false\n"
        "  mt:\n"
        "    enabled: false\n";
}

/**
 * @brief Return one kT-only slicing configuration for runner integration.
 * @return Complete hbt_pair_slicing YAML block with three kT bins.
 */
std::string kt_pair_slicing_text() {
    return
        "hbt_pair_slicing:\n"
        "  kt:\n"
        "    enabled: true\n"
        "    bin_edges_gev: [0.5, 0.6, 0.7, 0.8]\n"
        "  mt:\n"
        "    enabled: false\n";
}

/**
 * @brief Return one complete HBT configuration for two required species.
 * @param origin_mode Canonical configured origin-mode token.
 * @param pair_slicing Complete validated pair-slicing YAML block.
 * @return Complete valid HBT YAML text.
 */
std::string hbt_config_text(
    const std::string& origin_mode = "all",
    const std::string& pair_slicing = disabled_pair_slicing_text()
) {
    return
        "hbt_enabled_channels: \"pi_plus_pi_plus,k_plus_k_plus\"\n"
        "\n"
        "hbt_particle_acceptance:\n"
        "  longitudinal_variable: \"pseudorapidity\"\n"
        "  groups:\n"
        "    pions:\n"
        "      longitudinal_abs_max: 0.8\n"
        "      pt_min_gev: 0.14\n"
        "      pt_max_gev: 4.0\n"
        "    kaons:\n"
        "      longitudinal_abs_max: 0.8\n"
        "      pt_min_gev: 0.4\n"
        "      pt_max_gev: 1.4\n"
        "    nucleons:\n"
        "      longitudinal_abs_max: 0.8\n"
        "      pt_min_gev: 0.5\n"
        "      pt_max_gev: 4.05\n"
        "    sigmas:\n"
        "      longitudinal_abs_max: 0.8\n"
        "      pt_min_gev: 1.0\n"
        "      pt_max_gev: 10000.0\n"
        "    lambdas:\n"
        "      longitudinal_abs_max: 0.8\n"
        "      pt_min_gev: 0.3\n"
        "      pt_max_gev: 10000.0\n" +
        pair_slicing +
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
        "hbt_origin_mode: \"" + origin_mode + "\"\n";
}

/**
 * @brief Return one enabled-HBT global configuration for a fixture case.
 * @param number_of_events Configured outer-event count.
 * @param number_of_subevents Configured per-event subevent count.
 * @return Complete valid main YAML text.
 */
std::string main_config_text(
    std::size_t number_of_events,
    std::size_t number_of_subevents,
    std::size_t threads = 1U
) {
    return
        "events_path: \"../events\"\n"
        "output_path: \"output\"\n"
        "number_of_events: " + std::to_string(number_of_events) + "\n"
        "number_of_subevents: " +
        std::to_string(number_of_subevents) + "\n"
        "threads: " + std::to_string(threads) + "\n"
        "hbt_enabled: true\n"
        "hbt_config_path: \"hbt.yaml\"\n";
}

/**
 * @brief Write the common configuration pair for one enabled-HBT case.
 * @param case_root Root directory of the case.
 * @param events Configured outer-event count.
 * @param subevents Configured per-event subevent count.
 * @param origin_mode Canonical origin-mode token.
 * @param pair_slicing Complete pair-slicing YAML block.
 * @return `true` when both configuration files are written.
 */
bool write_configs(
    const std::filesystem::path& case_root,
    std::size_t events,
    std::size_t subevents,
    const std::string& origin_mode = "all",
    const std::string& pair_slicing = disabled_pair_slicing_text(),
    std::size_t threads = 1U
) {
    const std::filesystem::path config_dir = case_root / "config";

    return write_text_file(
               config_dir / "main.yaml",
               main_config_text(events, subevents, threads)
           ) &&
           write_text_file(
               config_dir / "hbt.yaml",
               hbt_config_text(origin_mode, pair_slicing)
           );
}

/**
 * @brief Write one outer event's complete Afterburner and Sampler files.
 * @param case_root Root directory of the case.
 * @param event_number One-based outer-event number.
 * @param afterburner Complete Afterburner file text.
 * @param sampler Complete Sampler file text.
 * @return `true` when both files are written.
 */
bool write_outer_event(
    const std::filesystem::path& case_root,
    std::size_t event_number,
    const std::string& afterburner,
    const std::string& sampler
) {
    const std::filesystem::path event_dir =
        case_root / "events" / std::to_string(event_number);

    return write_text_file(
               event_dir / "Afterburner" / "particle_lists.oscar",
               afterburner
           ) &&
           write_text_file(
               event_dir / "Sampler" / "particle_lists.oscar",
               sampler
           );
}

/**
 * @brief Find one species counter in a completed preparation summary.
 * @param summary Completed HBT preparation summary.
 * @param species Canonical species to find.
 * @return Pointer to its counters, or nullptr if absent.
 */
const app::HBTSpeciesPreparationCounts* find_species_counts(
    const app::HBTEventPreparationSummary& summary,
    hbt::SpeciesId species
) {
    for (const auto& entry : summary.species) {
        if (entry.species == species) {
            return &entry;
        }
    }

    return nullptr;
}

/**
 * @brief Verify one ordered two-channel pair-count summary.
 * @param summary Pair counts to inspect.
 * @param pion_pairs Expected pi_plus_pi_plus count.
 * @param kaon_pairs Expected k_plus_k_plus count.
 * @return true when both channels and counts match exactly.
 */
bool verify_two_channel_pair_counts(
    const hbt::PairCountSummary& summary,
    std::uint64_t pion_pairs,
    std::uint64_t kaon_pairs
) {
    return summary.channels.size() == 2U &&
           summary.channels[0].channel ==
               hbt::PrimitiveChannelId::PiPlusPiPlus &&
           summary.channels[0].pair_count == pion_pairs &&
           summary.channels[1].channel ==
               hbt::PrimitiveChannelId::KPlusKPlus &&
           summary.channels[1].pair_count == kaon_pairs;
}

/**
 * @brief Test whether every counter in one raw histogram set is zero.
 * @param histogram Histogram set to inspect.
 * @return true only when bins, underflows, and overflows are all zero.
 */
bool raw_histogram_set_is_zero(const hbt::RawHistogramSet& histogram) {
    const auto vector_is_zero = [](const std::vector<std::uint64_t>& values) {
        return std::all_of(
            values.begin(),
            values.end(),
            [](std::uint64_t value) { return value == 0U; }
        );
    };
    const auto array4_is_zero = [](
        const std::array<std::uint64_t, 4U>& values
    ) {
        return std::all_of(
            values.begin(),
            values.end(),
            [](std::uint64_t value) { return value == 0U; }
        );
    };
    const auto array3_is_zero = [](
        const std::array<std::uint64_t, 3U>& values
    ) {
        return std::all_of(
            values.begin(),
            values.end(),
            [](std::uint64_t value) { return value == 0U; }
        );
    };
    const auto array2_is_zero = [](
        const std::array<std::uint64_t, 2U>& values
    ) {
        return std::all_of(
            values.begin(),
            values.end(),
            [](std::uint64_t value) { return value == 0U; }
        );
    };

    return vector_is_zero(histogram.osl.bins) &&
           array4_is_zero(histogram.osl.underflow_counts) &&
           array4_is_zero(histogram.osl.overflow_counts) &&
           vector_is_zero(histogram.radial.bins) &&
           array2_is_zero(histogram.radial.underflow_counts) &&
           array2_is_zero(histogram.radial.overflow_counts) &&
           vector_is_zero(histogram.delta_t.bins) &&
           array3_is_zero(histogram.delta_t.underflow_counts) &&
           array3_is_zero(histogram.delta_t.overflow_counts);
}

/**
 * @brief Test whether one complete raw-histogram state contains only zeros.
 * @param state Complete raw-histogram state to inspect.
 * @return true only when every global and sliced destination is zero.
 */
bool raw_histograms_are_zero(const hbt::RawHistogramState& state) {
    for (const hbt::ProductRawHistogramState& product : state.products) {
        for (const hbt::OriginRawHistogramState& origin : product.origins) {
            if (!raw_histogram_set_is_zero(origin.global)) {
                return false;
            }
            for (const hbt::RawHistogramSet& slice : origin.slices) {
                if (!raw_histogram_set_is_zero(slice)) {
                    return false;
                }
            }
        }
    }
    return true;
}

/**
 * @brief Verify the complete ordered multi-event preparation pipeline.
 *
 * The fixture includes unsupported and unrequired particles, a required pion
 * that fails kinematic acceptance while lacking a Sampler key, all three
 * origin-membership patterns, all three emission-position sources, and fresh
 * Sampler routing across two outer events.
 *
 * @return `true` when all exact counters match expectation.
 */
bool verify_complete_pipeline() {
    const std::filesystem::path root = test_root() / "complete";

    if (!write_configs(root, 2U, 2U)) {
        return false;
    }

    std::string afterburner1 = afterburner_header();
    afterburner1 += opening_marker(0, 4U);
    afterburner1 += afterburner_row(211, 101, 1, 0, 0.5, 0.0, 4.0, 0, 0);
    afterburner1 += afterburner_row(321, 102, 1, 1, 0.7, 0.0, 4.0, 1, 2);
    afterburner1 += afterburner_row(2212, 103, 1, 0, 0.8, 0.0, 4.0, 0, 0);
    afterburner1 += afterburner_row(211, 104, 1, 0, 0.1, 0.0, 4.0, 0, 0);
    afterburner1 += closing_marker(0);
    afterburner1 += opening_marker(1, 3U);
    afterburner1 += afterburner_row(211, 105, 1, 1, 0.5, 0.0, 3.0, 1, 0);
    afterburner1 += afterburner_row(321, 106, 1, 0, 0.7, 0.0, 4.0, 0, 0);
    afterburner1 += afterburner_row(13, 107, -1, 0, 0.5, 0.0, 4.0, 0, 0);
    afterburner1 += closing_marker(1);

    std::string sampler1 = sampler_header();
    sampler1 += opening_marker(0, 1U);
    sampler1 += sampler_row(1.0, 101, 211, 1);
    sampler1 += closing_marker(0);
    sampler1 += opening_marker(1, 1U);
    sampler1 += sampler_row(2.0, 106, 321, 1);
    sampler1 += closing_marker(1);

    if (!write_outer_event(root, 1U, afterburner1, sampler1)) {
        return false;
    }

    std::string afterburner2 = afterburner_header();
    afterburner2 += opening_marker(0, 1U);
    afterburner2 += afterburner_row(211, 201, 1, 0, 0.5, 0.0, 4.0, 0, 0);
    afterburner2 += closing_marker(0);
    afterburner2 += opening_marker(1, 1U);
    afterburner2 += afterburner_row(321, 202, 1, 1, 0.7, 0.0, -1.0, 1, 0);
    afterburner2 += closing_marker(1);

    std::string sampler2 = sampler_header();
    sampler2 += opening_marker(0, 1U);
    sampler2 += sampler_row(3.0, 201, 211, 1);
    sampler2 += closing_marker(0);
    sampler2 += opening_marker(1, 0U);
    sampler2 += closing_marker(1);

    if (!write_outer_event(root, 2U, afterburner2, sampler2)) {
        return false;
    }

    try {
        const app::AnalysisRunSummary result = app::AnalysisRunner(
            root / "config" / "main.yaml"
        ).run();

        if (!result.hbt_event_preparation.has_value() ||
            !result.hbt_pair_processing.has_value()) {
            std::cerr
                << "event_preparation_integration_test: missing HBT summary.\n";
            return false;
        }

        const app::HBTEventPreparationSummary& summary =
            result.hbt_event_preparation.value();

        if (summary.outer_events_processed != 2U ||
            summary.subevents_processed != 4U ||
            summary.raw_particles != 9U ||
            summary.unsupported_species != 1U ||
            summary.unrequired_species != 1U ||
            summary.particle_acceptance_rejections != 1U ||
            !summary.numerical_rejections.empty() ||
            summary.origin_rejections != 0U ||
            summary.accepted_particles != 6U) {
            std::cerr
                << "event_preparation_integration_test: pipeline counts "
                << "differ.\n";
            return false;
        }

        if (summary.emission_points.sampler != 3U ||
            summary.emission_points.propagation != 2U ||
            summary.emission_points.afterburner != 1U) {
            std::cerr
                << "event_preparation_integration_test: source counts "
                << "differ.\n";
            return false;
        }

        const auto* pion = find_species_counts(summary, hbt::SpeciesId::PiPlus);
        const auto* kaon = find_species_counts(summary, hbt::SpeciesId::KPlus);

        if (pion == nullptr || kaon == nullptr ||
            pion->accepted != 3U || pion->primordial != 2U ||
            pion->primordial_rescattering != 2U ||
            pion->primordial_rescattering_decay != 3U ||
            kaon->accepted != 3U || kaon->primordial != 1U ||
            kaon->primordial_rescattering != 2U ||
            kaon->primordial_rescattering_decay != 3U) {
            std::cerr
                << "event_preparation_integration_test: origin counts "
                << "differ.\n";
            return false;
        }

        if (summary.subevents.size() != 4U ||
            summary.subevents[0].outer_event_number != 1U ||
            summary.subevents[0].subevent_id != 0 ||
            summary.subevents[0].accepted_particles != 2U ||
            summary.subevents[1].outer_event_number != 1U ||
            summary.subevents[1].subevent_id != 1 ||
            summary.subevents[1].accepted_particles != 2U ||
            summary.subevents[2].outer_event_number != 2U ||
            summary.subevents[2].subevent_id != 0 ||
            summary.subevents[2].accepted_particles != 1U ||
            summary.subevents[3].outer_event_number != 2U ||
            summary.subevents[3].subevent_id != 1 ||
            summary.subevents[3].accepted_particles != 1U) {
            std::cerr
                << "event_preparation_integration_test: subevent counts "
                << "differ.\n";
            return false;
        }

        const hbt::HBTPairProcessingSummary& pair_summary =
            result.hbt_pair_processing.value();

        if (!verify_two_channel_pair_counts(
                pair_summary.total_pair_counts, 0U, 0U) ||
            !verify_two_channel_pair_counts(
                pair_summary.total_valid_pair_counts, 0U, 0U) ||
            !verify_two_channel_pair_counts(
                pair_summary.total_numerical_rejection_counts, 0U, 0U) ||
            pair_summary.total_origin_route_counts.origin_mode !=
                hbt::OriginMode::All ||
            !verify_two_channel_pair_counts(
                pair_summary.total_origin_route_counts.routed_P, 0U, 0U) ||
            !verify_two_channel_pair_counts(
                pair_summary.total_origin_route_counts.routed_PR, 0U, 0U) ||
            !verify_two_channel_pair_counts(
                pair_summary.total_origin_route_counts.routed_PRD, 0U, 0U) ||
            !pair_summary.numerical_rejections.empty() ||
            pair_summary.subevents.size() != 4U) {
            std::cerr
                << "event_preparation_integration_test: zero-pair "
                << "integration counts differ.\n";
            return false;
        }

        return true;
    } catch (const std::exception& error) {
        std::cerr
            << "event_preparation_integration_test: complete case threw: "
            << error.what() << ".\n";
        return false;
    }
}

/**
 * @brief Verify pair enumeration crosses the runner boundary correctly.
 * @return true when local and run-total identical-species counts match.
 */
bool verify_pair_count_integration() {
    const std::filesystem::path root = test_root() / "pair_counts";

    if (!write_configs(root, 1U, 1U)) {
        return false;
    }

    std::string afterburner = afterburner_header();
    afterburner += opening_marker(0, 5U);
    afterburner +=
        afterburner_row(211, 401, 1, 1, 0.5, 0.0, 3.0, 1, 0);
    afterburner +=
        afterburner_row(211, 402, 1, 1, 0.6, 0.0, 3.0, 1, 0);
    afterburner +=
        afterburner_row(211, 403, 1, 1, 0.7, 0.0, 3.0, 1, 0);
    afterburner +=
        afterburner_row(321, 404, 1, 1, 0.7, 0.0, 3.0, 1, 0);
    afterburner +=
        afterburner_row(321, 405, 1, 1, 0.8, 0.0, 3.0, 1, 0);
    afterburner += closing_marker(0);

    std::string sampler = sampler_header();
    sampler += opening_marker(0, 0U);
    sampler += closing_marker(0);

    if (!write_outer_event(root, 1U, afterburner, sampler)) {
        return false;
    }

    try {
        const app::AnalysisRunSummary result = app::AnalysisRunner(
            root / "config" / "main.yaml"
        ).run();

        if (!result.hbt_event_preparation.has_value() ||
            !result.hbt_pair_processing.has_value()) {
            return false;
        }

        if (result.hbt_event_preparation->accepted_particles != 5U) {
            return false;
        }

        const hbt::HBTPairProcessingSummary& pair_summary =
            result.hbt_pair_processing.value();

        if (!verify_two_channel_pair_counts(
                pair_summary.total_pair_counts, 3U, 1U) ||
            !verify_two_channel_pair_counts(
                pair_summary.total_valid_pair_counts, 3U, 1U) ||
            !verify_two_channel_pair_counts(
                pair_summary.total_numerical_rejection_counts, 0U, 0U) ||
            pair_summary.total_origin_route_counts.origin_mode !=
                hbt::OriginMode::All ||
            !verify_two_channel_pair_counts(
                pair_summary.total_origin_route_counts.routed_P, 0U, 0U) ||
            !verify_two_channel_pair_counts(
                pair_summary.total_origin_route_counts.routed_PR, 0U, 0U) ||
            !verify_two_channel_pair_counts(
                pair_summary.total_origin_route_counts.routed_PRD, 3U, 1U) ||
            !pair_summary.total_pair_slice_counts.entries.empty() ||
            !pair_summary.numerical_rejections.empty()) {
            return false;
        }

        if (pair_summary.subevents.size() != 1U ||
            pair_summary.subevents[0].outer_event_number != 1U ||
            pair_summary.subevents[0].subevent_id != 0 ||
            !verify_two_channel_pair_counts(
                pair_summary.subevents[0].pair_counts, 3U, 1U) ||
            !verify_two_channel_pair_counts(
                pair_summary.subevents[0].valid_pair_counts, 3U, 1U) ||
            !verify_two_channel_pair_counts(
                pair_summary.subevents[0].numerical_rejection_counts,
                0U,
                0U) ||
            pair_summary.subevents[0].origin_route_counts.origin_mode !=
                hbt::OriginMode::All ||
            !verify_two_channel_pair_counts(
                pair_summary.subevents[0].origin_route_counts.routed_P,
                0U,
                0U) ||
            !verify_two_channel_pair_counts(
                pair_summary.subevents[0].origin_route_counts.routed_PR,
                0U,
                0U) ||
            !verify_two_channel_pair_counts(
                pair_summary.subevents[0].origin_route_counts.routed_PRD,
                3U,
                1U) ||
            !pair_summary.subevents[0].pair_slice_counts.entries.empty()) {
            return false;
        }

        return true;
    } catch (const std::exception& error) {
        std::cerr
            << "event_preparation_integration_test: pair-count case threw: "
            << error.what() << ".\n";
        return false;
    }
}

/**
 * @brief Verify configured kT slicing crosses pair processing and runner.
 * @return true when two subevents accumulate exact slice/origin/channel counts.
 */
bool verify_pair_slice_routing_crosses_runner_boundary() {
    const std::filesystem::path root = test_root() / "pair_slices";

    if (!write_configs(
            root,
            1U,
            2U,
            "all",
            kt_pair_slicing_text())) {
        return false;
    }

    std::string afterburner = afterburner_header();
    afterburner += opening_marker(0, 3U);
    afterburner +=
        afterburner_row(211, 601, 1, 1, 0.51, 0.0, 3.0, 1, 0);
    afterburner +=
        afterburner_row(211, 602, 1, 1, 0.63, 0.0, 3.0, 1, 0);
    afterburner +=
        afterburner_row(211, 603, 1, 1, 0.72, 0.0, 3.0, 1, 0);
    afterburner += closing_marker(0);
    afterburner += opening_marker(1, 2U);
    afterburner +=
        afterburner_row(321, 604, 1, 1, 0.72, 0.0, 3.0, 1, 0);
    afterburner +=
        afterburner_row(321, 605, 1, 1, 0.78, 0.0, 3.0, 1, 0);
    afterburner += closing_marker(1);

    std::string sampler = sampler_header();
    sampler += opening_marker(0, 0U);
    sampler += closing_marker(0);
    sampler += opening_marker(1, 0U);
    sampler += closing_marker(1);

    if (!write_outer_event(root, 1U, afterburner, sampler)) {
        return false;
    }

    try {
        const app::AnalysisRunSummary result = app::AnalysisRunner(
            root / "config" / "main.yaml"
        ).run();

        if (!result.hbt_pair_processing.has_value()) {
            return false;
        }

        const hbt::HBTPairProcessingSummary& pair_summary =
            result.hbt_pair_processing.value();
        const hbt::PairSliceCountSummary& slices =
            pair_summary.total_pair_slice_counts;

        if (!verify_two_channel_pair_counts(
                pair_summary.total_valid_pair_counts, 3U, 1U) ||
            !verify_two_channel_pair_counts(
                pair_summary.total_origin_route_counts.routed_P, 0U, 0U) ||
            !verify_two_channel_pair_counts(
                pair_summary.total_origin_route_counts.routed_PR, 0U, 0U) ||
            !verify_two_channel_pair_counts(
                pair_summary.total_origin_route_counts.routed_PRD, 3U, 1U) ||
            slices.origin_mode != hbt::OriginMode::All ||
            slices.entries.size() != 3U ||
            pair_summary.subevents.size() != 2U) {
            return false;
        }

        if (slices.entries[0].kt_slice_index != 0U ||
            slices.entries[1].kt_slice_index != 1U ||
            slices.entries[2].kt_slice_index != 2U ||
            slices.entries[0].mt_slice_index.has_value() ||
            slices.entries[1].mt_slice_index.has_value() ||
            slices.entries[2].mt_slice_index.has_value()) {
            return false;
        }

        const hbt::PairSliceOriginCounts& slice0 =
            slices.entries[0].origin_counts;
        const hbt::PairSliceOriginCounts& slice1 =
            slices.entries[1].origin_counts;
        const hbt::PairSliceOriginCounts& slice2 =
            slices.entries[2].origin_counts;

        if (!verify_two_channel_pair_counts(slice0.routed_P, 0U, 0U) ||
            !verify_two_channel_pair_counts(slice0.routed_PR, 0U, 0U) ||
            !verify_two_channel_pair_counts(slice0.routed_PRD, 1U, 0U) ||
            !verify_two_channel_pair_counts(slice1.routed_P, 0U, 0U) ||
            !verify_two_channel_pair_counts(slice1.routed_PR, 0U, 0U) ||
            !verify_two_channel_pair_counts(slice1.routed_PRD, 2U, 0U) ||
            !verify_two_channel_pair_counts(slice2.routed_P, 0U, 0U) ||
            !verify_two_channel_pair_counts(slice2.routed_PR, 0U, 0U) ||
            !verify_two_channel_pair_counts(slice2.routed_PRD, 0U, 1U)) {
            return false;
        }

        const hbt::PairSliceCountSummary& local0 =
            pair_summary.subevents[0].pair_slice_counts;
        const hbt::PairSliceCountSummary& local1 =
            pair_summary.subevents[1].pair_slice_counts;

        return local0.entries.size() == 3U &&
               local1.entries.size() == 3U &&
               verify_two_channel_pair_counts(
                   local0.entries[0].origin_counts.routed_PRD,
                   1U,
                   0U) &&
               verify_two_channel_pair_counts(
                   local0.entries[1].origin_counts.routed_PRD,
                   2U,
                   0U) &&
               verify_two_channel_pair_counts(
                   local1.entries[2].origin_counts.routed_PRD,
                   0U,
                   1U);
    } catch (const std::exception& error) {
        std::cerr
            << "event_preparation_integration_test: pair-slice case threw: "
            << error.what() << ".\n";
        return false;
    }
}

/**
 * @brief Verify configured origin mode reaches production pair routing.
 * @return true when P-P is routed only to the selected individual route.
 */
bool verify_origin_mode_crosses_runner_boundary() {
    const struct {
        const char* directory;
        const char* mode;
        std::uint64_t p;
        std::uint64_t pr;
        std::uint64_t prd;
    } cases[]{
        {"origin_mode_p", "primordial", 1U, 0U, 0U},
        {"origin_mode_pr", "primordial_rescattering", 0U, 1U, 0U},
        {
            "origin_mode_prd",
            "primordial_rescattering_decay",
            0U,
            0U,
            1U
        }
    };

    for (const auto& test_case : cases) {
        const std::filesystem::path root =
            test_root() / test_case.directory;

        if (!write_configs(root, 1U, 1U, test_case.mode)) {
            return false;
        }

        std::string afterburner = afterburner_header();
        afterburner += opening_marker(0, 2U);
        afterburner +=
            afterburner_row(211, 501, 1, 1, 0.5, 0.0, 3.0, 0, 0);
        afterburner +=
            afterburner_row(211, 502, 1, 1, 0.6, 0.0, 3.0, 0, 0);
        afterburner += closing_marker(0);

        std::string sampler = sampler_header();
        sampler += opening_marker(0, 0U);
        sampler += closing_marker(0);

        if (!write_outer_event(root, 1U, afterburner, sampler)) {
            return false;
        }

        try {
            const app::AnalysisRunSummary result = app::AnalysisRunner(
                root / "config" / "main.yaml"
            ).run();

            if (!result.hbt_pair_processing.has_value()) {
                return false;
            }

            const hbt::HBTPairProcessingSummary& pair_summary =
                result.hbt_pair_processing.value();

            if (!verify_two_channel_pair_counts(
                    pair_summary.total_valid_pair_counts, 1U, 0U) ||
                !verify_two_channel_pair_counts(
                    pair_summary.total_origin_route_counts.routed_P,
                    test_case.p,
                    0U) ||
                !verify_two_channel_pair_counts(
                    pair_summary.total_origin_route_counts.routed_PR,
                    test_case.pr,
                    0U) ||
                !verify_two_channel_pair_counts(
                    pair_summary.total_origin_route_counts.routed_PRD,
                    test_case.prd,
                    0U)) {
                std::cerr
                    << "event_preparation_integration_test: configured "
                    << "origin mode did not reach pair routing.\n";
                return false;
            }
        } catch (const std::exception& error) {
            std::cerr
                << "event_preparation_integration_test: origin-mode pair "
                << "routing case threw: " << error.what() << ".\n";
            return false;
        }
    }

    return true;
}

/**
 * @brief Test observer verifying callbacks remain on the orchestration thread.
 */
class RecordingProgressObserver final : public app::AnalysisProgressObserver {
public:
    RecordingProgressObserver()
        : owner_thread_(std::this_thread::get_id()) {
    }

    void begin(
        std::size_t total_events,
        std::size_t total_subevents
    ) override {
        record_thread();
        total_events_ = total_events;
        total_subevents_ = total_subevents;
        if (total_events == 0U ||
            total_subevents % total_events != 0U) {
            canonical_order_ = false;
            return;
        }
        subevents_per_event_ = total_subevents / total_events;
    }

    void subevent_completed(
        std::size_t completed_subevents,
        std::size_t outer_event_number,
        std::size_t subevent_number
    ) override {
        record_thread();
        if (completed_subevents != completed_subevents_ + 1U ||
            total_events_ == 0U ||
            total_subevents_ % total_events_ != 0U) {
            canonical_order_ = false;
        } else {
            const std::size_t subevents_per_event =
                total_subevents_ / total_events_;
            const std::size_t index = completed_subevents - 1U;
            if (outer_event_number != index / subevents_per_event + 1U ||
                subevent_number != index % subevents_per_event + 1U) {
                canonical_order_ = false;
            }
        }
        completed_subevents_ = completed_subevents;
        last_outer_event_ = outer_event_number;
        last_subevent_ = subevent_number;
    }

    void begin_postprocessing() override {
        record_thread();
    }

    void analysis_complete() override {
        record_thread();
    }

    void begin_output() override {
        record_thread();
    }

    void finish() override {
        record_thread();
    }

    void fail() noexcept override {
        if (std::this_thread::get_id() != owner_thread_) {
            wrong_thread_ = true;
        }
    }

    [[nodiscard]] bool complete_for(
        std::size_t total_events,
        std::size_t total_subevents
    ) const noexcept {
        return !wrong_thread_ &&
               canonical_order_ &&
               total_events_ == total_events &&
               total_subevents_ == total_subevents &&
               completed_subevents_ == total_subevents &&
               last_outer_event_ == total_events &&
               last_subevent_ == subevents_per_event_;
    }

private:
    void record_thread() noexcept {
        if (std::this_thread::get_id() != owner_thread_) {
            wrong_thread_ = true;
        }
    }

    std::thread::id owner_thread_;
    std::size_t total_events_{0U};
    std::size_t total_subevents_{0U};
    std::size_t subevents_per_event_{0U};
    std::size_t completed_subevents_{0U};
    std::size_t last_outer_event_{0U};
    std::size_t last_subevent_{0U};
    bool wrong_thread_{false};
    bool canonical_order_{true};
};

/**
 * @brief Observer that deliberately fails on the first subevent callback.
 */
class ThrowingProgressObserver final : public app::AnalysisProgressObserver {
public:
    void begin(std::size_t, std::size_t) override {
    }

    void subevent_completed(
        std::size_t,
        std::size_t,
        std::size_t
    ) override {
        throw std::runtime_error("test progress observer failure");
    }

    void begin_postprocessing() override {
    }

    void analysis_complete() override {
    }

    void begin_output() override {
    }

    void finish() override {
    }

    void fail() noexcept override {
    }
};

/**
 * @brief Serialize one complete in-memory run for exact comparison.
 * @param result Completed run result.
 * @return Stable analysis-output text.
 */
std::string serialized_analysis(const app::AnalysisRunSummary& result) {
    std::ostringstream output_stream;
    output::write_analysis_output(result, output_stream);
    return output_stream.str();
}

/**
 * @brief Snapshot every production-output file as exact bytes by relative path.
 * @param root Existing production-output root.
 * @return Ordered relative-path to byte-content mapping.
 * @throws std::runtime_error If any regular file cannot be opened.
 */
std::map<std::string, std::string> output_tree_snapshot(
    const std::filesystem::path& root
) {
    std::map<std::string, std::string> snapshot;
    for (const auto& entry :
         std::filesystem::recursive_directory_iterator(root)) {
        if (!entry.is_regular_file()) {
            continue;
        }
        std::ifstream input(entry.path(), std::ios::binary);
        if (!input) {
            throw std::runtime_error(
                "failed to open production output for exact comparison"
            );
        }
        std::ostringstream bytes;
        bytes << input.rdbuf();
        snapshot.emplace(
            std::filesystem::relative(entry.path(), root).generic_string(),
            bytes.str()
        );
    }
    return snapshot;
}

/**
 * @brief Verify exact serial/parallel equivalence and canonical merge order.
 * @return true for serial, explicit N, automatic, and capped worker counts.
 */
bool verify_parallel_event_equivalence() {
    const std::filesystem::path root = test_root() / "parallel_equivalence";
    if (!write_configs(root, 4U, 2U)) {
        return false;
    }

    const std::array<std::size_t, 4U> particle_counts{40U, 2U, 5U, 3U};
    for (std::size_t event = 0U; event < particle_counts.size(); ++event) {
        std::string afterburner = afterburner_header();
        afterburner += opening_marker(0, particle_counts[event]);
        for (std::size_t particle = 0U;
             particle < particle_counts[event];
             ++particle) {
            afterburner += afterburner_row(
                211,
                static_cast<int>(1000U + event * 100U + particle),
                1,
                1,
                0.5 + 0.01 * static_cast<double>(particle),
                0.0,
                3.0,
                1,
                0
            );
        }
        afterburner += closing_marker(0);
        afterburner += opening_marker(1, 0U);
        afterburner += closing_marker(1);

        std::string sampler = sampler_header();
        sampler += opening_marker(0, 0U);
        sampler += closing_marker(0);
        sampler += opening_marker(1, 0U);
        sampler += closing_marker(1);
        if (!write_outer_event(root, event + 1U, afterburner, sampler)) {
            return false;
        }
    }

    try {
        const std::filesystem::path main_path = root / "config" / "main.yaml";
        const app::AnalysisRunSummary serial =
            app::AnalysisRunner(main_path).run();
        const std::string serial_output = serialized_analysis(serial);
        output::write_production_output(serial);
        const auto serial_production = output_tree_snapshot(
            serial.startup.run_config.output_path
        );
        std::filesystem::remove_all(serial.startup.run_config.output_path);

        std::optional<app::AnalysisRunSummary> parallel;
        const std::array<std::size_t, 4U> thread_settings{2U, 4U, 0U, 16U};
        for (const std::size_t threads : thread_settings) {
            if (!write_text_file(
                    main_path,
                    main_config_text(4U, 2U, threads))) {
                return false;
            }

            if (threads == 4U) {
                RecordingProgressObserver progress;
                parallel.emplace(
                    app::AnalysisRunner(main_path).run(&progress)
                );
                if (!progress.complete_for(4U, 8U)) {
                    return false;
                }
            } else {
                parallel.emplace(app::AnalysisRunner(main_path).run());
            }

            if (serial_output != serialized_analysis(parallel.value())) {
                return false;
            }

            if (threads == 4U) {
                output::write_production_output(parallel.value());
                if (serial_production != output_tree_snapshot(
                        parallel->startup.run_config.output_path)) {
                    return false;
                }
                std::filesystem::remove_all(
                    parallel->startup.run_config.output_path
                );
            }
        }

        if (!parallel.has_value() ||
            !parallel->hbt_event_preparation.has_value()) {
            return false;
        }
        const auto& preparation =
            parallel->hbt_event_preparation.value();
        if (preparation.events.size() != 4U ||
            preparation.subevents.size() != 8U) {
            return false;
        }
        for (std::size_t event = 0U; event < 4U; ++event) {
            if (preparation.events[event].outer_event_number != event + 1U) {
                return false;
            }
            for (std::size_t subevent = 0U; subevent < 2U; ++subevent) {
                const std::size_t index = event * 2U + subevent;
                if (preparation.subevents[index].outer_event_number !=
                        event + 1U ||
                    preparation.subevents[index].subevent_id !=
                        static_cast<int>(subevent)) {
                    return false;
                }
            }
        }
        return true;
    } catch (const std::exception& error) {
        std::cerr
            << "event_preparation_integration_test: parallel equivalence "
            << "threw: " << error.what() << ".\n";
        return false;
    }
}

/**
 * @brief Verify an event-local input failure is skipped in parallel execution.
 * @return true when the failed event contributes zero and the run continues.
 */
bool verify_parallel_event_failure_recovers() {
    const std::filesystem::path root =
        test_root() / "parallel_event_failure";
    if (!write_configs(
            root,
            2U,
            1U,
            "all",
            disabled_pair_slicing_text(),
            2U)) {
        return false;
    }

    std::string afterburner = afterburner_header();
    afterburner += opening_marker(0, 0U);
    afterburner += closing_marker(0);
    std::string valid_sampler = sampler_header();
    valid_sampler += opening_marker(0, 0U);
    valid_sampler += closing_marker(0);
    if (!write_outer_event(root, 1U, afterburner, valid_sampler) ||
        !write_outer_event(root, 2U, afterburner, sampler_header())) {
        return false;
    }

    try {
        RecordingProgressObserver progress;
        const app::AnalysisRunSummary result = app::AnalysisRunner(
            root / "config" / "main.yaml"
        ).run(&progress);
        if (!progress.complete_for(2U, 2U) ||
            !result.hbt_event_preparation.has_value()) {
            return false;
        }

        const auto& summary = result.hbt_event_preparation.value();
        const std::string output_text = serialized_analysis(result);
        return summary.outer_events_processed == 2U &&
               summary.subevents_processed == 1U &&
               summary.event_status_counts.processed == 0U &&
               summary.event_status_counts.skipped_due_to_event_empty == 1U &&
               summary.event_status_counts
                       .skipped_due_to_event_failure == 1U &&
               summary.events.size() == 2U &&
               summary.events[1].status ==
                   app::EventStatus::SkippedDueToEventFailure &&
               summary.events[1].diagnostic.find(
                   "Sampler subevent count mismatch") != std::string::npos &&
               output_text.find("skipped_due_to_event_failure") !=
                   std::string::npos &&
               output_text.find("Sampler subevent count mismatch") !=
                   std::string::npos;
    } catch (const std::exception& error) {
        std::cerr
            << "event_preparation_integration_test: parallel event failure "
            << "threw: " << error.what() << ".\n";
        return false;
    }
}

/**
 * @brief Verify a subevent-local Sampler miss is skipped in parallel execution.
 * @return true when the failed subevent contributes zero and the run continues.
 */
bool verify_parallel_subevent_failure_recovers() {
    const std::filesystem::path root =
        test_root() / "parallel_subevent_failure";
    if (!write_configs(
            root,
            2U,
            1U,
            "all",
            disabled_pair_slicing_text(),
            2U)) {
        return false;
    }

    std::string empty_afterburner = afterburner_header();
    empty_afterburner += opening_marker(0, 0U);
    empty_afterburner += closing_marker(0);
    std::string missing_afterburner = afterburner_header();
    missing_afterburner += opening_marker(0, 1U);
    missing_afterburner += afterburner_row(
        211, 9001, 1, 0, 0.5, 0.0, 3.0, 0, 0
    );
    missing_afterburner += closing_marker(0);

    std::string sampler = sampler_header();
    sampler += opening_marker(0, 0U);
    sampler += closing_marker(0);
    if (!write_outer_event(root, 1U, empty_afterburner, sampler) ||
        !write_outer_event(root, 2U, missing_afterburner, sampler)) {
        return false;
    }

    try {
        RecordingProgressObserver progress;
        const app::AnalysisRunSummary result = app::AnalysisRunner(
            root / "config" / "main.yaml"
        ).run(&progress);
        if (!progress.complete_for(2U, 2U) ||
            !result.hbt_event_preparation.has_value()) {
            return false;
        }

        const auto& summary = result.hbt_event_preparation.value();
        const std::string output_text = serialized_analysis(result);
        return summary.outer_events_processed == 2U &&
               summary.subevents_processed == 2U &&
               summary.accepted_particles == 0U &&
               summary.event_status_counts.skipped_due_to_event_empty == 1U &&
               summary.event_status_counts.processed == 1U &&
               summary.event_status_counts
                       .skipped_due_to_event_failure == 0U &&
               summary.subevent_status_counts
                       .skipped_due_to_subevent_empty == 1U &&
               summary.subevent_status_counts
                       .skipped_due_to_subevent_failure == 1U &&
               summary.subevents.size() == 2U &&
               summary.subevents[1].status ==
                   app::SubeventStatus::SkippedDueToSubeventFailure &&
               summary.subevents[1].accepted_particles == 0U &&
               summary.subevents[1].diagnostic.find(
                   "missing mandatory Sampler match") != std::string::npos &&
               output_text.find("skipped_due_to_subevent_failure") !=
                   std::string::npos &&
               output_text.find("missing mandatory Sampler match") !=
                   std::string::npos;
    } catch (const std::exception& error) {
        std::cerr
            << "event_preparation_integration_test: parallel subevent "
            << "failure threw: " << error.what() << ".\n";
        return false;
    }
}

/**
 * @brief Verify observer exceptions cannot leave joinable workers behind.
 * @return true when the observer failure propagates normally from the run.
 */
bool verify_parallel_progress_exception_is_safe() {
    const std::filesystem::path root = test_root() / "progress_exception";
    if (!write_configs(
            root,
            2U,
            1U,
            "all",
            disabled_pair_slicing_text(),
            2U)) {
        return false;
    }

    std::string afterburner = afterburner_header();
    afterburner += opening_marker(0, 0U);
    afterburner += closing_marker(0);
    std::string sampler = sampler_header();
    sampler += opening_marker(0, 0U);
    sampler += closing_marker(0);
    if (!write_outer_event(root, 1U, afterburner, sampler) ||
        !write_outer_event(root, 2U, afterburner, sampler)) {
        return false;
    }

    ThrowingProgressObserver progress;
    try {
        static_cast<void>(app::AnalysisRunner(
            root / "config" / "main.yaml"
        ).run(&progress));
    } catch (const std::runtime_error& error) {
        return std::string(error.what()) == "test progress observer failure";
    } catch (...) {
        return false;
    }
    return false;
}

/**
 * @brief Verify explicit empty event/subevent status reporting.
 * @return true when empty inputs contribute zero and remain visible.
 */
bool verify_empty_execution_statuses() {
    const std::filesystem::path root = test_root() / "empty_statuses";
    if (!write_configs(root, 2U, 2U)) {
        return false;
    }

    std::string sampler = sampler_header();
    sampler += opening_marker(0, 0U);
    sampler += closing_marker(0);
    sampler += opening_marker(1, 0U);
    sampler += closing_marker(1);

    std::string empty_event = afterburner_header();
    empty_event += opening_marker(0, 0U);
    empty_event += closing_marker(0);
    empty_event += opening_marker(1, 0U);
    empty_event += closing_marker(1);

    std::string mixed_event = afterburner_header();
    mixed_event += opening_marker(0, 1U);
    mixed_event += afterburner_row(211, 2001, 1, 1, 0.5, 0.0, 3.0, 1, 0);
    mixed_event += closing_marker(0);
    mixed_event += opening_marker(1, 0U);
    mixed_event += closing_marker(1);

    if (!write_outer_event(root, 1U, empty_event, sampler) ||
        !write_outer_event(root, 2U, mixed_event, sampler)) {
        return false;
    }

    try {
        const app::AnalysisRunSummary result = app::AnalysisRunner(
            root / "config" / "main.yaml"
        ).run();
        if (!result.hbt_event_preparation.has_value()) {
            return false;
        }
        const auto& summary = result.hbt_event_preparation.value();
        const std::string output_text = serialized_analysis(result);
        if (output_text.find("skipped_due_to_event_empty") ==
                std::string::npos ||
            output_text.find("skipped_due_to_subevent_empty") ==
                std::string::npos) {
            return false;
        }
        return summary.event_status_counts.processed == 1U &&
               summary.event_status_counts.skipped_due_to_event_empty == 1U &&
               summary.event_status_counts.skipped_due_to_event_failure == 0U &&
               summary.subevent_status_counts.processed == 1U &&
               summary.subevent_status_counts
                       .skipped_due_to_subevent_empty == 3U &&
               summary.subevent_status_counts
                       .skipped_due_to_subevent_failure == 0U &&
               summary.events.size() == 2U &&
               summary.events[0].status ==
                   app::EventStatus::SkippedDueToEventEmpty &&
               summary.events[1].status == app::EventStatus::Processed &&
               summary.subevents.size() == 4U &&
               summary.subevents[0].status ==
                   app::SubeventStatus::SkippedDueToSubeventEmpty &&
               summary.subevents[1].status ==
                   app::SubeventStatus::SkippedDueToSubeventEmpty &&
               summary.subevents[2].status == app::SubeventStatus::Processed &&
               summary.subevents[3].status ==
                   app::SubeventStatus::SkippedDueToSubeventEmpty;
    } catch (const std::exception& error) {
        std::cerr
            << "event_preparation_integration_test: empty status case threw: "
            << error.what() << ".\n";
        return false;
    }
}

/**
 * @brief Verify disabled HBT never accesses the configured events path.
 * @return `true` when run succeeds with a missing event directory.
 */
bool verify_disabled_hbt_skips_event_access() {
    const std::filesystem::path root = test_root() / "disabled";
    const std::filesystem::path config_dir = root / "config";

    if (!write_text_file(
            config_dir / "main.yaml",
            "events_path: \"missing-events\"\n"
            "output_path: \"output\"\n"
            "number_of_events: 1\n"
            "number_of_subevents: 1\n"
            "hbt_enabled: false\n")) {
        return false;
    }

    try {
        const app::AnalysisRunSummary result = app::AnalysisRunner(
            config_dir / "main.yaml"
        ).run();

        return !result.hbt_event_preparation.has_value() &&
               !result.hbt_pair_processing.has_value();
    } catch (const std::exception& error) {
        std::cerr
            << "event_preparation_integration_test: disabled case threw: "
            << error.what() << ".\n";
        return false;
    }
}

/**
 * @brief Verify origin rejection occurs before any mandatory Sampler lookup.
 * @return `true` when a non-primordial particle is rejected without a miss.
 */
bool verify_origin_rejection_precedes_sampler() {
    const std::filesystem::path root = test_root() / "origin_reject";

    if (!write_configs(root, 1U, 1U, "primordial")) {
        return false;
    }

    std::string afterburner = afterburner_header();
    afterburner += opening_marker(0, 1U);
    afterburner += afterburner_row(211, 301, 1, 0, 0.5, 0.0, 4.0, 1, 2);
    afterburner += closing_marker(0);

    std::string sampler = sampler_header();
    sampler += opening_marker(0, 0U);
    sampler += closing_marker(0);

    if (!write_outer_event(root, 1U, afterburner, sampler)) {
        return false;
    }

    try {
        const app::AnalysisRunSummary result = app::AnalysisRunner(
            root / "config" / "main.yaml"
        ).run();
        const auto& summary = result.hbt_event_preparation.value();

        return summary.origin_rejections == 1U &&
               summary.accepted_particles == 0U &&
               summary.emission_points.sampler == 0U &&
               summary.emission_points.propagation == 0U &&
               summary.emission_points.afterburner == 0U;
    } catch (const std::exception& error) {
        std::cerr
            << "event_preparation_integration_test: origin case threw: "
            << error.what() << ".\n";
        return false;
    }
}

/**
 * @brief Verify non-positive energy is rejected before any Sampler lookup.
 * @return true when the particle is recorded and never reaches emission logic.
 */
bool verify_nonpositive_energy_rejection_precedes_sampler() {
    const std::filesystem::path root = test_root() / "energy_reject";

    if (!write_configs(root, 1U, 1U)) {
        return false;
    }

    std::string afterburner = afterburner_header();
    afterburner += opening_marker(0, 1U);
    afterburner += afterburner_mass_test_row(-2.0, 0.5, 350);
    afterburner += closing_marker(0);

    std::string sampler = sampler_header();
    sampler += opening_marker(0, 0U);
    sampler += closing_marker(0);

    if (!write_outer_event(root, 1U, afterburner, sampler)) {
        return false;
    }

    try {
        const app::AnalysisRunSummary result = app::AnalysisRunner(
            root / "config" / "main.yaml"
        ).run();
        const auto& summary = result.hbt_event_preparation.value();
        const auto& report = summary.numerical_rejections;

        if (summary.accepted_particles != 0U ||
            summary.particle_acceptance_rejections != 0U ||
            summary.origin_rejections != 0U || report.size() != 1U ||
            summary.emission_points.sampler != 0U ||
            summary.emission_points.propagation != 0U ||
            summary.emission_points.afterburner != 0U) {
            std::cerr
                << "event_preparation_integration_test: energy rejection "
                << "counts differ.\n";
            return false;
        }

        const auto& record = report.records().front();

        return
            record.outer_event_number == 1U &&
            record.subevent_id == 0 &&
            record.particle_id == 350 &&
            record.raw_pdg == 211 &&
            record.raw_charge == 1 &&
            record.species == hbt::SpeciesId::PiPlus &&
            record.reason == hbt::ParticleRejectionReason::
                NonPositiveEnergy &&
            record.diagnostic_value.has_value() &&
            record.diagnostic_value.value() == -2.0;
    } catch (const std::exception& error) {
        std::cerr
            << "event_preparation_integration_test: energy rejection threw: "
            << error.what() << ".\n";
        return false;
    }
}

/**
 * @brief Verify invalid invariant mass is reported before Sampler lookup.
 * @return true when one particle is rejected, recorded, and does not query
 *         its otherwise mandatory missing Sampler key.
 */
bool verify_mass_rejection_precedes_sampler() {
    const std::filesystem::path root = test_root() / "mass_reject";

    if (!write_configs(root, 1U, 1U)) {
        return false;
    }

    std::string afterburner = afterburner_header();
    afterburner += opening_marker(0, 1U);
    afterburner += afterburner_mass_test_row(0.4, 0.5, 351);
    afterburner += closing_marker(0);

    std::string sampler = sampler_header();
    sampler += opening_marker(0, 0U);
    sampler += closing_marker(0);

    if (!write_outer_event(root, 1U, afterburner, sampler)) {
        return false;
    }

    try {
        const app::AnalysisRunSummary result = app::AnalysisRunner(
            root / "config" / "main.yaml"
        ).run();
        const auto& summary = result.hbt_event_preparation.value();
        const auto& report = summary.numerical_rejections;

        if (summary.accepted_particles != 0U ||
            summary.particle_acceptance_rejections != 0U ||
            summary.origin_rejections != 0U || report.size() != 1U ||
            summary.emission_points.sampler != 0U ||
            summary.emission_points.propagation != 0U ||
            summary.emission_points.afterburner != 0U) {
            std::cerr
                << "event_preparation_integration_test: mass rejection "
                << "counts differ.\n";
            return false;
        }

        const auto& record = report.records().front();

        return
            record.outer_event_number == 1U &&
            record.subevent_id == 0 &&
            record.particle_id == 351 &&
            record.raw_pdg == 211 &&
            record.raw_charge == 1 &&
            record.species == hbt::SpeciesId::PiPlus &&
            record.reason == hbt::ParticleRejectionReason::
                NonPositiveInvariantMassSquared &&
            record.diagnostic_value.has_value() &&
            record.diagnostic_value.value() < 0.0;
    } catch (const std::exception& error) {
        std::cerr
            << "event_preparation_integration_test: mass rejection threw: "
            << error.what() << ".\n";
        return false;
    }
}

/**
 * @brief Verify non-finite propagation rejects only the affected particle.
 * @return true when the run succeeds and records the exact propagation fault.
 */
bool verify_nonfinite_propagation_is_reported() {
    const std::filesystem::path root = test_root() / "propagation_reject";

    if (!write_configs(root, 1U, 1U)) {
        return false;
    }

    std::string afterburner = afterburner_header();
    afterburner += opening_marker(0, 1U);
    afterburner += afterburner_propagation_overflow_row(352);
    afterburner += closing_marker(0);

    std::string sampler = sampler_header();
    sampler += opening_marker(0, 0U);
    sampler += closing_marker(0);

    if (!write_outer_event(root, 1U, afterburner, sampler)) {
        return false;
    }

    try {
        const app::AnalysisRunSummary result = app::AnalysisRunner(
            root / "config" / "main.yaml"
        ).run();
        const auto& summary = result.hbt_event_preparation.value();
        const auto& report = summary.numerical_rejections;

        if (summary.accepted_particles != 0U || report.size() != 1U ||
            summary.emission_points.sampler != 0U ||
            summary.emission_points.propagation != 0U ||
            summary.emission_points.afterburner != 0U) {
            std::cerr
                << "event_preparation_integration_test: propagation "
                << "rejection counts differ.\n";
            return false;
        }

        const auto& record = report.records().front();

        return
            record.outer_event_number == 1U &&
            record.subevent_id == 0 &&
            record.particle_id == 352 &&
            record.raw_pdg == 211 &&
            record.raw_charge == 1 &&
            record.species == hbt::SpeciesId::PiPlus &&
            record.reason == hbt::ParticleRejectionReason::
                NonFinitePropagationEmissionPosition &&
            record.diagnostic_position.has_value() &&
            !std::isfinite(record.diagnostic_position.value().x1);
    } catch (const std::exception& error) {
        std::cerr
            << "event_preparation_integration_test: propagation rejection "
            << "threw: " << error.what() << ".\n";
        return false;
    }
}

/**
 * @brief Verify a mandatory Sampler miss skips only its current subevent.
 * @return true when the failed subevent contributes zero and the next runs.
 */
bool verify_mandatory_sampler_miss_is_subevent_failure() {
    const std::filesystem::path root = test_root() / "missing_sampler";

    if (!write_configs(root, 1U, 2U)) {
        return false;
    }

    std::string afterburner = afterburner_header();
    afterburner += opening_marker(0, 1U);
    afterburner += afterburner_row(211, 401, 1, 0, 0.5, 0.0, 4.0, 0, 0);
    afterburner += closing_marker(0);
    afterburner += opening_marker(1, 1U);
    afterburner += afterburner_row(211, 402, 1, 1, 0.5, 0.0, 3.0, 1, 0);
    afterburner += closing_marker(1);

    std::string sampler = sampler_header();
    sampler += opening_marker(0, 0U);
    sampler += closing_marker(0);
    sampler += opening_marker(1, 0U);
    sampler += closing_marker(1);

    if (!write_outer_event(root, 1U, afterburner, sampler)) {
        return false;
    }

    try {
        RecordingProgressObserver progress;
        const app::AnalysisRunSummary result = app::AnalysisRunner(
            root / "config" / "main.yaml"
        ).run(&progress);
        if (!progress.complete_for(1U, 2U) ||
            !result.hbt_event_preparation.has_value() ||
            !result.hbt_pair_processing.has_value()) {
            return false;
        }

        const auto& summary = result.hbt_event_preparation.value();
        const auto& pairs = result.hbt_pair_processing.value();
        const std::string output_text = serialized_analysis(result);
        return summary.outer_events_processed == 1U &&
               summary.subevents_processed == 2U &&
               summary.raw_particles == 1U &&
               summary.accepted_particles == 1U &&
               summary.event_status_counts.processed == 1U &&
               summary.subevent_status_counts.processed == 1U &&
               summary.subevent_status_counts
                       .skipped_due_to_subevent_failure == 1U &&
               summary.subevents.size() == 2U &&
               summary.subevents[0].status ==
                   app::SubeventStatus::SkippedDueToSubeventFailure &&
               summary.subevents[0].accepted_particles == 0U &&
               summary.subevents[0].diagnostic.find(
                   "missing mandatory Sampler match") != std::string::npos &&
               summary.subevents[1].status ==
                   app::SubeventStatus::Processed &&
               summary.subevents[1].accepted_particles == 1U &&
               pairs.subevents.size() == 1U &&
               pairs.subevents[0].subevent_id == 1 &&
               output_text.find("missing mandatory Sampler match") !=
                   std::string::npos;
    } catch (const std::exception& error) {
        std::cerr
            << "event_preparation_integration_test: mandatory miss recovery "
            << "threw: " << error.what() << ".\n";
        return false;
    }
}

/**
 * @brief Verify a late event-local cardinality failure rolls the event back.
 * @return true when prior subevent science is discarded and the run continues.
 */
bool verify_event_failure_rolls_back_prior_subevents() {
    const std::filesystem::path root =
        test_root() / "event_failure_rollback";

    if (!write_configs(root, 1U, 2U)) {
        return false;
    }

    std::string afterburner = afterburner_header();
    afterburner += opening_marker(0, 2U);
    afterburner += afterburner_row(211, 501, 1, 1, 0.5, 0.0, 3.0, 1, 0);
    afterburner += afterburner_row(211, 502, 1, 1, 0.6, 0.0, 3.0, 1, 0);
    afterburner += closing_marker(0);

    std::string sampler = sampler_header();
    sampler += opening_marker(0, 0U);
    sampler += closing_marker(0);
    sampler += opening_marker(1, 0U);
    sampler += closing_marker(1);

    if (!write_outer_event(root, 1U, afterburner, sampler)) {
        return false;
    }

    try {
        RecordingProgressObserver progress;
        const app::AnalysisRunSummary result = app::AnalysisRunner(
            root / "config" / "main.yaml"
        ).run(&progress);
        if (!progress.complete_for(1U, 2U) ||
            !result.hbt_event_preparation.has_value() ||
            !result.hbt_pair_processing.has_value() ||
            !result.hbt_raw_histograms.has_value()) {
            return false;
        }

        const auto& summary = result.hbt_event_preparation.value();
        const auto& pairs = result.hbt_pair_processing.value();
        const std::string output_text = serialized_analysis(result);
        return summary.outer_events_processed == 1U &&
               summary.subevents_processed == 0U &&
               summary.raw_particles == 0U &&
               summary.accepted_particles == 0U &&
               summary.event_status_counts.processed == 0U &&
               summary.event_status_counts
                       .skipped_due_to_event_failure == 1U &&
               summary.subevent_status_counts.processed == 0U &&
               summary.subevents.empty() &&
               summary.events.size() == 1U &&
               summary.events[0].status ==
                   app::EventStatus::SkippedDueToEventFailure &&
               summary.events[0].diagnostic.find(
                   "Afterburner ended before configured") !=
                   std::string::npos &&
               verify_two_channel_pair_counts(
                   pairs.total_pair_counts,
                   0U,
                   0U
               ) &&
               pairs.subevents.empty() &&
               raw_histograms_are_zero(result.hbt_raw_histograms.value()) &&
               output_text.find("skipped_due_to_event_failure") !=
                   std::string::npos &&
               output_text.find("Afterburner ended before configured") !=
                   std::string::npos;
    } catch (const std::exception& error) {
        std::cerr
            << "event_preparation_integration_test: event rollback threw: "
            << error.what() << ".\n";
        return false;
    }
}

}  // namespace

/**
 * @brief Run the event-preparation and pair-count integration collection.
 * @return EXIT_SUCCESS when every integration case passes.
 */
int main() {
    if (!prepare_test_root()) {
        return EXIT_FAILURE;
    }

    bool success = true;
    success = verify_complete_pipeline() && success;
    success = verify_parallel_event_equivalence() && success;
    success = verify_parallel_event_failure_recovers() && success;
    success = verify_parallel_subevent_failure_recovers() && success;
    success = verify_parallel_progress_exception_is_safe() && success;
    success = verify_empty_execution_statuses() && success;
    success = verify_pair_count_integration() && success;
    success = verify_pair_slice_routing_crosses_runner_boundary() && success;
    success = verify_origin_mode_crosses_runner_boundary() && success;
    success = verify_disabled_hbt_skips_event_access() && success;
    success = verify_origin_rejection_precedes_sampler() && success;
    success = verify_nonpositive_energy_rejection_precedes_sampler() &&
              success;
    success = verify_mass_rejection_precedes_sampler() && success;
    success = verify_nonfinite_propagation_is_reported() && success;
    success = verify_mandatory_sampler_miss_is_subevent_failure() &&
              success;
    success = verify_event_failure_rolls_back_prior_subevents() &&
              success;

    cleanup_test_root();
    return success ? EXIT_SUCCESS : EXIT_FAILURE;
}
