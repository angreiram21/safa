/**
 * @file event_preparation_real_fixture_test.cpp
 * @brief Regression test for real-fixture HBT preparation and pair counts.
 */

#include "app/analysis_runner.h"

#include "hbt/channels/primitive_channel.h"
#include "hbt/config/hbt_config.h"
#include "hbt/species/species.h"

#include <cstdint>
#include <cstdlib>
#include <exception>
#include <filesystem>
#include <iostream>
#include <limits>
#include <string>
#include <vector>

namespace {

/**
 * @brief Compare one integer regression value and print a diagnostic on error.
 * @param label Human-readable name of the value being checked.
 * @param actual Value produced by the modular event-preparation run.
 * @param expected Frozen real-fixture reference value.
 * @return `true` when actual equals expected.
 */
bool expect_count(
    const std::string& label,
    std::size_t actual,
    std::size_t expected
) {
    if (actual == expected) {
        return true;
    }

    std::cerr
        << "event_preparation_real_fixture_test: "
        << label
        << " expected "
        << expected
        << ", found "
        << actual
        << ".\n";
    return false;
}

/**
 * @brief Compare one 64-bit pair-count regression value.
 * @param label Human-readable name of the value being checked.
 * @param actual Pair count produced by the modular run.
 * @param expected Frozen real-fixture pair count validated against the legacy.
 * @return `true` when actual equals expected.
 */
bool expect_pair_count(
    const std::string& label,
    std::uint64_t actual,
    std::uint64_t expected
) {
    if (actual == expected) {
        return true;
    }

    std::cerr
        << "event_preparation_real_fixture_test: "
        << label
        << " expected "
        << expected
        << ", found "
        << actual
        << ".\n";
    return false;
}

/**
 * @brief Compare one particle-acceptance cut group exactly.
 * @param label Human-readable particle-group name.
 * @param actual Resolved cuts loaded by AnalysisRunner.
 * @param longitudinal_abs_max Expected longitudinal absolute-value limit.
 * @param pt_min_gev Expected lower transverse-momentum limit in GeV.
 * @param pt_max_gev Expected upper transverse-momentum limit in GeV.
 * @return `true` when all three resolved values match the frozen config.
 */
bool expect_cuts(
    const char* label,
    const hbt::ParticleAcceptanceCuts& actual,
    double longitudinal_abs_max,
    double pt_min_gev,
    double pt_max_gev
) {
    if (actual.longitudinal_abs_max == longitudinal_abs_max &&
        actual.pt_min_gev == pt_min_gev &&
        actual.pt_max_gev == pt_max_gev) {
        return true;
    }

    std::cerr
        << "event_preparation_real_fixture_test: "
        << label
        << " acceptance config differs from the frozen fixture config.\n";
    return false;
}

/**
 * @brief Check whether one SpeciesId occurs in a species vector.
 * @param species Collection of canonical species identifiers.
 * @param expected Species identifier that must occur.
 * @return `true` when expected occurs in species.
 */
bool contains_species(
    const std::vector<hbt::SpeciesId>& species,
    hbt::SpeciesId expected
) {
    for (const hbt::SpeciesId value : species) {
        if (value == expected) {
            return true;
        }
    }

    return false;
}

/**
 * @brief Verify the frozen scientific and run configuration for the fixture.
 * @param result Completed application run containing resolved startup state.
 * @return `true` when the resolved configuration matches the regression setup.
 */
bool verify_startup_contract(const app::AnalysisRunSummary& result) {
    bool success = true;

    success = expect_count(
        "number_of_events",
        result.startup.run_config.number_of_events,
        1U
    ) && success;
    success = expect_count(
        "number_of_subevents",
        result.startup.run_config.number_of_subevents,
        2000U
    ) && success;

    if (!result.startup.run_config.hbt_enabled ||
        !result.startup.hbt_config.has_value() ||
        !result.startup.hbt_startup_state.has_value()) {
        std::cerr
            << "event_preparation_real_fixture_test: HBT startup is not "
            << "fully enabled and resolved.\n";
        return false;
    }

    const hbt::HBTConfig& config = result.startup.hbt_config.value();

    if (config.origin_mode != hbt::OriginMode::All) {
        std::cerr
            << "event_preparation_real_fixture_test: origin mode is not All.\n";
        success = false;
    }

    if (config.particle_acceptance.longitudinal_variable !=
        hbt::LongitudinalVariable::Pseudorapidity) {
        std::cerr
            << "event_preparation_real_fixture_test: longitudinal variable "
            << "is not pseudorapidity.\n";
        success = false;
    }

    success = expect_cuts(
        "pion",
        config.particle_acceptance.pions,
        0.8,
        0.14,
        4.0
    ) && success;
    success = expect_cuts(
        "kaon",
        config.particle_acceptance.kaons,
        0.8,
        0.4,
        1.4
    ) && success;
    success = expect_cuts(
        "nucleon",
        config.particle_acceptance.nucleons,
        0.8,
        0.5,
        4.05
    ) && success;
    success = expect_cuts(
        "sigma",
        config.particle_acceptance.sigmas,
        0.8,
        1.0,
        10000.0
    ) && success;
    success = expect_cuts(
        "lambda",
        config.particle_acceptance.lambdas,
        0.8,
        0.3,
        10000.0
    ) && success;

    const hbt::HBTStartupState& startup =
        result.startup.hbt_startup_state.value();

    const std::vector<hbt::PrimitiveChannelId> expected_channels{
        hbt::PrimitiveChannelId::PiPlusPiPlus,
        hbt::PrimitiveChannelId::PiMinusPiMinus,
        hbt::PrimitiveChannelId::ProtonProton,
        hbt::PrimitiveChannelId::ProtonBarProtonBar,
        hbt::PrimitiveChannelId::PiPlusProton,
        hbt::PrimitiveChannelId::PiMinusProtonBar,
        hbt::PrimitiveChannelId::KPlusProton,
        hbt::PrimitiveChannelId::KMinusProtonBar,
        hbt::PrimitiveChannelId::KPlusKPlus,
        hbt::PrimitiveChannelId::KMinusKMinus
    };

    if (startup.required_primitive_channels != expected_channels) {
        std::cerr
            << "event_preparation_real_fixture_test: required primitive "
            << "channels differ from the frozen fixture config.\n";
        success = false;
    }

    const std::vector<hbt::SpeciesId>& required_species =
        startup.required_species;

    success = expect_count(
        "required_species.size",
        required_species.size(),
        6U
    ) && success;

    const std::vector<hbt::SpeciesId> expected_species{
        hbt::SpeciesId::PiPlus,
        hbt::SpeciesId::PiMinus,
        hbt::SpeciesId::Proton,
        hbt::SpeciesId::ProtonBar,
        hbt::SpeciesId::KPlus,
        hbt::SpeciesId::KMinus
    };

    for (const hbt::SpeciesId species : expected_species) {
        if (!contains_species(required_species, species)) {
            std::cerr
                << "event_preparation_real_fixture_test: required species "
                << hbt::species_metadata(species).ascii_token
                << " is absent.\n";
            success = false;
        }
    }

    return success;
}

/**
 * @brief Find one accepted-particle species counter in a completed summary.
 * @param summary Completed HBT event-preparation diagnostics.
 * @param species Canonical species whose counters are requested.
 * @return Pointer to the counters, or nullptr when the species is absent.
 */
const app::HBTSpeciesPreparationCounts* find_species_counts(
    const app::HBTEventPreparationSummary& summary,
    hbt::SpeciesId species
) {
    for (const app::HBTSpeciesPreparationCounts& entry : summary.species) {
        if (entry.species == species) {
            return &entry;
        }
    }

    return nullptr;
}

/**
 * @brief Verify all frozen accepted-particle counts for one species.
 * @param summary Completed HBT event-preparation diagnostics.
 * @param species Canonical species being checked.
 * @param accepted Expected total accepted-particle count.
 * @param primordial Expected primordial membership count.
 * @param middle Expected primordial+rescattering membership count.
 * @param widest Expected widest nested-origin membership count.
 * @return `true` when all four species counts match exactly.
 */
bool verify_species_counts(
    const app::HBTEventPreparationSummary& summary,
    hbt::SpeciesId species,
    std::size_t accepted,
    std::size_t primordial,
    std::size_t middle,
    std::size_t widest
) {
    const app::HBTSpeciesPreparationCounts* counts =
        find_species_counts(summary, species);

    if (counts == nullptr) {
        std::cerr
            << "event_preparation_real_fixture_test: missing counters for "
            << hbt::species_metadata(species).ascii_token
            << ".\n";
        return false;
    }

    const std::string token{
        hbt::species_metadata(species).ascii_token
    };
    bool success = true;

    success = expect_count(
        token + ".accepted",
        counts->accepted,
        accepted
    ) && success;
    success = expect_count(
        token + ".primordial",
        counts->primordial,
        primordial
    ) && success;
    success = expect_count(
        token + ".primordial_rescattering",
        counts->primordial_rescattering,
        middle
    ) && success;
    success = expect_count(
        token + ".primordial_rescattering_decay",
        counts->primordial_rescattering_decay,
        widest
    ) && success;

    return success;
}

/**
 * @brief Verify frozen aggregate and per-species real-fixture counts.
 * @param summary Completed HBT event-preparation diagnostics.
 * @return `true` when every integer regression value matches exactly.
 */
bool verify_preparation_summary(
    const app::HBTEventPreparationSummary& summary
) {
    bool success = true;

    success = expect_count(
        "outer_events_processed",
        summary.outer_events_processed,
        1U
    ) && success;
    success = expect_count(
        "subevents_processed",
        summary.subevents_processed,
        2000U
    ) && success;
    success = expect_count(
        "raw_particles",
        summary.raw_particles,
        791392U
    ) && success;
    success = expect_count(
        "unsupported_species",
        summary.unsupported_species,
        41459U
    ) && success;
    success = expect_count(
        "unrequired_species",
        summary.unrequired_species,
        292628U
    ) && success;
    success = expect_count(
        "particle_acceptance_rejections",
        summary.particle_acceptance_rejections,
        391285U
    ) && success;
    success = expect_count(
        "numerical_rejections",
        summary.numerical_rejections.size(),
        0U
    ) && success;
    success = expect_count(
        "origin_rejections",
        summary.origin_rejections,
        0U
    ) && success;
    success = expect_count(
        "accepted_particles",
        summary.accepted_particles,
        66020U
    ) && success;

    success = expect_count(
        "emission_points.sampler",
        summary.emission_points.sampler,
        11759U
    ) && success;
    success = expect_count(
        "emission_points.propagation",
        summary.emission_points.propagation,
        54261U
    ) && success;
    success = expect_count(
        "emission_points.afterburner",
        summary.emission_points.afterburner,
        0U
    ) && success;

    success = expect_count(
        "species.size",
        summary.species.size(),
        6U
    ) && success;

    success = verify_species_counts(
        summary,
        hbt::SpeciesId::PiPlus,
        26728U,
        5497U,
        5890U,
        26728U
    ) && success;
    success = verify_species_counts(
        summary,
        hbt::SpeciesId::PiMinus,
        29709U,
        6383U,
        6842U,
        29709U
    ) && success;
    success = verify_species_counts(
        summary,
        hbt::SpeciesId::Proton,
        1671U,
        248U,
        278U,
        1671U
    ) && success;
    success = verify_species_counts(
        summary,
        hbt::SpeciesId::ProtonBar,
        1108U,
        179U,
        201U,
        1108U
    ) && success;
    success = verify_species_counts(
        summary,
        hbt::SpeciesId::KPlus,
        3481U,
        1372U,
        1465U,
        3481U
    ) && success;
    success = verify_species_counts(
        summary,
        hbt::SpeciesId::KMinus,
        3323U,
        1342U,
        1403U,
        3323U
    ) && success;

    success = expect_count(
        "subevents.size",
        summary.subevents.size(),
        2000U
    ) && success;

    if (summary.subevents.empty()) {
        return false;
    }

    std::size_t minimum = std::numeric_limits<std::size_t>::max();
    std::size_t maximum = 0U;
    std::size_t sum = 0U;

    for (const app::HBTSubeventPreparationSummary& subevent :
         summary.subevents) {
        if (subevent.outer_event_number != 1U) {
            std::cerr
                << "event_preparation_real_fixture_test: unexpected outer "
                << "event number in subevent summary.\n";
            success = false;
        }

        if (subevent.accepted_particles < minimum) {
            minimum = subevent.accepted_particles;
        }

        if (subevent.accepted_particles > maximum) {
            maximum = subevent.accepted_particles;
        }

        sum += subevent.accepted_particles;
    }

    success = expect_count(
        "accepted_per_subevent.minimum",
        minimum,
        9U
    ) && success;
    success = expect_count(
        "accepted_per_subevent.maximum",
        maximum,
        66U
    ) && success;
    success = expect_count(
        "accepted_per_subevent.sum",
        sum,
        66020U
    ) && success;

    return success;
}

/**
 * @brief Frozen pair-count expectation for one primitive channel.
 */
struct ExpectedPairCount {
    const char* label;                  ///< Canonical channel label.
    hbt::PrimitiveChannelId channel;    ///< Expected primitive channel.
    std::uint64_t routed_P;             ///< Frozen legacy primordial count.
    std::uint64_t routed_PR;            ///< Frozen legacy P+R count.
    std::uint64_t routed_PRD;           ///< Frozen legacy widest count.
};

/**
 * @brief Verify frozen formed, valid, rejected, and origin-routed pair counts.
 * @param result Completed application run containing pair-processing results.
 * @return `true` when all frozen real-fixture pair counts match exactly.
 */
bool verify_pair_processing_summary(const app::AnalysisRunSummary& result) {
    if (!result.hbt_pair_processing.has_value()) {
        std::cerr
            << "event_preparation_real_fixture_test: missing HBT "
            << "pair-processing summary.\n";
        return false;
    }

    const hbt::HBTPairProcessingSummary& summary =
        result.hbt_pair_processing.value();

    const std::vector<ExpectedPairCount> expected{
        {
            "pi_plus_pi_plus",
            hbt::PrimitiveChannelId::PiPlusPiPlus,
            7713U,
            8754U,
            178440U
        },
        {
            "pi_minus_pi_minus",
            hbt::PrimitiveChannelId::PiMinusPiMinus,
            10241U,
            11774U,
            221678U
        },
        {
            "p_p",
            hbt::PrimitiveChannelId::ProtonProton,
            17U,
            20U,
            686U
        },
        {
            "p_bar_p_bar",
            hbt::PrimitiveChannelId::ProtonBarProtonBar,
            6U,
            8U,
            335U
        },
        {
            "pi_plus_p",
            hbt::PrimitiveChannelId::PiPlusProton,
            693U,
            828U,
            22503U
        },
        {
            "pi_minus_p_bar",
            hbt::PrimitiveChannelId::PiMinusProtonBar,
            559U,
            682U,
            16654U
        },
        {
            "k_plus_p",
            hbt::PrimitiveChannelId::KPlusProton,
            164U,
            188U,
            2864U
        },
        {
            "k_minus_p_bar",
            hbt::PrimitiveChannelId::KMinusProtonBar,
            125U,
            143U,
            1852U
        },
        {
            "k_plus_k_plus",
            hbt::PrimitiveChannelId::KPlusKPlus,
            495U,
            569U,
            3070U
        },
        {
            "k_minus_k_minus",
            hbt::PrimitiveChannelId::KMinusKMinus,
            459U,
            494U,
            2759U
        }
    };

    bool success = expect_count(
        "pair_processing.total_pair_counts.channels.size",
        summary.total_pair_counts.channels.size(),
        expected.size()
    );
    success = expect_count(
        "pair_processing.total_valid_pair_counts.channels.size",
        summary.total_valid_pair_counts.channels.size(),
        expected.size()
    ) && success;
    success = expect_count(
        "pair_processing.total_numerical_rejection_counts.channels.size",
        summary.total_numerical_rejection_counts.channels.size(),
        expected.size()
    ) && success;
    success = expect_count(
        "pair_processing.numerical_rejections.size",
        summary.numerical_rejections.size(),
        0U
    ) && success;
    if (summary.total_origin_route_counts.origin_mode !=
        hbt::OriginMode::All) {
        std::cerr
            << "event_preparation_real_fixture_test: pair origin route "
            << "mode is not All.\n";
        success = false;
    }
    if (!summary.total_pair_slice_counts.entries.empty()) {
        std::cerr
            << "event_preparation_real_fixture_test: disabled slicing "
            << "created production slice entries.\n";
        success = false;
    }
    success = expect_count(
        "pair_processing.routed_P.channels.size",
        summary.total_origin_route_counts.routed_P.channels.size(),
        expected.size()
    ) && success;
    success = expect_count(
        "pair_processing.routed_PR.channels.size",
        summary.total_origin_route_counts.routed_PR.channels.size(),
        expected.size()
    ) && success;
    success = expect_count(
        "pair_processing.routed_PRD.channels.size",
        summary.total_origin_route_counts.routed_PRD.channels.size(),
        expected.size()
    ) && success;

    if (summary.total_pair_counts.channels.size() != expected.size() ||
        summary.total_valid_pair_counts.channels.size() != expected.size() ||
        summary.total_numerical_rejection_counts.channels.size() !=
            expected.size() ||
        summary.total_origin_route_counts.routed_P.channels.size() !=
            expected.size() ||
        summary.total_origin_route_counts.routed_PR.channels.size() !=
            expected.size() ||
        summary.total_origin_route_counts.routed_PRD.channels.size() !=
            expected.size()) {
        return false;
    }

    for (std::size_t i = 0U; i < expected.size(); ++i) {
        const hbt::PairChannelCount& formed =
            summary.total_pair_counts.channels[i];
        const hbt::PairChannelCount& valid =
            summary.total_valid_pair_counts.channels[i];
        const hbt::PairChannelCount& rejected =
            summary.total_numerical_rejection_counts.channels[i];
        const hbt::PairChannelCount& routed_P =
            summary.total_origin_route_counts.routed_P.channels[i];
        const hbt::PairChannelCount& routed_PR =
            summary.total_origin_route_counts.routed_PR.channels[i];
        const hbt::PairChannelCount& routed_PRD =
            summary.total_origin_route_counts.routed_PRD.channels[i];
        const ExpectedPairCount& reference = expected[i];

        if (formed.channel != reference.channel ||
            valid.channel != reference.channel ||
            rejected.channel != reference.channel ||
            routed_P.channel != reference.channel ||
            routed_PR.channel != reference.channel ||
            routed_PRD.channel != reference.channel) {
            std::cerr
                << "event_preparation_real_fixture_test: pair channel "
                << reference.label
                << " is not in the expected position.\n";
            success = false;
        }

        success = expect_pair_count(
            std::string("pair_count.") + reference.label,
            formed.pair_count,
            reference.routed_PRD
        ) && success;
        success = expect_pair_count(
            std::string("valid_pair_count.") + reference.label,
            valid.pair_count,
            reference.routed_PRD
        ) && success;
        success = expect_pair_count(
            std::string("numerical_pair_rejection_count.") +
                reference.label,
            rejected.pair_count,
            0U
        ) && success;
        success = expect_pair_count(
            std::string("routed_P.") + reference.label,
            routed_P.pair_count,
            reference.routed_P
        ) && success;
        success = expect_pair_count(
            std::string("routed_PR.") + reference.label,
            routed_PR.pair_count,
            reference.routed_PR
        ) && success;
        success = expect_pair_count(
            std::string("routed_PRD.") + reference.label,
            routed_PRD.pair_count,
            reference.routed_PRD
        ) && success;
    }

    return success;
}

}  // namespace

/**
 * @brief Run the frozen HBT regression against one real-fixture main.yaml.
 * @param argc Number of command-line arguments including the executable name.
 * @param argv Command-line argument array; argv[1] must name main.yaml.
 * @return EXIT_SUCCESS when every frozen regression value matches exactly.
 */
int main(int argc, char* argv[]) {
    if (argc != 2) {
        std::cerr
            << "usage: event_preparation_real_fixture_test <main.yaml>\n";
        return EXIT_FAILURE;
    }

    try {
        const app::AnalysisRunSummary result = app::AnalysisRunner(
            std::filesystem::path(argv[1])
        ).run();

        if (!result.hbt_event_preparation.has_value()) {
            std::cerr
                << "event_preparation_real_fixture_test: missing HBT "
                << "event-preparation summary.\n";
            return EXIT_FAILURE;
        }

        bool success = verify_startup_contract(result);
        success = verify_preparation_summary(
            result.hbt_event_preparation.value()
        ) && success;
        success = verify_pair_processing_summary(result) && success;

        return success ? EXIT_SUCCESS : EXIT_FAILURE;
    } catch (const std::exception& error) {
        std::cerr
            << "event_preparation_real_fixture_test: run failed: "
            << error.what()
            << '\n';
        return EXIT_FAILURE;
    }
}
