/**
 * @file analysis_output_writer_test.cpp
 * @brief Unit tests for the analysis-output writer.
 */

#include "output/analysis_output_writer.h"

#include <array>
#include <cstdint>
#include <cstdlib>
#include <filesystem>
#include <iomanip>
#include <iostream>
#include <limits>
#include <sstream>
#include <stdexcept>
#include <string>
#include <utility>

namespace {

/**
 * @brief Construct a complete HBT configuration for output-writer tests.
 * @param pair_slicing Validated pair-slicing configuration to preserve.
 * @param origin_mode Requested nested origin-routing mode.
 * @return Complete representative HBT configuration.
 */
hbt::HBTConfig make_hbt_config(
    const hbt::PairSlicingConfig& pair_slicing,
    hbt::OriginMode origin_mode
) {
    return hbt::HBTConfig{
        hbt::HBTSelection{{
            hbt::AnalysisProduct{{
                hbt::PrimitiveChannelId::PiPlusPiPlus
            }}
        }},
        hbt::ParticleAcceptanceConfig{
            hbt::LongitudinalVariable::Pseudorapidity,
            {0.8, 0.14, 4.0},
            {0.8, 0.4, 1.4},
            {0.8, 0.5, 4.05},
            {0.8, 1.0, 10000.0},
            {0.8, 0.3, 10000.0}
        },
        pair_slicing,
        hbt::HBTHistogramConfig{
            {10U, 0.0, 10.0, 1.0},
            {10U, 0.0, 10.0, 1.0},
            {20U, -10.0, 10.0, 1.0}
        },
        origin_mode,
        hbt::FitEstimatorMode::All
    };
}

/**
 * @brief Construct one-channel pair counts for compact output tests.
 * @param count Pair count assigned to the pi-plus identical channel.
 * @return One-entry pair-count summary.
 */
hbt::PairCountSummary make_pi_plus_counts(std::uint64_t count) {
    return hbt::PairCountSummary{{
        {hbt::PrimitiveChannelId::PiPlusPiPlus, count}
    }};
}

/**
 * @brief Construct nested origin counts for one slice and one channel.
 * @param primordial Primordial-route count.
 * @param primordial_rescattering Middle-route count.
 * @param primordial_rescattering_decay Widest-route count.
 * @return One-channel nested origin counts.
 */
hbt::PairSliceOriginCounts make_pi_plus_slice_origin_counts(
    std::uint64_t primordial,
    std::uint64_t primordial_rescattering,
    std::uint64_t primordial_rescattering_decay
) {
    return hbt::PairSliceOriginCounts{
        make_pi_plus_counts(primordial),
        make_pi_plus_counts(primordial_rescattering),
        make_pi_plus_counts(primordial_rescattering_decay)
    };
}

/**
 * @brief Verify complete rejected-particle serialization.
 * @return true when identifiers, reason, momentum, and diagnostic are present.
 */
bool verify_rejected_particle_serialization() {
    hbt::RejectedParticleReport report;
    report.add({
        184U,
        731,
        628,
        -211,
        -1,
        hbt::SpeciesId::PiMinus,
        {2519.0, 0.12, 0.30, 2519.0},
        {10.0, 1.0, 2.0, 3.0},
        0.138,
        0,
        4.0,
        hbt::ParticleRejectionReason::NonPositiveInvariantMassSquared,
        -0.003305,
        std::nullopt
    });

    std::ostringstream output;
    output::write_rejected_particle_report(report, output);
    const std::string text = output.str();

    const char* required[] = {
        "numerical_particle_rejections: 1",
        "outer_event_number: 184",
        "subevent_id: 731",
        "particle_id: 628",
        "pdg: -211",
        "charge: -1",
        "species: pi_minus",
        "raw_mass_gev: 0.138",
        "ncoll: 0",
        "time_last_coll_fm: 4",
        "momentum_gev:",
        "raw_position_fm:",
        "reason: non_positive_invariant_mass_squared",
        "diagnostic_value: -0.003305"
    };

    for (const char* token : required) {
        if (text.find(token) == std::string::npos) {
            std::cerr
                << "analysis_output_writer_test: missing token: "
                << token << ".\n";
            return false;
        }
    }

    return true;
}

/**
 * @brief Verify non-positive energy uses its stable output token.
 * @return true when the new particle-rejection reason is serialized exactly.
 */
bool verify_nonpositive_energy_serialization() {
    hbt::RejectedParticleReport report;
    report.add({
        1U,
        2,
        3,
        211,
        1,
        hbt::SpeciesId::PiPlus,
        {-2.0, 0.5, 0.0, 0.0},
        {10.0, 1.0, 2.0, 3.0},
        0.138,
        0,
        4.0,
        hbt::ParticleRejectionReason::NonPositiveEnergy,
        -2.0,
        std::nullopt
    });

    std::ostringstream output;
    output::write_rejected_particle_report(report, output);
    const std::string text = output.str();

    return text.find("reason: non_positive_energy") != std::string::npos &&
           text.find("diagnostic_value: -2") != std::string::npos;
}

/**
 * @brief Verify emission-position diagnostics preserve their full vector.
 * @return true when the source-specific reason and position are serialized.
 */
bool verify_emission_position_serialization() {
    hbt::RejectedParticleReport report;
    report.add({
        2U,
        19,
        77,
        211,
        1,
        hbt::SpeciesId::PiPlus,
        {2.0, -1.0, 0.0, 0.0},
        {10.0, 5.0, 7.0, 9.0},
        0.138,
        1,
        0.0,
        hbt::ParticleRejectionReason::
            NonFinitePropagationEmissionPosition,
        std::nullopt,
        common::FourVector{0.0, 1.0, 2.0, 3.0}
    });

    std::ostringstream output;
    output::write_rejected_particle_report(report, output);
    const std::string text = output.str();

    return
        text.find(
            "reason: non_finite_propagation_emission_position") !=
            std::string::npos &&
        text.find("diagnostic_value: none") != std::string::npos &&
        text.find("diagnostic_position_fm: [0, 1, 2, 3]") !=
            std::string::npos;
}

/**
 * @brief Verify complete rejected-pair serialization.
 * @return true when pair identity, particles, kinematics, and reason are
 *         present.
 */
bool verify_rejected_pair_serialization() {
    hbt::RejectedPairReport report;
    report.add({
        3U,
        17,
        hbt::PrimitiveChannelId::PiPlusPiPlus,
        42U,
        {
            hbt::SpeciesId::PiPlus,
            {2.0, 0.5, 0.1, 0.0},
            0.139,
            211,
            1
        },
        {
            hbt::SpeciesId::PiPlus,
            {2.0, 0.7, -0.1, 0.0},
            0.139,
            211,
            1
        },
        {
            {std::numeric_limits<double>::infinity(), 0.0, 0.0, 0.0},
            0.0,
            0.0,
            std::numeric_limits<double>::infinity(),
            1.0
        },
        hbt::PairRejectionReason::NonFiniteKt
    });

    std::ostringstream output;
    output::write_rejected_pair_report(report, output);
    const std::string text = output.str();

    const char* required[] = {
        "numerical_pair_rejections: 1",
        "non_finite_kt: 1",
        "non_finite_mt: 0",
        "outer_event_number: 3",
        "subevent_id: 17",
        "channel: pi_plus_pi_plus",
        "pair_ordinal_in_channel: 42",
        "particle_a_species: pi_plus",
        "particle_a_pdg: 211",
        "particle_b_species: pi_plus",
        "kt_gev: inf",
        "mt_gev: 1",
        "reason: non_finite_kt"
    };

    for (const char* token : required) {
        if (text.find(token) == std::string::npos) {
            std::cerr
                << "analysis_output_writer_test: missing pair token: "
                << token << ".\n";
            return false;
        }
    }

    return true;
}

/**
 * @brief Verify every pair-rejection reason has its exact stable output token.
 * @return true when all eleven count and record tokens are serialized.
 */
bool verify_all_pair_rejection_reason_tokens() {
    using Reason = hbt::PairRejectionReason;
    const std::array<Reason, 11U> reasons{
        Reason::NonFiniteKt,
        Reason::NonFiniteMt,
        Reason::NonFiniteDeltaTLab,
        Reason::NonFiniteDeltaTLcms,
        Reason::NonFiniteDeltaTPrf,
        Reason::NonFiniteROutLcms,
        Reason::NonFiniteROutPrf,
        Reason::NonFiniteRSide,
        Reason::NonFiniteRLong,
        Reason::NonFiniteRRadialLcms,
        Reason::NonFiniteRRadialPrf
    };
    const std::array<const char*, 11U> tokens{
        "non_finite_kt",
        "non_finite_mt",
        "non_finite_delta_t_lab",
        "non_finite_delta_t_lcms",
        "non_finite_delta_t_prf",
        "non_finite_r_out_lcms",
        "non_finite_r_out_prf",
        "non_finite_r_side",
        "non_finite_r_long",
        "non_finite_r_radial_lcms",
        "non_finite_r_radial_prf"
    };
    hbt::RejectedPairReport report;

    for (std::size_t index = 0U; index < reasons.size(); ++index) {
        report.add({
            3U,
            17,
            hbt::PrimitiveChannelId::PiPlusPiPlus,
            static_cast<std::uint64_t>(index + 1U),
            {hbt::SpeciesId::PiPlus, {2.0, 0.5, 0.1, 0.0},
             0.139, 211, 1},
            {hbt::SpeciesId::PiPlus, {2.0, 0.7, -0.1, 0.0},
             0.139, 211, 1},
            {{4.0, 1.2, 0.0, 0.0}, 0.6, 0.0, 0.6, 1.0},
            reasons[index]
        });
    }

    std::ostringstream output;
    output::write_rejected_pair_report(report, output);
    const std::string text = output.str();

    for (const char* token : tokens) {
        const std::string count_token =
            std::string("  ") + token + ": 1\n";
        const std::string record_token =
            std::string("    reason: ") + token + "\n";
        if (text.find(count_token) == std::string::npos ||
            text.find(record_token) == std::string::npos) {
            std::cerr
                << "analysis_output_writer_test: missing pair reason token: "
                << token << ".\n";
            return false;
        }
    }
    return true;
}

/**
 * @brief Verify the run writer handles disabled HBT without fake HBT output.
 * @return true when only the disabled run state is reported.
 */
bool verify_disabled_hbt_serialization() {
    app::AnalysisRunSummary result{
        app::AnalysisStartupState{
            config::RunConfig{
                "/tmp/events",
                "/tmp/output",
                1U,
                1U,
                false,
                std::nullopt
            },
            std::nullopt,
            std::nullopt
        },
        std::nullopt,
        std::nullopt
    };

    std::ostringstream output;
    output::write_analysis_output(result, output);
    const std::string text = output.str();

    return text.find("hbt_enabled: false") != std::string::npos &&
           text.find("raw_particles:") == std::string::npos &&
           text.find("pair_counts_by_primitive_channel:") ==
               std::string::npos;
}

/**
 * @brief Verify run-total pair counts use canonical channel names.
 * @return true when ordered primitive-channel counts are serialized exactly.
 */
bool verify_pair_count_serialization() {
    const hbt::PairSlicingConfig pair_slicing{
        {false, {}},
        {false, {}}
    };

    app::AnalysisRunSummary result{
        app::AnalysisStartupState{
            config::RunConfig{
                "/tmp/events",
                "/tmp/output",
                1U,
                1U,
                true,
                std::filesystem::path("/tmp/hbt.yaml")
            },
            make_hbt_config(pair_slicing, hbt::OriginMode::All),
            std::nullopt
        },
        std::nullopt,
        hbt::HBTPairProcessingSummary{
            hbt::PairCountSummary{{
                {hbt::PrimitiveChannelId::PiPlusPiPlus, 12U},
                {hbt::PrimitiveChannelId::KMinusProton, 34U}
            }},
            hbt::PairCountSummary{{
                {hbt::PrimitiveChannelId::PiPlusPiPlus, 12U},
                {hbt::PrimitiveChannelId::KMinusProton, 34U}
            }},
            hbt::PairCountSummary{{
                {hbt::PrimitiveChannelId::PiPlusPiPlus, 0U},
                {hbt::PrimitiveChannelId::KMinusProton, 0U}
            }},
            {
                hbt::OriginMode::All,
                hbt::PairCountSummary{{
                    {hbt::PrimitiveChannelId::PiPlusPiPlus, 5U},
                    {hbt::PrimitiveChannelId::KMinusProton, 10U}
                }},
                hbt::PairCountSummary{{
                    {hbt::PrimitiveChannelId::PiPlusPiPlus, 8U},
                    {hbt::PrimitiveChannelId::KMinusProton, 20U}
                }},
                hbt::PairCountSummary{{
                    {hbt::PrimitiveChannelId::PiPlusPiPlus, 12U},
                    {hbt::PrimitiveChannelId::KMinusProton, 34U}
                }}
            },
            {pair_slicing, hbt::OriginMode::All, {}},
            {},
            {}
        }
    };

    std::ostringstream output;
    output::write_analysis_output(result, output);
    const std::string text = output.str();

    return
        text.find("pair_counts_by_primitive_channel:\n") !=
            std::string::npos &&
        text.find("  pi_plus_pi_plus: 12\n") != std::string::npos &&
        text.find("  k_minus_p: 34\n") != std::string::npos &&
        text.find("valid_pair_counts_by_primitive_channel:\n") !=
            std::string::npos &&
        text.find("  pi_plus_pi_plus: 12\n") != std::string::npos &&
        text.find("numerical_pair_rejections_by_primitive_channel:\n") !=
            std::string::npos &&
        text.find("  pi_plus_pi_plus: 0\n") != std::string::npos &&
        text.find("pair_origin_route_mode: all\n") != std::string::npos &&
        text.find("routed_P_pair_counts_by_primitive_channel:\n"
                  "  pi_plus_pi_plus: 5\n") != std::string::npos &&
        text.find("routed_PR_pair_counts_by_primitive_channel:\n"
                  "  pi_plus_pi_plus: 8\n") != std::string::npos &&
        text.find("routed_PRD_pair_counts_by_primitive_channel:\n"
                  "  pi_plus_pi_plus: 12\n") != std::string::npos &&
        text.find("pair_slice_kt_enabled: false\n") !=
            std::string::npos &&
        text.find("pair_slice_mt_enabled: false\n") !=
            std::string::npos &&
        text.find("pair_slice_counts:\n") != std::string::npos &&
        text.find("  - kt_slice_index:") == std::string::npos &&
        text.find("numerical_pair_rejections: 0\n") != std::string::npos;
}

/**
 * @brief Verify production bins and Cartesian slice counts are serialized.
 * @return true when startup bins and every nested slice route are present.
 */
bool verify_pair_slice_serialization() {
    const hbt::PairSlicingConfig pair_slicing{
        {true, {0.125, 0.25, 0.5}},
        {true, {0.5, 0.75, 1.0}}
    };

    const hbt::PairSliceCountSummary slice_counts{
        pair_slicing,
        hbt::OriginMode::All,
        {
            {0U, 0U, make_pi_plus_slice_origin_counts(1U, 2U, 3U)},
            {0U, 1U, make_pi_plus_slice_origin_counts(4U, 5U, 6U)},
            {1U, 0U, make_pi_plus_slice_origin_counts(7U, 8U, 9U)},
            {1U, 1U, make_pi_plus_slice_origin_counts(10U, 11U, 12U)}
        }
    };

    app::AnalysisRunSummary result{
        {
            {
                "/tmp/events",
                "/tmp/output",
                1U,
                1U,
                true,
                std::filesystem::path("/tmp/hbt.yaml")
            },
            make_hbt_config(pair_slicing, hbt::OriginMode::All),
            std::nullopt
        },
        std::nullopt,
        hbt::HBTPairProcessingSummary{
            make_pi_plus_counts(30U),
            make_pi_plus_counts(30U),
            make_pi_plus_counts(0U),
            {
                hbt::OriginMode::All,
                make_pi_plus_counts(22U),
                make_pi_plus_counts(26U),
                make_pi_plus_counts(30U)
            },
            slice_counts,
            {},
            {}
        }
    };

    std::ostringstream output;
    output::write_analysis_output(result, output);
    const std::string text = output.str();

    const char* required[] = {
        "pair_slice_kt_enabled: true",
        "pair_slice_kt_bin_edges_gev: [0.125, 0.25, 0.5]",
        "pair_slice_mt_enabled: true",
        "pair_slice_mt_bin_edges_gev: [0.5, 0.75, 1]",
        "  - kt_slice_index: 0\n    mt_slice_index: 0",
        "  - kt_slice_index: 0\n    mt_slice_index: 1",
        "  - kt_slice_index: 1\n    mt_slice_index: 0",
        "  - kt_slice_index: 1\n    mt_slice_index: 1",
        "    routed_P_pair_counts_by_primitive_channel:\n"
        "      pi_plus_pi_plus: 1",
        "    routed_PR_pair_counts_by_primitive_channel:\n"
        "      pi_plus_pi_plus: 11",
        "    routed_PRD_pair_counts_by_primitive_channel:\n"
        "      pi_plus_pi_plus: 12"
    };

    for (const char* token : required) {
        if (text.find(token) == std::string::npos) {
            std::cerr
                << "analysis_output_writer_test: missing slice token: "
                << token << ".\n";
            return false;
        }
    }

    return true;
}

/**
 * @brief Verify kT-only slicing serializes a missing mT slice index.
 * @return true when configured run bins and `mt_slice_index: none` are present.
 */
bool verify_kt_only_pair_slice_serialization() {
    const hbt::PairSlicingConfig pair_slicing{
        {true, {0.125, 0.25}},
        {false, {0.5, 1.0}}
    };

    const hbt::PairSliceCountSummary slice_counts{
        pair_slicing,
        hbt::OriginMode::All,
        {
            {
                0U,
                std::nullopt,
                make_pi_plus_slice_origin_counts(1U, 2U, 3U)
            }
        }
    };

    app::AnalysisRunSummary result{
        {
            {
                "/tmp/events",
                "/tmp/output",
                1U,
                1U,
                true,
                std::filesystem::path("/tmp/hbt.yaml")
            },
            make_hbt_config(pair_slicing, hbt::OriginMode::All),
            std::nullopt
        },
        std::nullopt,
        hbt::HBTPairProcessingSummary{
            make_pi_plus_counts(3U),
            make_pi_plus_counts(3U),
            make_pi_plus_counts(0U),
            {
                hbt::OriginMode::All,
                make_pi_plus_counts(1U),
                make_pi_plus_counts(2U),
                make_pi_plus_counts(3U)
            },
            slice_counts,
            {},
            {}
        }
    };

    std::ostringstream output;
    output::write_analysis_output(result, output);
    const std::string text = output.str();

    return
        text.find("pair_slice_kt_enabled: true\n") !=
            std::string::npos &&
        text.find("pair_slice_kt_bin_edges_gev: [0.125, 0.25]\n") !=
            std::string::npos &&
        text.find("pair_slice_mt_enabled: false\n") !=
            std::string::npos &&
        text.find("pair_slice_mt_bin_edges_gev: [0.5, 1]\n") !=
            std::string::npos &&
        text.find(
            "  - kt_slice_index: 0\n"
            "    mt_slice_index: none\n"
        ) != std::string::npos;
}

/**
 * @brief Verify mT-only slicing serializes a missing kT slice index.
 * @return true when configured run bins and `kt_slice_index: none` are present.
 */
bool verify_mt_only_pair_slice_serialization() {
    const hbt::PairSlicingConfig pair_slicing{
        {false, {0.125, 0.25}},
        {true, {0.5, 1.0}}
    };

    const hbt::PairSliceCountSummary slice_counts{
        pair_slicing,
        hbt::OriginMode::All,
        {
            {
                std::nullopt,
                0U,
                make_pi_plus_slice_origin_counts(1U, 2U, 3U)
            }
        }
    };

    app::AnalysisRunSummary result{
        {
            {
                "/tmp/events",
                "/tmp/output",
                1U,
                1U,
                true,
                std::filesystem::path("/tmp/hbt.yaml")
            },
            make_hbt_config(pair_slicing, hbt::OriginMode::All),
            std::nullopt
        },
        std::nullopt,
        hbt::HBTPairProcessingSummary{
            make_pi_plus_counts(3U),
            make_pi_plus_counts(3U),
            make_pi_plus_counts(0U),
            {
                hbt::OriginMode::All,
                make_pi_plus_counts(1U),
                make_pi_plus_counts(2U),
                make_pi_plus_counts(3U)
            },
            slice_counts,
            {},
            {}
        }
    };

    std::ostringstream output;
    output::write_analysis_output(result, output);
    const std::string text = output.str();

    return
        text.find("pair_slice_kt_enabled: false\n") !=
            std::string::npos &&
        text.find("pair_slice_kt_bin_edges_gev: [0.125, 0.25]\n") !=
            std::string::npos &&
        text.find("pair_slice_mt_enabled: true\n") !=
            std::string::npos &&
        text.find("pair_slice_mt_bin_edges_gev: [0.5, 1]\n") !=
            std::string::npos &&
        text.find(
            "  - kt_slice_index: none\n"
            "    mt_slice_index: 0\n"
        ) != std::string::npos;
}

/**
 * @brief Verify disabled axes preserve configured edges without dummy slices.
 * @return true when retained edges are serialized and entries remain absent.
 */
bool verify_disabled_pair_slice_serialization() {
    const hbt::PairSlicingConfig pair_slicing{
        {false, {0.125, 0.25}},
        {false, {0.5, 1.0}}
    };

    app::AnalysisRunSummary result{
        {
            {
                "/tmp/events",
                "/tmp/output",
                1U,
                1U,
                true,
                std::filesystem::path("/tmp/hbt.yaml")
            },
            make_hbt_config(pair_slicing, hbt::OriginMode::All),
            std::nullopt
        },
        std::nullopt,
        hbt::HBTPairProcessingSummary{
            make_pi_plus_counts(0U),
            make_pi_plus_counts(0U),
            make_pi_plus_counts(0U),
            {
                hbt::OriginMode::All,
                make_pi_plus_counts(0U),
                make_pi_plus_counts(0U),
                make_pi_plus_counts(0U)
            },
            {pair_slicing, hbt::OriginMode::All, {}},
            {},
            {}
        }
    };

    std::ostringstream output;
    output::write_analysis_output(result, output);
    const std::string text = output.str();

    return
        text.find("pair_slice_kt_enabled: false\n") !=
            std::string::npos &&
        text.find("pair_slice_kt_bin_edges_gev: [0.125, 0.25]\n") !=
            std::string::npos &&
        text.find("pair_slice_mt_enabled: false\n") !=
            std::string::npos &&
        text.find("pair_slice_mt_bin_edges_gev: [0.5, 1]\n") !=
            std::string::npos &&
        text.find("pair_slice_counts:\n") != std::string::npos &&
        text.find("  - kt_slice_index:") == std::string::npos;
}

/**
 * @brief Verify startup bins cannot disagree with serialized slice counts.
 * @return true when a slicing-layout mismatch is rejected structurally.
 */
bool verify_pair_slice_config_mismatch_rejected() {
    const hbt::PairSlicingConfig startup_slicing{
        {true, {0.125, 0.25}},
        {false, {}}
    };
    const hbt::PairSlicingConfig summary_slicing{
        {true, {0.125, 0.5}},
        {false, {}}
    };

    app::AnalysisRunSummary result{
        {
            {
                "/tmp/events",
                "/tmp/output",
                1U,
                1U,
                true,
                std::filesystem::path("/tmp/hbt.yaml")
            },
            make_hbt_config(startup_slicing, hbt::OriginMode::All),
            std::nullopt
        },
        std::nullopt,
        hbt::HBTPairProcessingSummary{
            {},
            {},
            {},
            {hbt::OriginMode::All, {}, {}, {}},
            {summary_slicing, hbt::OriginMode::All, {}},
            {},
            {}
        }
    };

    std::ostringstream output;
    try {
        output::write_analysis_output(result, output);
    } catch (const std::logic_error&) {
        return true;
    }

    std::cerr
        << "analysis_output_writer_test: expected slicing mismatch error.\n";
    return false;
}

/**
 * @brief Verify slice counts cannot disagree with the pair origin mode.
 * @return true when an origin-routing mismatch is rejected structurally.
 */
bool verify_pair_slice_origin_mode_mismatch_rejected() {
    const hbt::PairSlicingConfig pair_slicing{
        {false, {}},
        {false, {}}
    };

    app::AnalysisRunSummary result{
        {
            {
                "/tmp/events",
                "/tmp/output",
                1U,
                1U,
                true,
                std::filesystem::path("/tmp/hbt.yaml")
            },
            make_hbt_config(pair_slicing, hbt::OriginMode::All),
            std::nullopt
        },
        std::nullopt,
        hbt::HBTPairProcessingSummary{
            {},
            {},
            {},
            {hbt::OriginMode::All, {}, {}, {}},
            {pair_slicing, hbt::OriginMode::Primordial, {}},
            {},
            {}
        }
    };

    std::ostringstream output;
    try {
        output::write_analysis_output(result, output);
    } catch (const std::logic_error&) {
        return true;
    }

    std::cerr
        << "analysis_output_writer_test: expected origin mismatch error.\n";
    return false;
}

/**
 * @brief Verify histogram range diagnostics preserve exact configured ranges.
 * @return true when affected identity, counts, and 17-digit range are present.
 */
bool verify_histogram_range_warning_serialization() {
    const hbt::PairSlicingConfig pair_slicing{
        {false, {}},
        {false, {}}
    };
    hbt::HBTConfig config =
        make_hbt_config(pair_slicing, hbt::OriginMode::All);
    constexpr double delta_t_min = -9.876543210987654;
    constexpr double delta_t_max = 9.123456789012345;
    config.histogram_config.delta_t = {
        20U,
        delta_t_min,
        delta_t_max,
        20.0 / (delta_t_max - delta_t_min)
    };
    hbt::RawHistogramState histograms =
        hbt::make_zero_raw_histogram_state(config);
    histograms.products[0].origins[1].global.osl.overflow_counts[2] = 7U;
    histograms.products[0].origins[2]
        .global.delta_t.underflow_counts[0] = 3U;

    app::AnalysisRunSummary result{
        {
            {
                "/tmp/events",
                "/tmp/output",
                1U,
                1U,
                true,
                std::filesystem::path("/tmp/hbt.yaml")
            },
            config,
            std::nullopt
        },
        std::nullopt,
        std::nullopt,
        std::move(histograms)
    };

    std::ostringstream output;
    output::write_analysis_output(result, output);
    const std::string text = output.str();
    const char* required[] = {
        "histogram_range_warnings:",
        "product_index: 0",
        "product: pi_plus_pi_plus",
        "origin: PR",
        "destination: global",
        "observable: r_side",
        "underflow: 0",
        "overflow: 7",
        "origin: PRD",
        "observable: delta_t_lab",
        "underflow: 3"
    };

    for (const char* token : required) {
        if (text.find(token) == std::string::npos) {
            std::cerr
                << "analysis_output_writer_test: missing histogram warning "
                << token << ".\n";
            return false;
        }
    }

    std::ostringstream expected_range;
    expected_range << std::setprecision(17)
                   << "range: [" << delta_t_min << ", "
                   << delta_t_max << ")";
    if (text.find(expected_range.str()) == std::string::npos) {
        std::cerr
            << "analysis_output_writer_test: histogram range lost precision.\n";
        return false;
    }
    return true;
}

/**
 * @brief Verify malformed raw histogram dimensions fail at output boundary.
 * @return true when product, origin, slice, and family mismatches all throw.
 */
bool verify_histogram_layout_mismatches_are_rejected() {
    const hbt::PairSlicingConfig pair_slicing{
        {true, {0.0, 1.0, 2.0}},
        {false, {}}
    };
    const hbt::HBTConfig config =
        make_hbt_config(pair_slicing, hbt::OriginMode::All);

    const auto output_rejects = [&](hbt::RawHistogramState state) {
        app::AnalysisRunSummary result{
            {
                {
                    "/tmp/events",
                    "/tmp/output",
                    1U,
                    1U,
                    true,
                    std::filesystem::path("/tmp/hbt.yaml")
                },
                config,
                std::nullopt
            },
            std::nullopt,
            std::nullopt,
            std::move(state)
        };
        std::ostringstream output;
        try {
            output::write_analysis_output(result, output);
        } catch (const std::logic_error&) {
            return true;
        }
        return false;
    };

    hbt::RawHistogramState missing_product =
        hbt::make_zero_raw_histogram_state(config);
    missing_product.products.clear();
    hbt::RawHistogramState missing_origin =
        hbt::make_zero_raw_histogram_state(config);
    missing_origin.products[0].origins.pop_back();
    hbt::RawHistogramState missing_slice =
        hbt::make_zero_raw_histogram_state(config);
    missing_slice.products[0].origins[0].slices.pop_back();
    hbt::RawHistogramState missing_family_bin =
        hbt::make_zero_raw_histogram_state(config);
    missing_family_bin.products[0].origins[0].global.osl.bins.pop_back();

    if (!output_rejects(std::move(missing_product)) ||
        !output_rejects(std::move(missing_origin)) ||
        !output_rejects(std::move(missing_slice)) ||
        !output_rejects(std::move(missing_family_bin))) {
        std::cerr
            << "analysis_output_writer_test: malformed histogram layout "
            << "was serialized.\n";
        return false;
    }
    return true;
}

/**
 * @brief Verify pair summaries require their resolved HBT configuration.
 * @return true when missing HBT startup configuration is rejected.
 */
bool verify_pair_summary_requires_hbt_config() {
    app::AnalysisRunSummary result{
        {
            {
                "/tmp/events",
                "/tmp/output",
                1U,
                1U,
                true,
                std::filesystem::path("/tmp/hbt.yaml")
            },
            std::nullopt,
            std::nullopt
        },
        std::nullopt,
        hbt::HBTPairProcessingSummary{}
    };

    std::ostringstream output;
    try {
        output::write_analysis_output(result, output);
    } catch (const std::logic_error&) {
        return true;
    }

    std::cerr
        << "analysis_output_writer_test: expected missing HBT config error.\n";
    return false;
}

}  // namespace

/**
 * @brief Run all analysis-output-writer unit tests.
 * @return EXIT_SUCCESS when every test passes, otherwise EXIT_FAILURE.
 */
int main() {
    bool success = true;
    success = verify_rejected_particle_serialization() && success;
    success = verify_nonpositive_energy_serialization() && success;
    success = verify_emission_position_serialization() && success;
    success = verify_rejected_pair_serialization() && success;
    success = verify_all_pair_rejection_reason_tokens() && success;
    success = verify_disabled_hbt_serialization() && success;
    success = verify_pair_count_serialization() && success;
    success = verify_pair_slice_serialization() && success;
    success = verify_kt_only_pair_slice_serialization() && success;
    success = verify_mt_only_pair_slice_serialization() && success;
    success = verify_disabled_pair_slice_serialization() && success;
    success = verify_pair_slice_config_mismatch_rejected() && success;
    success = verify_pair_slice_origin_mode_mismatch_rejected() && success;
    success = verify_pair_summary_requires_hbt_config() && success;
    success = verify_histogram_range_warning_serialization() && success;
    success = verify_histogram_layout_mismatches_are_rejected() && success;
    return success ? EXIT_SUCCESS : EXIT_FAILURE;
}
