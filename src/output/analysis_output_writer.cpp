/**
 * @file analysis_output_writer.cpp
 * @brief Implementation of analysis-output serialization.
 */

#include "output/analysis_output_writer.h"
#include "output/hbt_production_output.h"

#include "hbt/channels/channel_catalog.h"
#include "hbt/species/species.h"

#include <array>
#include <cstddef>
#include <iomanip>
#include <optional>
#include <ostream>
#include <stdexcept>
#include <string>
#include <string_view>

namespace output {
namespace {

/**
 * @brief Return a stable text token for one rejection reason.
 * @param reason Rejection reason to serialize.
 * @return Stable ASCII reason token.
 * @throws std::invalid_argument If reason is Count or invalid.
 */
std::string_view rejection_reason_token(
    hbt::ParticleRejectionReason reason
) {
    switch (reason) {
        case hbt::ParticleRejectionReason::NonFiniteMomentum:
            return "non_finite_momentum";
        case hbt::ParticleRejectionReason::NonPositiveEnergy:
            return "non_positive_energy";
        case hbt::ParticleRejectionReason::NonFiniteTransverseMomentum:
            return "non_finite_transverse_momentum";
        case hbt::ParticleRejectionReason::InvalidRapidityInput:
            return "invalid_rapidity_input";
        case hbt::ParticleRejectionReason::NonFiniteRapidity:
            return "non_finite_rapidity";
        case hbt::ParticleRejectionReason::InvalidPseudorapidityInput:
            return "invalid_pseudorapidity_input";
        case hbt::ParticleRejectionReason::NonFinitePseudorapidity:
            return "non_finite_pseudorapidity";
        case hbt::ParticleRejectionReason::NonFiniteInvariantMassSquared:
            return "non_finite_invariant_mass_squared";
        case hbt::ParticleRejectionReason::NonPositiveInvariantMassSquared:
            return "non_positive_invariant_mass_squared";
        case hbt::ParticleRejectionReason::NonFiniteInvariantMass:
            return "non_finite_invariant_mass";
        case hbt::ParticleRejectionReason::
                NonFiniteSamplerEmissionPosition:
            return "non_finite_sampler_emission_position";
        case hbt::ParticleRejectionReason::
                NonFinitePropagationEmissionPosition:
            return "non_finite_propagation_emission_position";
        case hbt::ParticleRejectionReason::
                NonFiniteAfterburnerEmissionPosition:
            return "non_finite_afterburner_emission_position";
        case hbt::ParticleRejectionReason::Count:
            break;
    }

    throw std::invalid_argument(
        "analysis output: invalid particle rejection reason"
    );
}

/**
 * @brief Serialize aggregate rejection counts by stable reason token.
 * @param report Rejected-particle report to inspect.
 * @param output Destination output stream.
 */
void write_rejection_counts(
    const hbt::RejectedParticleReport& report,
    std::ostream& output
) {
    const std::size_t count = static_cast<std::size_t>(
        hbt::ParticleRejectionReason::Count
    );

    for (std::size_t index = 0U; index < count; ++index) {
        const auto reason =
            static_cast<hbt::ParticleRejectionReason>(index);
        output
            << "  " << rejection_reason_token(reason) << ": "
            << report.count(reason) << '\n';
    }
}

/**
 * @brief Return a stable text token for one pair-rejection reason.
 * @param reason Pair-rejection reason to serialize.
 * @return Stable ASCII reason token.
 * @throws std::invalid_argument If reason is Count or invalid.
 */
std::string_view pair_rejection_reason_token(
    hbt::PairRejectionReason reason
) {
    switch (reason) {
        case hbt::PairRejectionReason::NonFiniteKt:
            return "non_finite_kt";
        case hbt::PairRejectionReason::NonFiniteMt:
            return "non_finite_mt";
        case hbt::PairRejectionReason::NonFiniteDeltaTLab:
            return "non_finite_delta_t_lab";
        case hbt::PairRejectionReason::NonFiniteDeltaTLcms:
            return "non_finite_delta_t_lcms";
        case hbt::PairRejectionReason::NonFiniteDeltaTPrf:
            return "non_finite_delta_t_prf";
        case hbt::PairRejectionReason::NonFiniteROutLcms:
            return "non_finite_r_out_lcms";
        case hbt::PairRejectionReason::NonFiniteROutPrf:
            return "non_finite_r_out_prf";
        case hbt::PairRejectionReason::NonFiniteRSide:
            return "non_finite_r_side";
        case hbt::PairRejectionReason::NonFiniteRLong:
            return "non_finite_r_long";
        case hbt::PairRejectionReason::NonFiniteRRadialLcms:
            return "non_finite_r_radial_lcms";
        case hbt::PairRejectionReason::NonFiniteRRadialPrf:
            return "non_finite_r_radial_prf";
        case hbt::PairRejectionReason::Count:
            break;
    }

    throw std::invalid_argument(
        "analysis output: invalid pair rejection reason"
    );
}

/**
 * @brief Serialize aggregate pair-rejection counts by stable reason token.
 * @param report Rejected-pair report to inspect.
 * @param output Destination output stream.
 */
void write_pair_rejection_counts(
    const hbt::RejectedPairReport& report,
    std::ostream& output
) {
    const std::size_t count = static_cast<std::size_t>(
        hbt::PairRejectionReason::Count
    );

    for (std::size_t index = 0U; index < count; ++index) {
        const auto reason = static_cast<hbt::PairRejectionReason>(index);
        output
            << "  " << pair_rejection_reason_token(reason) << ": "
            << report.count(reason) << '\n';
    }
}

/**
 * @brief Serialize one ordered primitive-channel pair-count summary.
 * @param label Stable section label.
 * @param counts Ordered primitive-channel counts.
 * @param output Destination output stream.
 */
void write_pair_count_summary(
    std::string_view label,
    const hbt::PairCountSummary& counts,
    std::ostream& output
) {
    output << label << ":\n";

    for (const hbt::PairChannelCount& entry : counts.channels) {
        const hbt::PrimitiveChannel& channel =
            hbt::primitive_channel_definition(entry.channel);
        output
            << "  " << channel.canonical_name << ": "
            << entry.pair_count << '\n';
    }
}

/**
 * @brief Compare two validated pair-slicing axis definitions exactly.
 * @param lhs First axis definition.
 * @param rhs Second axis definition.
 * @return true when activation and bin edges are identical.
 */
bool same_pair_slicing_axis(
    const hbt::PairSlicingAxisConfig& lhs,
    const hbt::PairSlicingAxisConfig& rhs
) {
    return lhs.enabled == rhs.enabled &&
           lhs.bin_edges_gev == rhs.bin_edges_gev;
}

/**
 * @brief Compare two validated pair-slicing configurations exactly.
 * @param lhs First slicing configuration.
 * @param rhs Second slicing configuration.
 * @return true when both axis definitions are identical.
 */
bool same_pair_slicing_config(
    const hbt::PairSlicingConfig& lhs,
    const hbt::PairSlicingConfig& rhs
) {
    return same_pair_slicing_axis(lhs.kt, rhs.kt) &&
           same_pair_slicing_axis(lhs.mt, rhs.mt);
}

/**
 * @brief Serialize one configured kinetic-axis definition.
 * @param axis_name Stable lowercase axis token.
 * @param axis Validated production axis definition.
 * @param output Destination output stream.
 */
void write_pair_slicing_axis(
    std::string_view axis_name,
    const hbt::PairSlicingAxisConfig& axis,
    std::ostream& output
) {
    output << "pair_slice_" << axis_name << "_enabled: "
           << (axis.enabled ? "true" : "false") << '\n';
    output << "pair_slice_" << axis_name << "_bin_edges_gev: [";

    for (std::size_t index = 0U;
         index < axis.bin_edges_gev.size();
         ++index) {
        if (index != 0U) {
            output << ", ";
        }
        output << axis.bin_edges_gev[index];
    }

    output << "]\n";
}

/**
 * @brief Serialize one optional slice index using a stable missing token.
 * @param prefix Text written before the serialized value.
 * @param index Zero-based slice index, or std::nullopt for a disabled axis.
 * @param output Destination output stream.
 */
void write_optional_slice_index(
    std::string_view prefix,
    const std::optional<std::size_t>& index,
    std::ostream& output
) {
    output << prefix;
    if (index.has_value()) {
        output << index.value();
    } else {
        output << "none";
    }
    output << '\n';
}

/**
 * @brief Serialize one nested primitive-channel pair-count summary.
 * @param label Stable nested section label.
 * @param counts Ordered primitive-channel counts.
 * @param output Destination output stream.
 */
void write_nested_pair_count_summary(
    std::string_view label,
    const hbt::PairCountSummary& counts,
    std::ostream& output
) {
    output << "    " << label << ":\n";

    for (const hbt::PairChannelCount& entry : counts.channels) {
        const hbt::PrimitiveChannel& channel =
            hbt::primitive_channel_definition(entry.channel);
        output << "      " << channel.canonical_name << ": "
               << entry.pair_count << '\n';
    }
}

/**
 * @brief Validate output identities joining startup bins and slice counts.
 * @param slicing Validated production slicing configuration from startup.
 * @param pair_summary Completed pair-processing summary.
 * @throws std::logic_error If summary routing state does not match startup.
 */
void validate_pair_slice_output_identity(
    const hbt::PairSlicingConfig& slicing,
    const hbt::HBTPairProcessingSummary& pair_summary
) {
    const hbt::PairSliceCountSummary& slice_counts =
        pair_summary.total_pair_slice_counts;

    if (!same_pair_slicing_config(slicing, slice_counts.pair_slicing)) {
        throw std::logic_error(
            "analysis output: pair-slicing configuration mismatch"
        );
    }
    if (slice_counts.origin_mode !=
        pair_summary.total_origin_route_counts.origin_mode) {
        throw std::logic_error(
            "analysis output: pair-slice origin-mode mismatch"
        );
    }
}

/**
 * @brief Serialize production slicing bins and run-total slice counts.
 * @param slicing Validated pair-slicing configuration used by this run.
 * @param pair_summary Completed run-total pair-processing summary.
 * @param output Destination output stream.
 */
void write_pair_slice_counts(
    const hbt::PairSlicingConfig& slicing,
    const hbt::HBTPairProcessingSummary& pair_summary,
    std::ostream& output
) {
    validate_pair_slice_output_identity(slicing, pair_summary);

    const std::streamsize previous_precision = output.precision();
    const std::ios::fmtflags previous_flags = output.flags();
    output << std::setprecision(17);

    write_pair_slicing_axis("kt", slicing.kt, output);
    write_pair_slicing_axis("mt", slicing.mt, output);
    output << "pair_slice_counts:\n";

    for (const hbt::PairSliceCountEntry& entry :
         pair_summary.total_pair_slice_counts.entries) {
        write_optional_slice_index(
            "  - kt_slice_index: ",
            entry.kt_slice_index,
            output
        );
        write_optional_slice_index(
            "    mt_slice_index: ",
            entry.mt_slice_index,
            output
        );
        write_nested_pair_count_summary(
            "routed_P_pair_counts_by_primitive_channel",
            entry.origin_counts.routed_P,
            output
        );
        write_nested_pair_count_summary(
            "routed_PR_pair_counts_by_primitive_channel",
            entry.origin_counts.routed_PR,
            output
        );
        write_nested_pair_count_summary(
            "routed_PRD_pair_counts_by_primitive_channel",
            entry.origin_counts.routed_PRD,
            output
        );
    }

    output.flags(previous_flags);
    output.precision(previous_precision);
}

/**
 * @brief Return the stable token for one configured origin-routing mode.
 * @param mode Origin mode to serialize.
 * @return Stable configuration token.
 * @throws std::invalid_argument If mode is not a valid OriginMode value.
 */
std::string_view origin_mode_token(hbt::OriginMode mode) {
    switch (mode) {
    case hbt::OriginMode::Primordial:
        return "primordial";
    case hbt::OriginMode::PrimordialRescattering:
        return "primordial_rescattering";
    case hbt::OriginMode::PrimordialRescatteringDecay:
        return "primordial_rescattering_decay";
    case hbt::OriginMode::All:
        return "all";
    }

    throw std::invalid_argument("analysis output: invalid origin mode");
}

/**
 * @brief Return a stable token for one stored histogram origin.
 * @param origin Physical stored origin identity.
 * @return Stable short origin token.
 * @throws std::invalid_argument If @p origin is invalid.
 */
std::string_view histogram_origin_token(hbt::HistogramOrigin origin) {
    switch (origin) {
        case hbt::HistogramOrigin::Primordial:
            return "P";
        case hbt::HistogramOrigin::PrimordialRescattering:
            return "PR";
        case hbt::HistogramOrigin::PrimordialRescatteringDecay:
            return "PRD";
    }
    throw std::invalid_argument("analysis output: invalid histogram origin");
}

/**
 * @brief Serialize one final-product composition without allocating a name.
 * @param product Final HBT analysis product.
 * @param output Destination output stream.
 */
void write_product_name(
    const hbt::AnalysisProduct& product,
    std::ostream& output
) {
    for (std::size_t index = 0U;
         index < product.primitive_channels.size();
         ++index) {
        if (index != 0U) {
            output << '+';
        }
        output << hbt::primitive_channel_definition(
            product.primitive_channels[index]
        ).canonical_name;
    }
}

/**
 * @brief Serialize one affected histogram-range diagnostic.
 * @param product_index Final-product index in configured order.
 * @param product Final-product composition.
 * @param origin Stored physical origin identity.
 * @param slice_index Flat slice index, or std::nullopt for global state.
 * @param observable Stable logical observable token.
 * @param binning Configured histogram range.
 * @param underflow Underflow count for this logical histogram.
 * @param overflow Overflow count for this logical histogram.
 * @param output Destination output stream.
 */
void write_histogram_range_warning(
    std::size_t product_index,
    const hbt::AnalysisProduct& product,
    hbt::HistogramOrigin origin,
    std::optional<std::size_t> slice_index,
    std::string_view observable,
    const hbt::HistogramBinningConfig& binning,
    std::uint64_t underflow,
    std::uint64_t overflow,
    std::ostream& output
) {
    output << "  - product_index: " << product_index << '\n';
    output << "    product: ";
    write_product_name(product, output);
    output << '\n';
    output << "    origin: " << histogram_origin_token(origin) << '\n';
    if (slice_index.has_value()) {
        output << "    destination: slice\n";
        output << "    flat_slice_index: "
               << slice_index.value() << '\n';
    } else {
        output << "    destination: global\n";
        output << "    flat_slice_index: none\n";
    }
    output << "    observable: " << observable << '\n';
    output << "    range: [" << binning.minimum << ", "
           << binning.maximum << ")\n";
    output << "    underflow: " << underflow << '\n';
    output << "    overflow: " << overflow << '\n';
}

/**
 * @brief Serialize range warnings for one nine-histogram destination.
 * @param product_index Final-product index in configured order.
 * @param product Final-product composition.
 * @param origin Stored physical origin identity.
 * @param slice_index Flat slice index, or std::nullopt for global state.
 * @param set Raw histogram destination to inspect.
 * @param config Validated histogram binning configuration.
 * @param output Destination output stream.
 * @return Number of affected logical histograms serialized.
 */
std::size_t write_histogram_set_range_warnings(
    std::size_t product_index,
    const hbt::AnalysisProduct& product,
    hbt::HistogramOrigin origin,
    std::optional<std::size_t> slice_index,
    const hbt::RawHistogramSet& set,
    const hbt::HBTHistogramConfig& config,
    std::ostream& output
) {
    static constexpr std::array<std::string_view, 4U> osl_tokens{
        "r_out_lcms",
        "r_out_prf",
        "r_side",
        "r_long"
    };
    static constexpr std::array<std::string_view, 2U> radial_tokens{
        "r_radial_lcms",
        "r_radial_prf"
    };
    static constexpr std::array<std::string_view, 3U> delta_t_tokens{
        "delta_t_lab",
        "delta_t_lcms",
        "delta_t_prf"
    };

    std::size_t warning_count = 0U;
    for (std::size_t slot = 0U; slot < osl_tokens.size(); ++slot) {
        const std::uint64_t underflow = set.osl.underflow_counts[slot];
        const std::uint64_t overflow = set.osl.overflow_counts[slot];
        if (underflow != 0U || overflow != 0U) {
            write_histogram_range_warning(
                product_index,
                product,
                origin,
                slice_index,
                osl_tokens[slot],
                config.osl,
                underflow,
                overflow,
                output
            );
            ++warning_count;
        }
    }
    for (std::size_t slot = 0U; slot < radial_tokens.size(); ++slot) {
        const std::uint64_t underflow = set.radial.underflow_counts[slot];
        const std::uint64_t overflow = set.radial.overflow_counts[slot];
        if (underflow != 0U || overflow != 0U) {
            write_histogram_range_warning(
                product_index,
                product,
                origin,
                slice_index,
                radial_tokens[slot],
                config.radial,
                underflow,
                overflow,
                output
            );
            ++warning_count;
        }
    }
    for (std::size_t slot = 0U; slot < delta_t_tokens.size(); ++slot) {
        const std::uint64_t underflow = set.delta_t.underflow_counts[slot];
        const std::uint64_t overflow = set.delta_t.overflow_counts[slot];
        if (underflow != 0U || overflow != 0U) {
            write_histogram_range_warning(
                product_index,
                product,
                origin,
                slice_index,
                delta_t_tokens[slot],
                config.delta_t,
                underflow,
                overflow,
                output
            );
            ++warning_count;
        }
    }
    return warning_count;
}

/**
 * @brief Serialize one aggregate summary of all histogram range warnings.
 * @param config Complete validated HBT configuration.
 * @param state Complete-sample raw histogram state.
 * @param output Destination output stream.
 * @throws std::logic_error If histogram dimensions disagree with config.
 */
void write_histogram_range_warnings(
    const hbt::HBTConfig& config,
    const hbt::RawHistogramState& state,
    std::ostream& output
) {
    hbt::require_raw_histogram_state_layout(config, state);

    const std::streamsize previous_precision = output.precision();
    const std::ios::fmtflags previous_flags = output.flags();
    output << std::setprecision(17);

    output << "histogram_range_warnings:\n";
    std::size_t warning_count = 0U;

    for (std::size_t product_index = 0U;
         product_index < state.products.size();
         ++product_index) {
        const hbt::ProductRawHistogramState& product_state =
            state.products[product_index];
        for (std::size_t origin_index = 0U;
             origin_index < product_state.origins.size();
             ++origin_index) {
            const hbt::HistogramOrigin origin = hbt::raw_histogram_origin_at(
                config.origin_mode,
                origin_index
            );
            const hbt::OriginRawHistogramState& origin_state =
                product_state.origins[origin_index];
            warning_count += write_histogram_set_range_warnings(
                product_index,
                config.selection.products[product_index],
                origin,
                std::nullopt,
                origin_state.global,
                config.histogram_config,
                output
            );
            for (std::size_t slice_index = 0U;
                 slice_index < origin_state.slices.size();
                 ++slice_index) {
                warning_count += write_histogram_set_range_warnings(
                    product_index,
                    config.selection.products[product_index],
                    origin,
                    slice_index,
                    origin_state.slices[slice_index],
                    config.histogram_config,
                    output
                );
            }
        }
    }

    if (warning_count == 0U) {
        output << "  none\n";
    }

    output.flags(previous_flags);
    output.precision(previous_precision);
}

/**
 * @brief Serialize run-total primitive-channel pair counts.
 * @param summary Completed HBT pair-processing summary.
 * @param output Destination output stream.
 */
void write_pair_counts(
    const hbt::HBTPairProcessingSummary& summary,
    std::ostream& output
) {
    write_pair_count_summary(
        "pair_counts_by_primitive_channel",
        summary.total_pair_counts,
        output
    );
    write_pair_count_summary(
        "valid_pair_counts_by_primitive_channel",
        summary.total_valid_pair_counts,
        output
    );
    write_pair_count_summary(
        "numerical_pair_rejections_by_primitive_channel",
        summary.total_numerical_rejection_counts,
        output
    );
    output << "pair_origin_route_mode: "
           << origin_mode_token(
                  summary.total_origin_route_counts.origin_mode
              )
           << '\n';
    write_pair_count_summary(
        "routed_P_pair_counts_by_primitive_channel",
        summary.total_origin_route_counts.routed_P,
        output
    );
    write_pair_count_summary(
        "routed_PR_pair_counts_by_primitive_channel",
        summary.total_origin_route_counts.routed_PR,
        output
    );
    write_pair_count_summary(
        "routed_PRD_pair_counts_by_primitive_channel",
        summary.total_origin_route_counts.routed_PRD,
        output
    );
}

}  // namespace

void write_rejected_particle_report(
    const hbt::RejectedParticleReport& report,
    std::ostream& output
) {
    const std::streamsize previous_precision = output.precision();
    const std::ios::fmtflags previous_flags = output.flags();
    output << std::setprecision(17);

    output << "numerical_particle_rejections: " << report.size() << '\n';
    output << "numerical_particle_rejections_by_reason:\n";
    write_rejection_counts(report, output);
    output << "rejected_particles:\n";

    for (const hbt::RejectedParticleRecord& record : report.records()) {
        output << "  - outer_event_number: "
               << record.outer_event_number << '\n';
        output << "    subevent_id: " << record.subevent_id << '\n';
        output << "    particle_id: " << record.particle_id << '\n';
        output << "    pdg: " << record.raw_pdg << '\n';
        output << "    charge: " << record.raw_charge << '\n';
        output << "    species: "
               << hbt::species_metadata(record.species).ascii_token << '\n';
        output << "    raw_mass_gev: " << record.raw_mass_gev << '\n';
        output << "    ncoll: " << record.ncoll << '\n';
        output << "    time_last_coll_fm: "
               << record.time_last_coll << '\n';
        output << "    momentum_gev: ["
               << record.momentum.x0 << ", "
               << record.momentum.x1 << ", "
               << record.momentum.x2 << ", "
               << record.momentum.x3 << "]\n";
        output << "    raw_position_fm: ["
               << record.raw_position.x0 << ", "
               << record.raw_position.x1 << ", "
               << record.raw_position.x2 << ", "
               << record.raw_position.x3 << "]\n";
        output << "    reason: "
               << rejection_reason_token(record.reason) << '\n';

        if (record.diagnostic_value.has_value()) {
            output << "    diagnostic_value: "
                   << record.diagnostic_value.value() << '\n';
        } else {
            output << "    diagnostic_value: none\n";
        }

        if (record.diagnostic_position.has_value()) {
            const common::FourVector& position =
                record.diagnostic_position.value();
            output << "    diagnostic_position_fm: ["
                   << position.x0 << ", "
                   << position.x1 << ", "
                   << position.x2 << ", "
                   << position.x3 << "]\n";
        } else {
            output << "    diagnostic_position_fm: none\n";
        }
    }

    output.flags(previous_flags);
    output.precision(previous_precision);
}

void write_rejected_pair_report(
    const hbt::RejectedPairReport& report,
    std::ostream& output
) {
    const std::streamsize previous_precision = output.precision();
    const std::ios::fmtflags previous_flags = output.flags();
    output << std::setprecision(17);

    output << "numerical_pair_rejections: " << report.size() << '\n';
    output << "numerical_pair_rejections_by_reason:\n";
    write_pair_rejection_counts(report, output);
    output << "rejected_pairs:\n";

    for (const hbt::RejectedPairRecord& record : report.records()) {
        const hbt::PrimitiveChannel& channel =
            hbt::primitive_channel_definition(record.channel);

        output << "  - outer_event_number: "
               << record.outer_event_number << '\n';
        output << "    subevent_id: " << record.subevent_id << '\n';
        output << "    channel: " << channel.canonical_name << '\n';
        output << "    pair_ordinal_in_channel: "
               << record.pair_ordinal_in_channel << '\n';
        const hbt::SpeciesMetadata& species_a =
            hbt::species_metadata(record.particle_a.species);
        const hbt::SpeciesMetadata& species_b =
            hbt::species_metadata(record.particle_b.species);

        output << "    particle_a_species: "
               << species_a.ascii_token << '\n';
        output << "    particle_a_pdg: "
               << record.particle_a.raw_pdg << '\n';
        output << "    particle_a_charge: "
               << record.particle_a.raw_charge << '\n';
        output << "    particle_a_invariant_mass_gev: "
               << record.particle_a.invariant_mass_gev << '\n';
        output << "    particle_a_momentum_gev: ["
               << record.particle_a.momentum.x0 << ", "
               << record.particle_a.momentum.x1 << ", "
               << record.particle_a.momentum.x2 << ", "
               << record.particle_a.momentum.x3 << "]\n";
        output << "    particle_b_species: "
               << species_b.ascii_token << '\n';
        output << "    particle_b_pdg: "
               << record.particle_b.raw_pdg << '\n';
        output << "    particle_b_charge: "
               << record.particle_b.raw_charge << '\n';
        output << "    particle_b_invariant_mass_gev: "
               << record.particle_b.invariant_mass_gev << '\n';
        output << "    particle_b_momentum_gev: ["
               << record.particle_b.momentum.x0 << ", "
               << record.particle_b.momentum.x1 << ", "
               << record.particle_b.momentum.x2 << ", "
               << record.particle_b.momentum.x3 << "]\n";
        output << "    kt_gev: " << record.kinematics.kt_gev << '\n';
        output << "    mt_gev: " << record.kinematics.mt_gev << '\n';
        output << "    reason: "
               << pair_rejection_reason_token(record.reason) << '\n';
    }

    output.flags(previous_flags);
    output.precision(previous_precision);
}

/**
 * @brief Return the stable output token for one event execution status.
 * @param status Event status to serialize.
 * @return Stable ASCII diagnostic token.
 * @throws std::invalid_argument If @p status is invalid.
 */
std::string_view event_status_token(app::EventStatus status) {
    switch (status) {
        case app::EventStatus::Processed:
            return "processed";
        case app::EventStatus::SkippedDueToEventEmpty:
            return "skipped_due_to_event_empty";
        case app::EventStatus::SkippedDueToEventFailure:
            return "skipped_due_to_event_failure";
    }
    throw std::invalid_argument("analysis output: invalid event status");
}

/**
 * @brief Return the stable output token for one subevent execution status.
 * @param status Subevent status to serialize.
 * @return Stable ASCII diagnostic token.
 * @throws std::invalid_argument If @p status is invalid.
 */
std::string_view subevent_status_token(app::SubeventStatus status) {
    switch (status) {
        case app::SubeventStatus::Processed:
            return "processed";
        case app::SubeventStatus::SkippedDueToSubeventEmpty:
            return "skipped_due_to_subevent_empty";
        case app::SubeventStatus::SkippedDueToSubeventFailure:
            return "skipped_due_to_subevent_failure";
    }
    throw std::invalid_argument("analysis output: invalid subevent status");
}

/**
 * @brief Serialize aggregate execution-status counts and non-processed details.
 * @param summary Completed preparation summary in canonical event-major order.
 * @param output Destination output stream.
 */
void write_execution_statuses(
    const app::HBTEventPreparationSummary& summary,
    std::ostream& output
) {
    output << "event_status_counts:\n";
    output << "  processed: "
           << summary.event_status_counts.processed << '\n';
    output << "  skipped_due_to_event_empty: "
           << summary.event_status_counts.skipped_due_to_event_empty
           << '\n';
    output << "  skipped_due_to_event_failure: "
           << summary.event_status_counts.skipped_due_to_event_failure
           << '\n';
    output << "subevent_status_counts:\n";
    output << "  processed: "
           << summary.subevent_status_counts.processed << '\n';
    output << "  skipped_due_to_subevent_empty: "
           << summary.subevent_status_counts
                  .skipped_due_to_subevent_empty
           << '\n';
    output << "  skipped_due_to_subevent_failure: "
           << summary.subevent_status_counts
                  .skipped_due_to_subevent_failure
           << '\n';

    output << "event_status_diagnostics:\n";
    bool wrote_event = false;
    for (const app::HBTEventExecutionSummary& event : summary.events) {
        if (event.status == app::EventStatus::Processed) {
            continue;
        }
        wrote_event = true;
        output << "  - outer_event_number: "
               << event.outer_event_number << '\n';
        output << "    status: " << event_status_token(event.status) << '\n';
        if (event.status == app::EventStatus::SkippedDueToEventFailure &&
            !event.diagnostic.empty()) {
            output << "    reason: " << std::quoted(event.diagnostic) << '\n';
        }
    }
    if (!wrote_event) {
        output << "  []\n";
    }

    output << "subevent_status_diagnostics:\n";
    bool wrote_subevent = false;
    for (const app::HBTSubeventPreparationSummary& subevent :
         summary.subevents) {
        if (subevent.status == app::SubeventStatus::Processed) {
            continue;
        }
        wrote_subevent = true;
        output << "  - outer_event_number: "
               << subevent.outer_event_number << '\n';
        output << "    subevent_id: " << subevent.subevent_id << '\n';
        output << "    status: "
               << subevent_status_token(subevent.status) << '\n';
        if (subevent.status ==
                app::SubeventStatus::SkippedDueToSubeventFailure &&
            !subevent.diagnostic.empty()) {
            output << "    reason: "
                   << std::quoted(subevent.diagnostic) << '\n';
        }
    }
    if (!wrote_subevent) {
        output << "  []\n";
    }

}

void write_analysis_output(
    const app::AnalysisRunSummary& result,
    std::ostream& output
) {
    output << "analysis_run_summary:\n";
    output << "hbt_enabled: "
           << (result.startup.run_config.hbt_enabled ? "true" : "false")
           << '\n';

    if (result.hbt_event_preparation.has_value()) {
        const app::HBTEventPreparationSummary& summary =
            result.hbt_event_preparation.value();

        output << "outer_events_processed: "
               << summary.outer_events_processed << '\n';
        output << "subevents_processed: "
               << summary.subevents_processed << '\n';
        output << "raw_particles: " << summary.raw_particles << '\n';
        output << "unsupported_species: "
               << summary.unsupported_species << '\n';
        output << "unrequired_species: "
               << summary.unrequired_species << '\n';
        output << "particle_acceptance_rejections: "
               << summary.particle_acceptance_rejections << '\n';
        output << "origin_rejections: "
               << summary.origin_rejections << '\n';
        output << "accepted_particles: "
               << summary.accepted_particles << '\n';
        output << "emission_points_sampler: "
               << summary.emission_points.sampler << '\n';
        output << "emission_points_propagation: "
               << summary.emission_points.propagation << '\n';
        output << "emission_points_afterburner: "
               << summary.emission_points.afterburner << '\n';
        if (!summary.events.empty()) {
            write_execution_statuses(summary, output);
        }
        write_rejected_particle_report(summary.numerical_rejections, output);
    }

    if (result.hbt_pair_processing.has_value()) {
        if (!result.startup.hbt_config.has_value()) {
            throw std::logic_error(
                "analysis output: pair summary requires HBT configuration"
            );
        }

        const hbt::HBTPairProcessingSummary& pair_summary =
            result.hbt_pair_processing.value();
        write_pair_counts(pair_summary, output);
        write_pair_slice_counts(
            result.startup.hbt_config->pair_slicing,
            pair_summary,
            output
        );
        write_rejected_pair_report(
            pair_summary.numerical_rejections,
            output
        );
    }

    if (result.hbt_raw_histograms.has_value()) {
        if (!result.startup.hbt_config.has_value()) {
            throw std::logic_error(
                "analysis output: histogram state requires HBT configuration"
            );
        }
        write_histogram_range_warnings(
            result.startup.hbt_config.value(),
            result.hbt_raw_histograms.value(),
            output
        );
    }
}

void write_production_output(
    const app::AnalysisRunSummary& result
) {
    write_hbt_production_output(
        result,
        result.startup.run_config.output_path
    );
}

}  // namespace output
