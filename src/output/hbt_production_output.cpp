/**
 * @file hbt_production_output.cpp
 * @brief Filesystem boundary for post-sample HBT production output.
 */

#include "output/hbt_production_output.h"

#include "hbt/channels/channel_catalog.h"
#include "hbt/fits/histogram_analysis.h"
#include "hbt/fits/fit_results.h"
#include "hbt/histograms/raw_histograms.h"

#include <array>
#include <cstddef>
#include <filesystem>
#include <fstream>
#include <iomanip>
#include <optional>
#include <stdexcept>
#include <string>

namespace output {
namespace {

/**
 * @brief Canonical presentation location for one logical histogram slot.
 */
struct ObservableLocation {
    const char* frame;      ///< LAB, LCMS, or PRF directory token.
    const char* family;     ///< osl, radial, or dt directory token.
    const char* observable; ///< Physical observable directory token.
};

/**
 * @brief Return the canonical presentation location of one OSL slot.
 * @param slot Stable raw-histogram OSL slot index.
 * @return Frame/family/observable directory identity.
 * @throws std::invalid_argument If @p slot is outside the four OSL slots.
 */
ObservableLocation osl_location(std::size_t slot) {
    switch (static_cast<hbt::OSLHistogramSlot>(slot)) {
        case hbt::OSLHistogramSlot::ROutLcms:
            return {"LCMS", "osl", "r_out"};
        case hbt::OSLHistogramSlot::ROutPrf:
            return {"PRF", "osl", "r_out"};
        case hbt::OSLHistogramSlot::RSide:
            return {"LCMS", "osl", "r_side"};
        case hbt::OSLHistogramSlot::RLong:
            return {"LCMS", "osl", "r_long"};
        case hbt::OSLHistogramSlot::Count:
            break;
    }
    throw std::invalid_argument("output: invalid OSL histogram slot");
}

/**
 * @brief Return the canonical presentation location of one radial slot.
 * @param slot Stable raw-histogram radial slot index.
 * @return Frame/family/observable directory identity.
 * @throws std::invalid_argument If @p slot is outside the two radial slots.
 */
ObservableLocation radial_location(std::size_t slot) {
    switch (static_cast<hbt::RadialHistogramSlot>(slot)) {
        case hbt::RadialHistogramSlot::RadialLcms:
            return {"LCMS", "radial", "r_radial"};
        case hbt::RadialHistogramSlot::RadialPrf:
            return {"PRF", "radial", "r_radial"};
        case hbt::RadialHistogramSlot::Count:
            break;
    }
    throw std::invalid_argument("output: invalid radial histogram slot");
}

/**
 * @brief Return the canonical presentation location of one delta-t slot.
 * @param slot Stable raw-histogram delta-t slot index.
 * @return Frame/family/observable directory identity.
 * @throws std::invalid_argument If @p slot is outside the three time slots.
 */
ObservableLocation delta_t_location(std::size_t slot) {
    switch (static_cast<hbt::DeltaTHistogramSlot>(slot)) {
        case hbt::DeltaTHistogramSlot::Lab:
            return {"LAB", "dt", "delta_t"};
        case hbt::DeltaTHistogramSlot::Lcms:
            return {"LCMS", "dt", "delta_t"};
        case hbt::DeltaTHistogramSlot::Prf:
            return {"PRF", "dt", "delta_t"};
        case hbt::DeltaTHistogramSlot::Count:
            break;
    }
    throw std::invalid_argument("output: invalid delta-t histogram slot");
}

/**
 * @brief Return stable output token for one physical origin.
 * @param origin Physical nested origin identity.
 * @return P, PR, or PRD static token.
 * @throws std::invalid_argument If @p origin is invalid.
 */
const char* origin_token(hbt::HistogramOrigin origin) {
    switch (origin) {
        case hbt::HistogramOrigin::Primordial:
            return "P";
        case hbt::HistogramOrigin::PrimordialRescattering:
            return "PR";
        case hbt::HistogramOrigin::PrimordialRescatteringDecay:
            return "PRD";
    }
    throw std::invalid_argument("output: invalid HistogramOrigin");
}

/**
 * @brief Construct one product directory token from its stable index.
 * @param product_index Zero-based configured final-product index.
 * @return Canonical serialization-only product directory token.
 */
std::string product_token(std::size_t product_index) {
    return "product_" + std::to_string(product_index);
}

/**
 * @brief Construct one global/slice directory token.
 * @param slice_index std::nullopt for global, otherwise flat slice index.
 * @return Canonical serialization-only scope directory token.
 */
std::string scope_token(std::optional<std::size_t> slice_index) {
    if (!slice_index.has_value()) {
        return "global";
    }
    return "slice_" + std::to_string(slice_index.value());
}

/**
 * @brief Return one complete canonical observable directory.
 * @param root Explicit output root.
 * @param product_index Stable configured product index.
 * @param slice_index Global or stable flat-slice index.
 * @param origin Physical origin token.
 * @param location Frame/family/observable serialization identity.
 * @return Path ending at the physical observable directory.
 */
std::filesystem::path observable_directory(
    const std::filesystem::path& root,
    std::size_t product_index,
    std::optional<std::size_t> slice_index,
    const char* origin,
    const ObservableLocation& location
) {
    return root /
        product_token(product_index) /
        scope_token(slice_index) /
        origin /
        location.frame /
        location.family /
        location.observable;
}

/**
 * @brief Open one canonical output file after creating its parent directory.
 * @param directory Observable directory to materialize.
 * @param filename Canonical filename inside the observable directory.
 * @return Writable truncating file stream.
 * @throws std::runtime_error If the file cannot be opened.
 * @throws std::filesystem::filesystem_error If directory creation fails.
 */
std::ofstream open_output_file(
    const std::filesystem::path& directory,
    const char* filename
) {
    std::filesystem::create_directories(directory);
    const std::filesystem::path path = directory / filename;
    std::ofstream output(path, std::ios::out | std::ios::trunc);
    if (!output) {
        throw std::runtime_error(
            "output: failed to open production file '" +
            path.string() + "'"
        );
    }
    output << std::setprecision(17);
    return output;
}

/**
 * @brief Require one completed file stream remained writable.
 * @param output Stream to check after serialization.
 * @param directory Owning observable directory for diagnostics.
 * @param filename Filename written in @p directory.
 * @throws std::runtime_error If a write failed.
 */
void require_written(
    const std::ofstream& output,
    const std::filesystem::path& directory,
    const char* filename
) {
    if (!output) {
        throw std::runtime_error(
            "output: failed while writing production file '" +
            (directory / filename).string() + "'"
        );
    }
}

/**
 * @brief Serialize resolved product/channel/species run metadata.
 * @param root Explicit production output root.
 * @param config Validated HBT configuration with configured product identity.
 * @throws std::logic_error If configured product metadata is unavailable or
 *         unsafe for the unquoted CSV representation.
 */
void write_product_catalog(
    const std::filesystem::path& root,
    const hbt::HBTConfig& config
) {
    std::ofstream output = open_output_file(root, "product_catalog.csv");
    output << "product_index,configured_expression,canonical_expression,"
           << "channel_index,primitive_channel,species_a,species_b\n";

    for (std::size_t product_index = 0U;
         product_index < config.selection.products.size();
         ++product_index) {
        const hbt::AnalysisProduct& product =
            config.selection.products[product_index];
        if (product.configured_expression.empty() ||
            product.configured_expression.find(',') != std::string::npos ||
            product.configured_expression.find('\n') != std::string::npos ||
            product.configured_expression.find('\r') != std::string::npos) {
            throw std::logic_error(
                "output: configured product expression is unavailable or "
                "not CSV-safe"
            );
        }

        std::string canonical_expression;
        for (std::size_t channel_index = 0U;
             channel_index < product.primitive_channels.size();
             ++channel_index) {
            if (channel_index != 0U) {
                canonical_expression += '+';
            }
            canonical_expression += hbt::primitive_channel_definition(
                product.primitive_channels[channel_index]
            ).canonical_name;
        }

        for (std::size_t channel_index = 0U;
             channel_index < product.primitive_channels.size();
             ++channel_index) {
            const hbt::PrimitiveChannel& channel =
                hbt::primitive_channel_definition(
                    product.primitive_channels[channel_index]
                );
            output << product_index << ','
                   << product.configured_expression << ','
                   << canonical_expression << ','
                   << channel_index << ','
                   << channel.canonical_name << ','
                   << hbt::species_metadata(channel.species_a).ascii_token
                   << ','
                   << hbt::species_metadata(channel.species_b).ascii_token
                   << '\n';
        }
    }
    require_written(output, root, "product_catalog.csv");
}

/**
 * @brief Serialize an optional double as one CSV field.
 * @param output Destination stream.
 * @param value Optional finite scientific or diagnostic value.
 */
void write_optional_double(
    std::ostream& output,
    const std::optional<double>& value
) {
    if (value.has_value()) {
        output << value.value();
    }
}

/**
 * @brief Serialize one boolean as a stable lowercase token.
 * @param output Destination stream.
 * @param value Boolean value to serialize.
 */
void write_bool(std::ostream& output, bool value) {
    output << (value ? "true" : "false");
}

/**
 * @brief Serialize all stable identity fields shared by local result rows.
 * @param output Destination CSV stream.
 * @param product_index Stable configured final-product index.
 * @param origin Physical P/PR/PRD token.
 * @param slice_index Global or flat-slice identity.
 * @param location Frame/family/observable serialization identity.
 */
void write_identity_fields(
    std::ostream& output,
    std::size_t product_index,
    const char* origin,
    std::optional<std::size_t> slice_index,
    const ObservableLocation& location
) {
    output << product_index << ',' << origin << ',';
    if (slice_index.has_value()) {
        output << "slice," << slice_index.value();
    } else {
        output << "global,";
    }
    output << ',' << location.frame << ',' << location.family << ','
           << location.observable;
}

/**
 * @brief Write MINOS side status fields for one fitted parameter.
 * @param output Destination CSV stream.
 * @param diagnostic Stable parameter MINOS diagnostic.
 */
void write_minos_fields(
    std::ostream& output,
    const hbt::MinosDiagnostic& diagnostic
) {
    output << ',';
    write_bool(output, diagnostic.attempted);
    output << ',';
    write_bool(output, diagnostic.lower_valid);
    output << ',';
    write_bool(output, diagnostic.upper_valid);
    output << ',';
    write_bool(output, diagnostic.at_lower_limit);
    output << ',';
    write_bool(output, diagnostic.at_upper_limit);
    output << ',';
    write_bool(output, diagnostic.lower_call_limit);
    output << ',';
    write_bool(output, diagnostic.upper_call_limit);
    output << ',';
    write_bool(output, diagnostic.lower_new_minimum);
    output << ',';
    write_bool(output, diagnostic.upper_new_minimum);
}

/**
 * @brief Write MIGRAD fields for one independent start.
 * @param output Destination CSV stream.
 * @param diagnostic Stable MIGRAD diagnostic.
 */
void write_migrad_fields(
    std::ostream& output,
    const hbt::MigradDiagnostic& diagnostic
) {
    output << ',';
    write_bool(output, diagnostic.attempted);
    output << ',';
    write_bool(output, diagnostic.valid);
    output << ',';
    write_bool(output, diagnostic.reached_call_limit);
    output << ',';
    write_bool(output, diagnostic.above_max_edm);
    output << ',';
    write_bool(output, diagnostic.objective_failure);
    output << ',' << diagnostic.function_calls << ',';
    write_optional_double(output, diagnostic.q_min);
}

/**
 * @brief Write the common parameter-table header.
 * @param output Destination CSV stream.
 */
void write_parameter_header(std::ostream& output) {
    output
        << "product_index,origin,scope,slice_index,frame,family,observable,"
        << "model,parameter,fully_valid,failure_reason,q_min,value,"
        << "error_minus,error_plus,migrad_attempted,migrad_valid,"
        << "migrad_call_limit,migrad_above_max_edm,objective_failure,"
        << "migrad_function_calls,migrad_q_min,"
        << "minos_attempted,minos_lower_valid,minos_upper_valid,"
        << "minos_at_lower_limit,minos_at_upper_limit,"
        << "minos_lower_call_limit,minos_upper_call_limit,"
        << "minos_lower_new_minimum,minos_upper_new_minimum,"
        << "selected_start,exact_q_tie,tail_below_core,"
        << "start_a_attempted,start_a_valid,start_a_call_limit,"
        << "start_a_above_max_edm,start_a_objective_failure,"
        << "start_a_function_calls,start_a_q_min,"
        << "start_b_attempted,start_b_valid,start_b_call_limit,"
        << "start_b_above_max_edm,start_b_objective_failure,"
        << "start_b_function_calls,start_b_q_min\n";
}

/**
 * @brief Return one unattempted start diagnostic for Gaussian CSV padding.
 * @return Stable unattempted MIGRAD state.
 */
hbt::MigradDiagnostic unattempted_migrad() {
    return {false, false, false, false, false, 0, std::nullopt};
}

/**
 * @brief Write one complete Gaussian parameter/status row.
 * @param output Destination CSV stream.
 * @param product_index Stable product index.
 * @param origin Physical origin token.
 * @param slice_index Global or flat-slice identity.
 * @param location Observable presentation identity.
 * @param result Completed independent Gaussian fit result.
 */
void write_gaussian_row(
    std::ostream& output,
    std::size_t product_index,
    const char* origin,
    std::optional<std::size_t> slice_index,
    const ObservableLocation& location,
    const hbt::GaussianFitResult& result
) {
    write_identity_fields(
        output,
        product_index,
        origin,
        slice_index,
        location
    );
    output << ",gaussian,R,";
    write_bool(output, result.fully_valid);
    output << ',' << hbt::fit_failure_reason_token(result.failure_reason)
           << ',';
    write_optional_double(output, result.q_min);
    output << ',';
    if (result.fully_valid && result.radius.has_value()) {
        output << result.radius->value << ','
               << result.radius->lower_error << ','
               << result.radius->upper_error;
    } else {
        output << ",,";
    }
    write_migrad_fields(output, result.migrad);
    write_minos_fields(output, result.minos_radius);
    output << ",,,";
    write_migrad_fields(output, result.migrad);
    write_migrad_fields(output, unattempted_migrad());
    output << '\n';
}

/**
 * @brief Write mixed-model fields shared by each parameter row.
 * @param output Destination CSV stream.
 * @param result Completed independent mixed-model fit result.
 */
void write_mixed_shared_fields(
    std::ostream& output,
    const hbt::MixedFitResult& result
) {
    output << ',';
    if (result.selected_start.has_value()) {
        output << (result.selected_start.value() == 0U ? 'A' : 'B');
    }
    output << ',';
    if (result.selected_start.has_value()) {
        write_bool(output, result.exact_q_tie);
    }
    output << ',';
    if (result.selected_start.has_value()) {
        write_bool(output, result.tail_below_core);
    }
    write_migrad_fields(output, result.starts[0U]);
    write_migrad_fields(output, result.starts[1U]);
}

/**
 * @brief Write one mixed-model parameter/status row.
 * @param output Destination CSV stream.
 * @param product_index Stable product index.
 * @param origin Physical origin token.
 * @param slice_index Global or flat-slice identity.
 * @param location Observable presentation identity.
 * @param result Completed independent mixed-model fit result.
 * @param parameter Stable physical parameter name.
 * @param estimate Published physical parameter estimate when fully valid.
 * @param minos Matching parameter MINOS diagnostics.
 */
void write_mixed_row(
    std::ostream& output,
    std::size_t product_index,
    const char* origin,
    std::optional<std::size_t> slice_index,
    const ObservableLocation& location,
    const hbt::MixedFitResult& result,
    const char* parameter,
    const std::optional<hbt::FitParameterEstimate>& estimate,
    const hbt::MinosDiagnostic& minos
) {
    write_identity_fields(
        output,
        product_index,
        origin,
        slice_index,
        location
    );
    output << ",gaussian_plus_exponential," << parameter << ',';
    write_bool(output, result.fully_valid);
    output << ',' << hbt::fit_failure_reason_token(result.failure_reason)
           << ',';
    write_optional_double(output, result.q_min);
    output << ',';
    if (estimate.has_value()) {
        output << estimate->value << ','
               << estimate->lower_error << ','
               << estimate->upper_error;
    } else {
        output << ",,";
    }
    write_migrad_fields(output, result.selected_migrad);
    write_minos_fields(output, minos);
    write_mixed_shared_fields(output, result);
    output << '\n';
}

/**
 * @brief Serialize one shape histogram distribution and fit-state table.
 * @param root Explicit production root.
 * @param product_index Stable configured final-product index.
 * @param origin Physical origin token.
 * @param slice_index Global or flat-slice identity.
 * @param location Frame/family/observable serialization identity.
 * @param result Completed derived shape result.
 */
void write_shape_result(
    const std::filesystem::path& root,
    std::size_t product_index,
    const char* origin,
    std::optional<std::size_t> slice_index,
    const ObservableLocation& location,
    const hbt::ShapeHistogramResult& result
) {
    const std::filesystem::path directory = observable_directory(
        root,
        product_index,
        slice_index,
        origin,
        location
    );

    if (!result.normalized_bins.empty()) {
        std::ofstream distribution = open_output_file(
            directory,
            "distribution.csv"
        );
        distribution
            << "bin_index,lower_edge,upper_edge,center,pdf,d_pdf";
        if (result.gaussian.fully_valid) {
            distribution << ",gaussian_fit_pdf";
        }
        if (result.mixed.fully_valid) {
            distribution << ",mixed_fit_pdf";
        }
        distribution << '\n';

        for (std::size_t index = 0U;
             index < result.normalized_bins.size();
             ++index) {
            const hbt::NormalizedHistogramBin& bin =
                result.normalized_bins[index];
            distribution << bin.bin_index << ','
                         << bin.lower_edge << ','
                         << bin.upper_edge << ','
                         << bin.center << ','
                         << bin.pdf << ','
                         << bin.counting_error_pdf;
            if (result.gaussian.fully_valid) {
                distribution << ',' << result.gaussian.fitted_pdf[index];
            }
            if (result.mixed.fully_valid) {
                distribution << ',' << result.mixed.fitted_pdf[index];
            }
            distribution << '\n';
        }
        require_written(distribution, directory, "distribution.csv");
    }

    std::ofstream parameters = open_output_file(
        directory,
        "fit_parameters.csv"
    );
    write_parameter_header(parameters);
    write_gaussian_row(
        parameters,
        product_index,
        origin,
        slice_index,
        location,
        result.gaussian
    );
    const std::optional<hbt::FitParameterEstimate> core_radius =
        result.mixed.fully_valid ? result.mixed.core_radius : std::nullopt;
    const std::optional<hbt::FitParameterEstimate> tail_radius =
        result.mixed.fully_valid ? result.mixed.tail_radius : std::nullopt;
    const std::optional<hbt::FitParameterEstimate> core_fraction =
        result.mixed.fully_valid ? result.mixed.core_fraction : std::nullopt;
    write_mixed_row(
        parameters,
        product_index,
        origin,
        slice_index,
        location,
        result.mixed,
        "R_core",
        core_radius,
        result.mixed.minos_core_radius
    );
    write_mixed_row(
        parameters,
        product_index,
        origin,
        slice_index,
        location,
        result.mixed,
        "R_tail",
        tail_radius,
        result.mixed.minos_tail_radius
    );
    write_mixed_row(
        parameters,
        product_index,
        origin,
        slice_index,
        location,
        result.mixed,
        "f_core",
        core_fraction,
        result.mixed.minos_core_fraction
    );
    require_written(parameters, directory, "fit_parameters.csv");
}

/**
 * @brief Serialize one delta-t distribution and moment/status table.
 * @param root Explicit production root.
 * @param product_index Stable configured final-product index.
 * @param origin Physical origin token.
 * @param slice_index Global or flat-slice identity.
 * @param location Frame/family/observable serialization identity.
 * @param result Completed derived signed delta-t result.
 */
void write_delta_t_result(
    const std::filesystem::path& root,
    std::size_t product_index,
    const char* origin,
    std::optional<std::size_t> slice_index,
    const ObservableLocation& location,
    const hbt::DeltaTHistogramResult& result
) {
    const std::filesystem::path directory = observable_directory(
        root,
        product_index,
        slice_index,
        origin,
        location
    );

    if (!result.normalized_bins.empty()) {
        std::ofstream distribution = open_output_file(
            directory,
            "distribution.csv"
        );
        distribution
            << "bin_index,lower_edge,upper_edge,center,pdf,d_pdf\n";
        for (const hbt::NormalizedHistogramBin& bin :
             result.normalized_bins) {
            distribution << bin.bin_index << ','
                         << bin.lower_edge << ','
                         << bin.upper_edge << ','
                         << bin.center << ','
                         << bin.pdf << ','
                         << bin.counting_error_pdf << '\n';
        }
        require_written(distribution, directory, "distribution.csv");
    }

    std::ofstream statistics = open_output_file(
        directory,
        "statistics.csv"
    );
    statistics
        << "product_index,origin,scope,slice_index,frame,family,observable,"
        << "status,N_selected,mean_dt,sigma_dt,error_sigma_dt\n";
    write_identity_fields(
        statistics,
        product_index,
        origin,
        slice_index,
        location
    );
    statistics << ',' << hbt::delta_t_status_token(result.status)
               << ',' << result.selected_count << ',';
    write_optional_double(statistics, result.mean);
    statistics << ',';
    write_optional_double(statistics, result.sigma);
    statistics << ',';
    write_optional_double(statistics, result.sigma_error);
    statistics << '\n';
    require_written(statistics, directory, "statistics.csv");
}

/**
 * @brief Serialize all nine canonical observables of one derived destination.
 * @param root Explicit production root.
 * @param product_index Stable configured product index.
 * @param origin Physical P/PR/PRD token.
 * @param slice_index Global or flat-slice identity.
 * @param derived Completed nine-histogram derived state.
 */
void write_derived_set(
    const std::filesystem::path& root,
    std::size_t product_index,
    const char* origin,
    std::optional<std::size_t> slice_index,
    const hbt::DerivedHistogramSet& derived
) {
    for (std::size_t slot = 0U; slot < derived.osl.size(); ++slot) {
        write_shape_result(
            root,
            product_index,
            origin,
            slice_index,
            osl_location(slot),
            derived.osl[slot]
        );
    }
    for (std::size_t slot = 0U; slot < derived.radial.size(); ++slot) {
        write_shape_result(
            root,
            product_index,
            origin,
            slice_index,
            radial_location(slot),
            derived.radial[slot]
        );
    }
    for (std::size_t slot = 0U; slot < derived.delta_t.size(); ++slot) {
        write_delta_t_result(
            root,
            product_index,
            origin,
            slice_index,
            delta_t_location(slot),
            derived.delta_t[slot]
        );
    }
}

/**
 * @brief Require explicit output root to be absent or empty before writing.
 * @param root Caller-supplied output root.
 * @throws std::invalid_argument If root is empty.
 * @throws std::runtime_error If an existing root contains entries.
 * @throws std::filesystem::filesystem_error If inspection fails.
 */
void require_clean_output_root(const std::filesystem::path& root) {
    if (root.empty()) {
        throw std::invalid_argument(
            "output: production root must not be empty"
        );
    }
    if (!std::filesystem::exists(root)) {
        return;
    }
    if (!std::filesystem::is_directory(root)) {
        throw std::runtime_error(
            "output: production root exists and is not a directory"
        );
    }
    if (std::filesystem::directory_iterator(root) !=
        std::filesystem::directory_iterator{}) {
        throw std::runtime_error(
            "output: production root must be absent or empty"
        );
    }
}

}  // namespace

void write_hbt_production_output(
    const app::AnalysisRunSummary& result,
    const std::filesystem::path& output_root
) {
    if (!result.startup.run_config.hbt_enabled) {
        return;
    }
    if (!result.startup.hbt_config.has_value() ||
        !result.hbt_raw_histograms.has_value() ||
        !result.hbt_histogram_analysis.has_value()) {
        throw std::logic_error(
            "output: enabled HBT production state is incomplete"
        );
    }

    hbt::require_histogram_analysis_layout(
        result.startup.hbt_config.value(),
        result.hbt_raw_histograms.value(),
        result.hbt_histogram_analysis.value()
    );
    require_clean_output_root(output_root);

    const hbt::HBTConfig& config = result.startup.hbt_config.value();
    const hbt::HistogramAnalysisState& derived =
        result.hbt_histogram_analysis.value();

    write_product_catalog(output_root, config);

    for (std::size_t product = 0U;
         product < derived.products.size();
         ++product) {
        const hbt::ProductDerivedHistogramState& product_state =
            derived.products[product];
        for (std::size_t origin = 0U;
             origin < product_state.origins.size();
             ++origin) {
            const char* token = origin_token(
                hbt::raw_histogram_origin_at(config.origin_mode, origin)
            );
            const hbt::OriginDerivedHistogramState& origin_state =
                product_state.origins[origin];
            write_derived_set(
                output_root,
                product,
                token,
                std::nullopt,
                origin_state.global
            );
            for (std::size_t slice = 0U;
                 slice < origin_state.slices.size();
                 ++slice) {
                write_derived_set(
                    output_root,
                    product,
                    token,
                    slice,
                    origin_state.slices[slice]
                );
            }
        }
    }
}

}  // namespace output
