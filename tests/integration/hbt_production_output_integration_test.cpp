/**
 * @file hbt_production_output_integration_test.cpp
 * @brief Cross-module raw-histogram to production-output integration.
 */

#include "hbt/fits/binned_models.h"
#include "hbt/fits/histogram_analysis.h"
#include "hbt/histograms/raw_histograms.h"
#include "hbt/startup/hbt_startup_builder.h"
#include "output/analysis_output_writer.h"

#include <cmath>
#include <cstddef>
#include <cstdint>
#include <cstdlib>
#include <filesystem>
#include <fstream>
#include <iostream>
#include <optional>
#include <iterator>
#include <sstream>
#include <stdexcept>
#include <string>
#include <utility>
#include <vector>

namespace {

/**
 * @brief Report one failed post-sample integration condition.
 * @param message Description of the violated contract.
 * @return Always false.
 */
bool fail(const char* message) {
    std::cerr << "hbt_production_output_integration_test: " << message << ".\n";
    return false;
}

/**
 * @brief Construct one compact validated-style HBT configuration.
 * @return One product, one origin, configured mT slices, and explicit binning.
 */
hbt::HBTConfig make_config() {
    return {
        {{hbt::AnalysisProduct{
            {
                hbt::PrimitiveChannelId::PiPlusPiPlus,
                hbt::PrimitiveChannelId::PiMinusPiMinus
            },
            "pi_plus_pi_plus + pi_minus_pi_minus"
        }}},
        {
            hbt::LongitudinalVariable::Pseudorapidity,
            {0.8, 0.14, 4.0},
            {0.8, 0.4, 1.4},
            {0.8, 0.5, 4.05},
            {0.8, 1.0, 10000.0},
            {0.8, 0.3, 10000.0}
        },
        {{false, {}}, {true, {0.5, 0.7, 0.9}}},
        {
            {20U, 0.0, 10.0, 2.0},
            {20U, 0.0, 10.0, 2.0},
            {8U, -4.0, 4.0, 1.0}
        },
        hbt::OriginMode::Primordial
    };
}

/**
 * @brief Convert normalized probabilities to positive deterministic counts.
 * @param probabilities Exact-bin probabilities.
 * @param total Nominal count scale.
 * @return Integer raw counts with no artificial empty tail.
 */
std::vector<std::uint64_t> make_counts(
    const std::vector<double>& probabilities,
    std::uint64_t total
) {
    std::vector<std::uint64_t> counts;
    counts.reserve(probabilities.size());
    for (const double probability : probabilities) {
        const auto rounded = static_cast<std::uint64_t>(
            std::llround(probability * static_cast<double>(total))
        );
        counts.push_back(rounded == 0U ? 1U : rounded);
    }
    return counts;
}

/**
 * @brief Read one small generated text file completely.
 * @param path Existing file to read.
 * @return Complete file text.
 * @throws std::runtime_error If the file cannot be opened.
 */
std::string read_text(const std::filesystem::path& path) {
    std::ifstream input(path);
    if (!input) {
        throw std::runtime_error("failed to open generated integration file");
    }
    return std::string(
        std::istreambuf_iterator<char>(input),
        std::istreambuf_iterator<char>()
    );
}

/**
 * @brief Count simple unquoted CSV fields in one line.
 * @param line CSV line written by the production writer.
 * @return Number of comma-separated fields.
 */
std::size_t csv_field_count(const std::string& line) {
    std::size_t fields = 1U;
    for (const char character : line) {
        if (character == ',') {
            ++fields;
        }
    }
    return fields;
}

/**
 * @brief Verify every CSV row has the same cardinality as its header.
 * @param path Generated CSV file.
 * @return true when all non-empty lines are column-aligned.
 */
bool csv_is_rectangular(const std::filesystem::path& path) {
    std::ifstream input(path);
    std::string line;
    if (!std::getline(input, line)) {
        return false;
    }
    const std::size_t fields = csv_field_count(line);
    while (std::getline(input, line)) {
        if (!line.empty() && csv_field_count(line) != fields) {
            return false;
        }
    }
    return true;
}

/**
 * @brief Verify the complete post-sample analysis and output ownership tree.
 * @return true when raw counts are reused unchanged and output is canonical.
 */
bool verify_analysis_and_output() {
    const hbt::HBTConfig config = make_config();
    hbt::RawHistogramState raw = hbt::make_zero_raw_histogram_state(config);
    hbt::RawHistogramSet& global = raw.products[0U].origins[0U].global;

    const hbt::StatisticalRegion model_region{0U, 19U, 1U};
    const std::vector<double> probabilities = hbt::mixed_bin_probabilities(
        hbt::FitObservableFamily::OSL,
        config.histogram_config.osl,
        model_region,
        0.85,
        2.4,
        0.68
    );
    const std::vector<std::uint64_t> shape_counts =
        make_counts(probabilities, 800000U);
    for (std::size_t bin = 0U; bin < shape_counts.size(); ++bin) {
        global.osl.bins[bin] = shape_counts[bin];
    }
    const std::vector<std::uint64_t> delta_counts{
        0U, 0U, 3U, 7U, 7U, 3U, 0U, 0U
    };
    for (std::size_t bin = 0U; bin < delta_counts.size(); ++bin) {
        global.delta_t.bins[bin] = delta_counts[bin];
    }

    // Populate one radial mT slice below the historical 10000-count threshold.
    // With the retained threshold machinery set to zero, this slice must now
    // attempt all Gaussian and mixed fits rather than being vetoed pre-fit.
    hbt::RawHistogramSet& low_stat_slice =
        raw.products[0U].origins[0U].slices[0U];
    const std::vector<double> radial_probabilities =
        hbt::gaussian_bin_probabilities(
            hbt::FitObservableFamily::Radial,
            config.histogram_config.radial,
            model_region,
            1.2
        );
    const std::vector<std::uint64_t> radial_counts =
        make_counts(radial_probabilities, 5000U);
    for (std::size_t bin = 0U; bin < radial_counts.size(); ++bin) {
        low_stat_slice.radial.bins[bin] = radial_counts[bin];
    }

    const std::vector<std::uint64_t> original_osl = global.osl.bins;
    const std::vector<std::uint64_t> original_delta_t = global.delta_t.bins;
    const hbt::HistogramAnalysisState derived =
        hbt::analyze_histograms(config, raw);
    if (global.osl.bins != original_osl ||
        global.delta_t.bins != original_delta_t) {
        return fail("post-sample analysis modified completed raw counts");
    }

    const hbt::ShapeHistogramResult& shape =
        derived.products[0U].origins[0U].global.osl[0U];
    if (!shape.region.has_value() || shape.normalized_bins.size() != 20U ||
        !shape.gaussian.fully_valid || !shape.gaussian_neyman.fully_valid ||
        !shape.gaussian_pearson.fully_valid || !shape.mixed.fully_valid ||
        !shape.mixed_neyman.fully_valid || !shape.mixed_pearson.fully_valid ||
        shape.gaussian.estimator != hbt::FitEstimator::Poisson ||
        shape.gaussian_neyman.estimator != hbt::FitEstimator::Neyman ||
        shape.gaussian_pearson.estimator != hbt::FitEstimator::Pearson ||
        shape.mixed.estimator != hbt::FitEstimator::Poisson ||
        shape.mixed_neyman.estimator != hbt::FitEstimator::Neyman ||
        shape.mixed_pearson.estimator != hbt::FitEstimator::Pearson) {
        return fail("post-sample shape analysis did not produce valid state");
    }
    double normalized_integral = 0.0;
    for (const hbt::NormalizedHistogramBin& bin : shape.normalized_bins) {
        normalized_integral += bin.pdf * (bin.upper_edge - bin.lower_edge);
    }
    if (std::fabs(normalized_integral - 1.0) > 1.0e-12) {
        return fail("final normalized distribution does not integrate to one");
    }

    const hbt::ShapeHistogramResult& slice_osl =
        derived.products[0U].origins[0U].slices[0U].osl[0U];
    if (slice_osl.gaussian.failure_reason !=
            hbt::FitFailureReason::NotApplicable ||
        slice_osl.gaussian_neyman.failure_reason !=
            hbt::FitFailureReason::NotApplicable ||
        slice_osl.gaussian_pearson.failure_reason !=
            hbt::FitFailureReason::NotApplicable ||
        slice_osl.mixed.failure_reason !=
            hbt::FitFailureReason::NotApplicable ||
        slice_osl.mixed_neyman.failure_reason !=
            hbt::FitFailureReason::NotApplicable ||
        slice_osl.mixed_pearson.failure_reason !=
            hbt::FitFailureReason::NotApplicable ||
        slice_osl.gaussian.migrad.attempted ||
        slice_osl.gaussian_neyman.migrad.attempted ||
        slice_osl.gaussian_pearson.migrad.attempted ||
        slice_osl.mixed.starts_attempted != 0U ||
        slice_osl.mixed_neyman.starts_attempted != 0U ||
        slice_osl.mixed_pearson.starts_attempted != 0U) {
        return fail("OSL kinetic slice unexpectedly executed a fit");
    }

    const hbt::ShapeHistogramResult& low_stat_radial =
        derived.products[0U].origins[0U].slices[0U].radial[0U];
    if (!low_stat_radial.region.has_value() ||
        low_stat_radial.selected_count >= 10000U ||
        low_stat_radial.gaussian.failure_reason ==
            hbt::FitFailureReason::InsufficientStatistics ||
        low_stat_radial.gaussian_neyman.failure_reason ==
            hbt::FitFailureReason::InsufficientStatistics ||
        low_stat_radial.gaussian_pearson.failure_reason ==
            hbt::FitFailureReason::InsufficientStatistics ||
        low_stat_radial.mixed.failure_reason ==
            hbt::FitFailureReason::InsufficientStatistics ||
        low_stat_radial.mixed_neyman.failure_reason ==
            hbt::FitFailureReason::InsufficientStatistics ||
        low_stat_radial.mixed_pearson.failure_reason ==
            hbt::FitFailureReason::InsufficientStatistics ||
        low_stat_radial.gaussian.starts_attempted !=
            hbt::GaussianFitResult::kStartCount ||
        low_stat_radial.gaussian_neyman.starts_attempted !=
            hbt::GaussianFitResult::kStartCount ||
        low_stat_radial.gaussian_pearson.starts_attempted !=
            hbt::GaussianFitResult::kStartCount ||
        low_stat_radial.mixed.starts_attempted !=
            hbt::MixedFitResult::kCoreStartCount ||
        low_stat_radial.mixed_neyman.starts_attempted !=
            hbt::MixedFitResult::kCoreStartCount ||
        low_stat_radial.mixed_pearson.starts_attempted !=
            hbt::MixedFitResult::kCoreStartCount) {
        return fail("radial mT threshold zero did not allow the low-count fits");
    }

    const std::filesystem::path root =
        std::filesystem::temp_directory_path() /
        "hbt_production_output_integration_test";
    std::filesystem::remove_all(root);

    const config::RunConfig run_config{
        {},
        root,
        1U,
        1U,
        true,
        std::filesystem::path{"unused_hbt.yaml"}
    };
    const app::AnalysisStartupState startup{
        run_config,
        config,
        hbt::build_hbt_startup_state(config)
    };
    app::AnalysisRunSummary summary{
        startup,
        std::nullopt,
        std::nullopt,
        raw,
        derived
    };

    output::write_production_output(summary);

    const std::filesystem::path catalog_path =
        root / "product_catalog.csv";
    if (!std::filesystem::exists(catalog_path)) {
        std::filesystem::remove_all(root);
        return fail("configured product catalog was not serialized");
    }
    const std::string product_catalog = read_text(catalog_path);
    if (product_catalog.find(
            "0,pi_plus_pi_plus + pi_minus_pi_minus,"
            "pi_plus_pi_plus+pi_minus_pi_minus,0,pi_plus_pi_plus,"
            "pi_plus,pi_plus"
        ) == std::string::npos ||
        product_catalog.find(
            "0,pi_plus_pi_plus + pi_minus_pi_minus,"
            "pi_plus_pi_plus+pi_minus_pi_minus,1,pi_minus_pi_minus,"
            "pi_minus,pi_minus"
        ) == std::string::npos ||
        !csv_is_rectangular(catalog_path)) {
        std::filesystem::remove_all(root);
        return fail("product catalog lost configured or resolved identity");
    }

    const std::filesystem::path base =
        root / "product_0" / "global" / "P";
    const std::filesystem::path lcms_out =
        base / "LCMS" / "osl" / "r_out";
    const std::filesystem::path lcms_side =
        base / "LCMS" / "osl" / "r_side";
    if (!std::filesystem::exists(lcms_out / "distribution.csv") ||
        !std::filesystem::exists(lcms_out / "fit_parameters.csv") ||
        std::filesystem::exists(lcms_side / "distribution.csv") ||
        !std::filesystem::exists(lcms_side / "fit_parameters.csv")) {
        std::filesystem::remove_all(root);
        return fail("valid and invalid shape files violate output contract");
    }

    if (std::filesystem::exists(base / "PRF" / "osl" / "r_side") ||
        std::filesystem::exists(base / "PRF" / "osl" / "r_long") ||
        std::filesystem::exists(base / "LAB" / "osl") ||
        std::filesystem::exists(base / "LAB" / "radial") ||
        std::filesystem::exists(base / "shared")) {
        std::filesystem::remove_all(root);
        return fail("production tree duplicated or invented physical branches");
    }

    const std::string valid_distribution =
        read_text(lcms_out / "distribution.csv");
    if (valid_distribution.find("gaussian_fit_pdf") == std::string::npos ||
        valid_distribution.find("gaussian_fit_pdf_neyman") == std::string::npos ||
        valid_distribution.find("gaussian_fit_pdf_pearson") == std::string::npos ||
        valid_distribution.find("mixed_fit_pdf") == std::string::npos ||
        valid_distribution.find("mixed_fit_pdf_neyman") == std::string::npos ||
        valid_distribution.find("mixed_fit_pdf_pearson") == std::string::npos) {
        std::filesystem::remove_all(root);
        return fail("valid fit curves were not serialized with their PDF");
    }
    const std::filesystem::path lab_delta_t =
        base / "LAB" / "dt" / "delta_t";
    if (!std::filesystem::exists(lab_delta_t / "distribution.csv") ||
        !std::filesystem::exists(lab_delta_t / "statistics.csv") ||
        std::filesystem::exists(lab_delta_t / "fit_parameters.csv")) {
        std::filesystem::remove_all(root);
        return fail("delta-t production files violate the no-fit contract");
    }
    const std::string delta_statistics =
        read_text(lab_delta_t / "statistics.csv");
    if (delta_statistics.find("valid") == std::string::npos ||
        delta_statistics.find("sigma_dt") == std::string::npos ||
        delta_statistics.find("error_sigma_dt") == std::string::npos) {
        std::filesystem::remove_all(root);
        return fail("delta-t statistics were not serialized completely");
    }

    const std::string invalid_parameters =
        read_text(lcms_side / "fit_parameters.csv");
    if (invalid_parameters.find("empty_histogram") == std::string::npos ||
        invalid_parameters.find("R_core") == std::string::npos ||
        invalid_parameters.find("R_tail") == std::string::npos ||
        invalid_parameters.find("f_core") == std::string::npos) {
        std::filesystem::remove_all(root);
        return fail("invalid fit did not retain explicit model diagnostics");
    }
    const std::string valid_parameters =
        read_text(lcms_out / "fit_parameters.csv");
    if (valid_parameters.find("minos_lower_valid") == std::string::npos ||
        valid_parameters.find("error_method") == std::string::npos ||
        valid_parameters.find("hesse_attempted") == std::string::npos ||
        valid_parameters.find("hesse_valid_covariance") == std::string::npos ||
        valid_parameters.find("fit_upper_edge") == std::string::npos ||
        valid_parameters.find("R_G_core") == std::string::npos ||
        valid_parameters.find("estimator") == std::string::npos ||
        valid_parameters.find("poisson") == std::string::npos ||
        valid_parameters.find("neyman") == std::string::npos ||
        valid_parameters.find("pearson") == std::string::npos ||
        valid_parameters.find("consensus_size") == std::string::npos ||
        valid_parameters.find("core_start0_valid") == std::string::npos ||
        valid_parameters.find("core_start4_valid") == std::string::npos ||
        valid_parameters.find("core_start35_valid") == std::string::npos ||
        valid_parameters.find("tail_below_core") != std::string::npos ||
        valid_parameters.find("start_a_valid") != std::string::npos ||
        valid_parameters.find("start_b_valid") != std::string::npos) {
        std::filesystem::remove_all(root);
        return fail("new fit diagnostics are missing or legacy A/B fields remain");
    }
    if (!csv_is_rectangular(lcms_out / "fit_parameters.csv") ||
        !csv_is_rectangular(lcms_side / "fit_parameters.csv")) {
        std::filesystem::remove_all(root);
        return fail("fit parameter CSV columns are misaligned");
    }

    const std::filesystem::path slice0 =
        root / "product_0" / "mT_slice0_0.5-0.7" / "P";
    const std::filesystem::path slice1 =
        root / "product_0" / "mT_slice1_0.7-0.9" / "P";
    if (!std::filesystem::exists(
            slice0 / "LCMS" / "radial" / "r_radial" /
            "fit_parameters.csv") ||
        !std::filesystem::exists(
            slice1 / "PRF" / "radial" / "r_radial" /
            "fit_parameters.csv") ||
        std::filesystem::exists(slice0 / "LCMS" / "osl") ||
        std::filesystem::exists(root / "product_0" / "slice_0")) {
        std::filesystem::remove_all(root);
        return fail("configured slice ranges were not used as directory names");
    }

    const std::string low_stat_parameters = read_text(
        slice0 / "LCMS" / "radial" / "r_radial" / "fit_parameters.csv"
    );
    if (low_stat_parameters.find("insufficient_statistics") !=
        std::string::npos) {
        std::filesystem::remove_all(root);
        return fail("threshold-zero radial slice was still rejected by count");
    }

    for (const auto& entry :
         std::filesystem::recursive_directory_iterator(root)) {
        if (entry.is_directory() &&
            std::filesystem::directory_iterator(entry.path()) ==
                std::filesystem::directory_iterator{}) {
            std::filesystem::remove_all(root);
            return fail("production output created an empty directory");
        }
    }

    std::filesystem::remove_all(root);
    return true;
}

/**
 * @brief Verify kT/mT directory tokens follow arbitrary configured edges.
 * @return true when Cartesian slice names preserve both configured ranges.
 *
 * This case intentionally uses 0.25-GeV intervals to prove the serializer is
 * independent of the 0.20-GeV production mT choice. Empty histograms avoid
 * introducing additional fit behavior into this naming-only regression.
 */
bool verify_cartesian_slice_directory_names() {
    hbt::HBTConfig config = make_config();
    config.pair_slicing = {
        {true, {0.10, 0.35, 0.60}},
        {true, {0.50, 0.75, 1.00}}
    };

    const hbt::RawHistogramState raw =
        hbt::make_zero_raw_histogram_state(config);
    const hbt::HistogramAnalysisState derived =
        hbt::analyze_histograms(config, raw);

    const std::filesystem::path root =
        std::filesystem::temp_directory_path() /
        "hbt_production_output_cartesian_slice_test";
    std::filesystem::remove_all(root);

    const config::RunConfig run_config{
        {},
        root,
        1U,
        1U,
        true,
        std::filesystem::path{"unused_hbt.yaml"}
    };
    const app::AnalysisStartupState startup{
        run_config,
        config,
        hbt::build_hbt_startup_state(config)
    };
    app::AnalysisRunSummary summary{
        startup,
        std::nullopt,
        std::nullopt,
        raw,
        derived
    };
    output::write_production_output(summary);

    const std::filesystem::path expected =
        root / "product_0" /
        "kT_slice1_0.35-0.6__mT_slice1_0.75-1" / "P" /
        "LCMS" / "radial" / "r_radial" / "fit_parameters.csv";
    const std::filesystem::path forbidden_osl =
        root / "product_0" /
        "kT_slice0_0.1-0.35__mT_slice0_0.5-0.75" / "P" /
        "LCMS" / "osl";
    const bool valid = std::filesystem::exists(expected) &&
        !std::filesystem::exists(forbidden_osl) &&
        !std::filesystem::exists(root / "product_0" / "slice_3");
    std::filesystem::remove_all(root);
    if (!valid) {
        return fail(
            "Cartesian kT/mT slice directories ignored configured edge values"
        );
    }
    return true;
}

}  // namespace

/**
 * @brief Run the raw-count to analysis to production-output integration
 *        check.
 * @return EXIT_SUCCESS when every checked contract holds, otherwise
 *         EXIT_FAILURE.
 */
int main() {
    try {
        if (!verify_analysis_and_output() ||
            !verify_cartesian_slice_directory_names()) {
            return EXIT_FAILURE;
        }
    } catch (const std::exception& error) {
        std::cerr
            << "hbt_production_output_integration_test: unexpected exception: "
                  << error.what() << ".\n";
        return EXIT_FAILURE;
    }
    return EXIT_SUCCESS;
}
