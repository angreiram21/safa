/**
 * @file statistical_analysis_test.cpp
 * @brief Focused tests for post-sample regions, normalization, and delta-t.
 */

#include "hbt/fits/statistical_analysis.h"

#include <cmath>
#include <cstddef>
#include <cstdint>
#include <cstdlib>
#include <iostream>
#include <optional>
#include <vector>

namespace {

/**
 * @brief Report one failed post-sample statistical-analysis condition.
 * @param message Description of the violated contract.
 * @return Always false.
 */
bool fail(const char* message) {
    std::cerr << "statistical_analysis_test: " << message << ".\n";
    return false;
}

/**
 * @brief Compare two finite doubles with a test-only relative tolerance.
 * @param actual Calculated value.
 * @param expected Reference value.
 * @return true when the values agree for the unit test.
 */
bool close(double actual, double expected) {
    const double scale = 1.0 + std::fabs(expected);
    return std::fabs(actual - expected) <= 1.0e-12 * scale;
}

/**
 * @brief Verify shape selection keeps leading zeros and cuts later islands.
 * @return true when the exact contiguous-region contract is satisfied.
 */
bool verify_shape_region() {
    const hbt::HistogramBinningConfig binning{8U, 0.0, 8.0, 1.0};
    const std::vector<std::uint64_t> bins{0U, 0U, 2U, 5U, 5U, 3U, 0U, 4U};
    const std::optional<hbt::StatisticalRegion> region =
        hbt::select_shape_region(bins, 0U, binning);
    if (!region.has_value()) {
        return fail("non-empty shape histogram had no selected region");
    }
    if (region->first_bin != 0U || region->last_bin != 5U ||
        region->selected_count != 15U) {
        return fail("shape region did not preserve origin-side empty bins");
    }
    return true;
}

/**
 * @brief Verify shape selection uses the full range when no tail zero exists.
 * @return true when no artificial right-tail cut is introduced.
 */
bool verify_shape_region_without_tail_cut() {
    const hbt::HistogramBinningConfig binning{6U, 0.0, 6.0, 1.0};
    const std::vector<std::uint64_t> bins{0U, 1U, 4U, 3U, 2U, 1U};
    const std::optional<hbt::StatisticalRegion> region =
        hbt::select_shape_region(bins, 0U, binning);
    if (!region.has_value() || region->first_bin != 0U ||
        region->last_bin != 5U || region->selected_count != 11U) {
        return fail("shape region invented a tail cut without an empty bin");
    }
    return true;
}

/**
 * @brief Verify empty shape and delta-t histograms have no selected region.
 * @return true when an all-zero histogram remains explicitly regionless.
 */
bool verify_empty_regions() {
    const hbt::HistogramBinningConfig shape_binning{4U, 0.0, 4.0, 1.0};
    const hbt::HistogramBinningConfig time_binning{4U, -2.0, 2.0, 1.0};
    const std::vector<std::uint64_t> bins(4U, 0U);
    if (hbt::select_shape_region(bins, 0U, shape_binning).has_value() ||
        hbt::select_delta_t_region(bins, 0U, time_binning).has_value()) {
        return fail("empty histogram acquired an artificial selected region");
    }
    return true;
}

/**
 * @brief Verify signed delta-t selection stops at first empty bin per side.
 * @return true when both delimiters and external islands are excluded.
 */
bool verify_delta_t_region() {
    const hbt::HistogramBinningConfig binning{8U, -4.0, 4.0, 1.0};
    const std::vector<std::uint64_t> bins{4U, 0U, 2U, 6U, 6U, 3U, 0U, 5U};
    const std::optional<hbt::StatisticalRegion> region =
        hbt::select_delta_t_region(bins, 0U, binning);
    if (!region.has_value()) {
        return fail("non-empty delta-t histogram had no selected region");
    }
    if (region->first_bin != 2U || region->last_bin != 5U ||
        region->selected_count != 17U) {
        return fail("delta-t bilateral empty-bin cut is incorrect");
    }
    return true;
}

/**
 * @brief Verify final normalization includes selected zero-count bins.
 * @return true when PDF and counting uncertainty follow the closed formulas.
 */
bool verify_normalization() {
    const hbt::HistogramBinningConfig binning{8U, 0.0, 8.0, 1.0};
    const std::vector<std::uint64_t> bins{0U, 0U, 2U, 5U, 5U, 3U, 0U, 4U};
    const hbt::StatisticalRegion region{0U, 5U, 15U};
    const std::vector<hbt::NormalizedHistogramBin> normalized =
        hbt::normalize_region(bins, 0U, binning, region);
    if (normalized.size() != 6U || normalized[0U].pdf != 0.0 ||
        normalized[0U].counting_error_pdf != 0.0) {
        return fail("selected empty bins were removed or assigned fake errors");
    }
    double integral = 0.0;
    for (const hbt::NormalizedHistogramBin& bin : normalized) {
        integral += bin.pdf * (bin.upper_edge - bin.lower_edge);
    }
    if (!close(integral, 1.0) ||
        !close(normalized[3U].pdf, 5.0 / 15.0) ||
        !close(normalized[3U].counting_error_pdf,
               std::sqrt(5.0) / 15.0)) {
        return fail("published PDF or dPDF formula is incorrect");
    }
    return true;
}

/**
 * @brief Verify delta-t moments use raw counts and selected bin centers.
 * @return true when mean, sigma, and required sigma error are exact.
 */
bool verify_delta_t_statistics() {
    const hbt::HistogramBinningConfig binning{8U, -4.0, 4.0, 1.0};
    const std::vector<std::uint64_t> bins{4U, 0U, 2U, 6U, 6U, 3U, 0U, 5U};
    const hbt::StatisticalRegion region{2U, 5U, 17U};
    const hbt::DeltaTHistogramResult result =
        hbt::calculate_delta_t_statistics(bins, 0U, binning, region);

    const double mean = 1.5 / 17.0;
    const double second_moment = 14.25 / 17.0;
    const double sigma = std::sqrt(second_moment - mean * mean);
    const double error = sigma / std::sqrt(32.0);
    if (result.status != hbt::DeltaTStatisticsStatus::Valid ||
        !result.mean.has_value() || !result.sigma.has_value() ||
        !result.sigma_error.has_value() ||
        !close(result.mean.value(), mean) ||
        !close(result.sigma.value(), sigma) ||
        !close(result.sigma_error.value(), error)) {
        return fail("delta-t raw-count moments differ from the contract");
    }
    return true;
}

/**
 * @brief Verify a non-finite delta-t variance is diagnosed, never clamped.
 * @return true when arithmetic overflow yields InvalidVariance explicitly.
 */
bool verify_invalid_delta_t_variance() {
    const hbt::HistogramBinningConfig binning{
        2U,
        -1.0e154,
        1.0e154,
        1.0e154
    };
    const std::vector<std::uint64_t> bins{10U, 10U};
    const hbt::StatisticalRegion region{0U, 1U, 20U};
    const hbt::DeltaTHistogramResult result =
        hbt::calculate_delta_t_statistics(bins, 0U, binning, region);
    if (result.status != hbt::DeltaTStatisticsStatus::InvalidVariance ||
        result.sigma.has_value() || result.sigma_error.has_value()) {
        return fail("invalid delta-t variance was hidden or clamped");
    }
    return true;
}

/**
 * @brief Verify N <= 1 is reported instead of inventing a sigma error.
 * @return true when the explicit insufficient-count state is preserved.
 */
bool verify_insufficient_delta_t_count() {
    const hbt::HistogramBinningConfig binning{3U, -1.5, 1.5, 1.0};
    const std::vector<std::uint64_t> bins{0U, 1U, 0U};
    const hbt::StatisticalRegion region{1U, 1U, 1U};
    const hbt::DeltaTHistogramResult result =
        hbt::calculate_delta_t_statistics(bins, 0U, binning, region);
    if (result.status != hbt::DeltaTStatisticsStatus::InsufficientCount ||
        !result.mean.has_value() || !result.sigma.has_value() ||
        result.sigma_error.has_value()) {
        return fail("N <= 1 did not remain an explicit invalid error state");
    }
    return true;
}

}  // namespace

/**
 * @brief Run focused statistical-region, normalization, and delta-t checks.
 * @return EXIT_SUCCESS when every checked contract holds, otherwise
 *         EXIT_FAILURE.
 */
int main() {
    if (!verify_shape_region() || !verify_shape_region_without_tail_cut() ||
        !verify_empty_regions() || !verify_delta_t_region() ||
        !verify_normalization() || !verify_delta_t_statistics() ||
        !verify_invalid_delta_t_variance() ||
        !verify_insufficient_delta_t_count()) {
        return EXIT_FAILURE;
    }
    return EXIT_SUCCESS;
}
