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
 * @brief Verify shape selection retains the full configured histogram range.
 * @return true when empty bins and later occupied islands are all retained.
 */
bool verify_shape_region() {
    const hbt::HistogramBinningConfig binning{8U, 0.0, 8.0, 1.0};
    const std::vector<std::uint64_t> bins{0U, 0U, 2U, 5U, 5U, 3U, 0U, 4U};
    const std::optional<hbt::StatisticalRegion> region =
        hbt::select_shape_region(bins, 0U, binning);
    if (!region.has_value()) {
        return fail("non-empty shape histogram had no selected region");
    }
    if (region->first_bin != 0U || region->last_bin != 7U ||
        region->selected_count != 19U) {
        return fail("shape region did not retain the full histogram range");
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
 * @brief Verify radial Gaussian-core selection uses only the right branch.
 * @return true when PAVA removes a recrossing and the threshold bin is excluded.
 */
bool verify_radial_gaussian_core_region() {
    const hbt::HistogramBinningConfig binning{10U, 0.0, 10.0, 1.0};
    const std::vector<std::uint64_t> bins{
        1U, 3U, 10U, 10U, 8U, 6U, 4U, 1U, 2U, 0U
    };
    const std::optional<hbt::StatisticalRegion> full =
        hbt::select_shape_region(bins, 0U, binning);
    if (!full.has_value() || full->last_bin != 9U) {
        return fail("radial core fixture lost its full statistical region");
    }
    const std::optional<hbt::StatisticalRegion> core =
        hbt::select_gaussian_core_region(
            hbt::FitObservableFamily::Radial,
            bins,
            0U,
            binning,
            full.value(),
            0.20
        );
    // PAVA pools [1,2] at bins 7,8 to 1.5. The 20% threshold is 2.0,
    // therefore bin 7 is the first excluded bin and the fit ends at bin 6.
    if (!core.has_value() || core->first_bin != 0U ||
        core->last_bin != 6U || core->selected_count != 42U) {
        return fail("radial Gaussian-core PAVA/threshold selection is wrong");
    }
    return true;
}

/**
 * @brief Verify OSL core selection starts at x=0 rather than an interior mode.
 * @return true when the monotonic reference is the PAVA level at bin zero.
 */
bool verify_osl_gaussian_core_region() {
    const hbt::HistogramBinningConfig binning{8U, 0.0, 8.0, 1.0};
    const std::vector<std::uint64_t> bins{9U, 10U, 8U, 6U, 3U, 1U, 2U, 0U};
    const std::optional<hbt::StatisticalRegion> full =
        hbt::select_shape_region(bins, 0U, binning);
    if (!full.has_value() || full->last_bin != 7U) {
        return fail("OSL core fixture lost its full statistical region");
    }
    const std::optional<hbt::StatisticalRegion> core =
        hbt::select_gaussian_core_region(
            hbt::FitObservableFamily::OSL,
            bins,
            0U,
            binning,
            full.value(),
            0.20
        );
    // PAVA pools the first two bins to 9.5 and the 1,2 recrossing to 1.5.
    // 20% of 9.5 is 1.9, so bin 5 is excluded and bin 4 is the last fit bin.
    if (!core.has_value() || core->first_bin != 0U ||
        core->last_bin != 4U || core->selected_count != 36U) {
        return fail("OSL Gaussian-core selection did not start from x=0");
    }
    return true;
}

/**
 * @brief Verify OSL half-maximum width is converted to Gaussian R.
 * @return true when the folded |x| crossing uses FWHM = 2*x_half.
 */
bool verify_osl_half_maximum_radius_seed() {
    const hbt::HistogramBinningConfig binning{5U, 0.0, 5.0, 1.0};
    const std::vector<std::uint64_t> bins{100U, 50U, 20U, 5U, 0U};
    const hbt::StatisticalRegion region{0U, 3U, 175U};
    const std::optional<double> seed = hbt::half_maximum_radius_seed(
        hbt::FitObservableFamily::OSL,
        bins,
        0U,
        binning,
        region
    );
    const double expected = 1.5 / (2.0 * std::sqrt(std::log(2.0)));
    if (!seed.has_value() || !close(seed.value(), expected)) {
        return fail("OSL half-maximum width was not converted to model R");
    }
    return true;
}

/**
 * @brief Verify radial FWHM is converted to the radial Gaussian R.
 * @return true when both measured half-height crossings are used.
 */
bool verify_radial_half_maximum_radius_seed() {
    const hbt::HistogramBinningConfig binning{6U, 0.0, 6.0, 1.0};
    const std::vector<std::uint64_t> bins{20U, 50U, 100U, 50U, 20U, 0U};
    const hbt::StatisticalRegion region{0U, 4U, 240U};
    const std::optional<double> seed = hbt::half_maximum_radius_seed(
        hbt::FitObservableFamily::Radial,
        bins,
        0U,
        binning,
        region
    );
    constexpr double radial_fwhm_over_radius = 2.3098847205021675;
    const double expected = 2.0 / radial_fwhm_over_radius;
    if (!seed.has_value() || !close(seed.value(), expected)) {
        return fail("radial half-maximum width was not converted to model R");
    }
    return true;
}


/**
 * @brief Verify OSL R_HM uses the non-increasing PAVA envelope.
 * @return true when an upward recrossing is pooled before locating half height.
 */
bool verify_osl_half_maximum_pava_recrossing() {
    const hbt::HistogramBinningConfig binning{6U, 0.0, 6.0, 1.0};
    const std::vector<std::uint64_t> bins{100U, 80U, 90U, 40U, 20U, 0U};
    const hbt::StatisticalRegion region{0U, 4U, 330U};
    const std::optional<double> seed = hbt::half_maximum_radius_seed(
        hbt::FitObservableFamily::OSL,
        bins,
        0U,
        binning,
        region
    );

    // PAVA pools 80 and 90 to 85. Half height is 50, so the crossing is
    // interpolated between centers 2.5 (85) and 3.5 (40).
    const double half_crossing = 2.5 + (35.0 / 45.0);
    const double expected =
        half_crossing / (2.0 * std::sqrt(std::log(2.0)));
    if (!seed.has_value() || !close(seed.value(), expected)) {
        return fail("OSL half-maximum seed ignored its PAVA envelope");
    }
    return true;
}

/**
 * @brief Verify radial R_HM ignores an isolated first-bin maximum.
 * @return true when unimodal PAVA preserves the broad interior physical peak.
 *
 * The first bin is deliberately one count above the broad peak. A raw-global-
 * maximum algorithm would place the radial mode at the left endpoint and fail
 * because no left half-height crossing exists. Least-squares unimodal PAVA
 * pools the isolated endpoint excess into the rising branch and selects the
 * interior peak instead.
 */
bool verify_radial_half_maximum_rejects_endpoint_spike() {
    const hbt::HistogramBinningConfig binning{11U, 0.0, 11.0, 1.0};
    const std::vector<std::uint64_t> bins{
        101U, 5U, 20U, 50U, 80U, 100U, 80U, 50U, 20U, 5U, 0U
    };
    const hbt::StatisticalRegion region{0U, 9U, 511U};
    const std::optional<double> seed = hbt::half_maximum_radius_seed(
        hbt::FitObservableFamily::Radial,
        bins,
        0U,
        binning,
        region
    );

    constexpr double radial_fwhm_over_radius = 2.3098847205021675;
    // The unimodal regression is [42,42,42,50,80,100,80,50,20,5].
    // Its half-height crossings are therefore 3.5 and 7.5.
    const double expected = 4.0 / radial_fwhm_over_radius;
    if (!seed.has_value() || !close(seed.value(), expected)) {
        return fail("radial half-maximum seed followed an endpoint spike");
    }
    return true;
}

/**
 * @brief Verify the core selector falls back to the full safety region.
 * @return true when no threshold crossing invents an earlier cut.
 */
bool verify_gaussian_core_fallback() {
    const hbt::HistogramBinningConfig binning{5U, 0.0, 5.0, 1.0};
    const std::vector<std::uint64_t> bins{10U, 9U, 8U, 7U, 6U};
    const hbt::StatisticalRegion full{0U, 4U, 40U};
    const auto core = hbt::select_gaussian_core_region(
        hbt::FitObservableFamily::OSL,
        bins,
        0U,
        binning,
        full,
        0.10
    );
    if (!core.has_value() || core->last_bin != 4U ||
        core->selected_count != 40U) {
        return fail("Gaussian-core selector failed its full-region fallback");
    }
    return true;
}

/**
 * @brief Verify numerical same-basin equivalence and connected grouping.
 * @return true when four nearby endpoints form one diagnostic component.
 *
 * Q is intentionally absent: this helper diagnoses repeated convergence only
 * and no longer decides which mixed minimum is accepted.
 */
bool verify_mixed_basin_grouping() {
    const std::vector<hbt::MixedBasinPoint> endpoints{
        {0.700, 2.000, 0.700},
        {0.704, 2.006, 0.696},
        {0.695, 1.995, 0.705},
        {0.702, 2.003, 0.701},
        {2.200, 0.500, 0.080}
    };
    if (!hbt::same_mixed_basin(endpoints[0U], endpoints[1U]) ||
        hbt::same_mixed_basin(endpoints[0U], endpoints[4U])) {
        return fail("mixed basin numerical equivalence tolerance is wrong");
    }

    const std::vector<std::size_t> valid_indices{0U, 1U, 2U, 3U, 4U};
    const std::vector<std::size_t> consensus =
        hbt::largest_mixed_basin_group(endpoints, valid_indices);
    if (consensus.size() != 4U) {
        return fail("mixed basin grouping did not preserve the four-point component");
    }
    for (const std::size_t index : consensus) {
        if (index >= 4U) {
            return fail("isolated alternative endpoint entered the grouped basin");
        }
    }
    return true;
}

/**
 * @brief Verify every origin prefers the most populated admissible mixed basin.
 * @return true when a three-start basin wins over a closer-to-R_HM two-start
 *         basin and the lowest q inside the selected basin is returned.
 */
bool verify_largest_mixed_basin_selection() {
    const std::vector<hbt::MixedBasinPoint> endpoints{
        // Two-start basin around the historical R_HM anchor.
        {std::log(4.00), std::log(6.00), 0.650},
        {std::log(4.02), std::log(6.03), 0.646},
        // Three-start remote basin.
        {std::log(14.00), std::log(3.50), 0.280},
        {std::log(14.04), std::log(3.52), 0.284},
        {std::log(13.96), std::log(3.48), 0.277}
    };
    const std::vector<double> q_values{20.0, 10.0, 90.0, 80.0, 70.0};
    const std::vector<std::size_t> valid_indices{0U, 1U, 2U, 3U, 4U};

    for (const auto policy : {
            hbt::MixedCoreFractionPolicy::RequireCoreAndTail,
            hbt::MixedCoreFractionPolicy::RejectPureGaussian}) {
        const std::optional<std::size_t> selected =
            hbt::select_mixed_start_by_largest_basin(
                endpoints,
                q_values,
                valid_indices,
                policy
            );
        if (!selected.has_value() || selected.value() != 4U) {
            return fail(
                "mixed basin selection did not prefer the largest admissible basin"
            );
        }
        if (!(q_values[1U] < q_values[selected.value()])) {
            return fail(
                "largest-basin test lacks the intended lower-q smaller basin"
            );
        }
        const std::vector<std::size_t> ranked =
            hbt::rank_mixed_starts_in_selected_basin(
                endpoints,
                q_values,
                valid_indices,
                policy
            );
        const std::vector<std::size_t> expected{4U, 3U, 2U};
        if (ranked != expected) {
            return fail(
                "selected mixed basin endpoints were not ranked by q"
            );
        }
    }
    return true;
}

/**
 * @brief Verify degenerate-f basins are removed before size ranking.
 * @return true when a realistic non-degenerate basin is selected even though
 *         degenerate basins may have a lower q.
 */
bool verify_degenerate_fraction_basin_filter() {
    const std::vector<hbt::MixedBasinPoint> endpoints{
        {std::log(10.40), std::log(3.03), 0.055},
        {std::log(10.44), std::log(3.04), 0.058},
        {std::log(2.83), std::log(3.81), 0.326},
        {std::log(2.85), std::log(3.79), 0.329},
        {std::log(3.92), std::log(3.19), 1.0e-9},
        {std::log(3.95), std::log(3.18), 8.0e-9}
    };
    const std::vector<double> q_values{
        653.0, 654.0, 684.0, 683.0, 500.0, 499.0
    };
    const std::vector<std::size_t> valid_indices{0U, 1U, 2U, 3U, 4U, 5U};

    const std::optional<std::size_t> selected =
        hbt::select_mixed_start_by_largest_basin(
            endpoints,
            q_values,
            valid_indices,
            hbt::MixedCoreFractionPolicy::RequireCoreAndTail
        );
    if (!selected.has_value() || selected.value() != 3U) {
        return fail(
            "degenerate-f basin filter did not select the physical mixed basin"
        );
    }
    return true;
}

/**
 * @brief Verify the PRD physical f_core basin bounds remain strict.
 * @return true when PRD basins at 0.1 or 0.99 yield no selectable solution.
 */
bool verify_no_nondegenerate_fraction_basin() {
    const std::vector<hbt::MixedBasinPoint> endpoints{
        {std::log(3.0), std::log(4.0), 0.100},
        {std::log(3.01), std::log(4.01), 0.100},
        {std::log(3.5), std::log(4.5), 0.990},
        {std::log(3.51), std::log(4.51), 0.990}
    };
    const std::vector<double> q_values{10.0, 9.0, 8.0, 7.0};
    const std::vector<std::size_t> valid_indices{0U, 1U, 2U, 3U};

    const std::optional<std::size_t> selected =
        hbt::select_mixed_start_by_largest_basin(
            endpoints,
            q_values,
            valid_indices,
            hbt::MixedCoreFractionPolicy::RequireCoreAndTail
        );
    if (selected.has_value()) {
        return fail("strict PRD f_core basin bounds admitted a degenerate basin");
    }
    return true;
}

/**
 * @brief Verify the common mixed-fraction bounds reject degeneracies.
 * @return true when PRD and P/PR reject f_core >= 0.99 but admit 0.98, while
 *         the lower 0.1 bound remains exclusive.
 */
bool verify_origin_specific_fraction_basin_policy() {
    const std::vector<hbt::MixedBasinPoint> endpoints{
        {std::log(3.00), std::log(4.00), 1.000},
        {std::log(6.00), std::log(3.00), 0.500}
    };
    const std::vector<double> q_values{20.0, 10.0};
    const std::vector<std::size_t> valid_indices{0U, 1U};

    const auto prd_selected = hbt::select_mixed_start_by_largest_basin(
        endpoints,
        q_values,
        valid_indices,
        hbt::MixedCoreFractionPolicy::RequireCoreAndTail
    );
    const auto p_pr_selected = hbt::select_mixed_start_by_largest_basin(
        endpoints,
        q_values,
        valid_indices,
        hbt::MixedCoreFractionPolicy::RejectPureGaussian
    );
    if (!prd_selected.has_value() || prd_selected.value() != 1U ||
        !p_pr_selected.has_value() || p_pr_selected.value() != 1U) {
        return fail("origin policy did not reject the pure-Gaussian basin");
    }

    const std::vector<hbt::MixedBasinPoint> p_pr_only{
        {std::log(3.0), std::log(4.0), 0.980}
    };
    const std::vector<double> single_q{1.0};
    const std::vector<std::size_t> single_index{0U};
    if (!hbt::select_mixed_start_by_largest_basin(
            p_pr_only,
            single_q,
            single_index,
            hbt::MixedCoreFractionPolicy::RequireCoreAndTail
        ).has_value() ||
        !hbt::select_mixed_start_by_largest_basin(
            p_pr_only,
            single_q,
            single_index,
            hbt::MixedCoreFractionPolicy::RejectPureGaussian
        ).has_value()) {
        return fail("common upper f_core bound did not admit f_core=0.98");
    }

    const std::vector<hbt::MixedBasinPoint> upper_boundary{
        {std::log(3.0), std::log(4.0), 0.990}
    };
    if (hbt::select_mixed_start_by_largest_basin(
            upper_boundary,
            single_q,
            single_index,
            hbt::MixedCoreFractionPolicy::RequireCoreAndTail
        ).has_value() ||
        hbt::select_mixed_start_by_largest_basin(
            upper_boundary,
            single_q,
            single_index,
            hbt::MixedCoreFractionPolicy::RejectPureGaussian
        ).has_value()) {
        return fail("policy admitted the exclusive f_core=0.99 boundary");
    }

    const std::vector<hbt::MixedBasinPoint> lower_boundary{
        {std::log(3.0), std::log(4.0), 0.100}
    };
    if (hbt::select_mixed_start_by_largest_basin(
            lower_boundary,
            single_q,
            single_index,
            hbt::MixedCoreFractionPolicy::RequireCoreAndTail
        ).has_value() ||
        hbt::select_mixed_start_by_largest_basin(
            lower_boundary,
            single_q,
            single_index,
            hbt::MixedCoreFractionPolicy::RejectPureGaussian
        ).has_value()) {
        return fail(
            "origin-specific policy admitted the exclusive f_core=0.1 boundary"
        );
    }
    return true;
}

/**
 * @brief Verify q breaks ties between equally populated admissible basins.
 * @return true when both origin policies select the lower-q equal-size basin
 *         and then the lower-q endpoint inside that basin.
 */
bool verify_equal_size_basin_q_tie_break() {
    const std::vector<hbt::MixedBasinPoint> endpoints{
        // First two-start basin.
        {std::log(4.00), std::log(7.00), 0.650},
        {std::log(4.02), std::log(7.03), 0.646},
        // Second two-start basin with smaller basin q.
        {std::log(8.00), std::log(3.50), 0.500},
        {std::log(8.02), std::log(3.52), 0.504}
    };
    const std::vector<double> q_values{120.0, 110.0, 90.0, 80.0};
    const std::vector<std::size_t> valid_indices{0U, 1U, 2U, 3U};

    for (const auto policy : {
            hbt::MixedCoreFractionPolicy::RequireCoreAndTail,
            hbt::MixedCoreFractionPolicy::RejectPureGaussian}) {
        const auto selected = hbt::select_mixed_start_by_largest_basin(
            endpoints,
            q_values,
            valid_indices,
            policy
        );
        if (!selected.has_value() || selected.value() != 3U) {
            return fail(
                "equal-size mixed basins were not ranked by their smallest q"
            );
        }
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
 * @brief Verify signed delta-t selection retains the full configured range.
 * @return true when internal zeros and external occupied islands are retained.
 */
bool verify_delta_t_region() {
    const hbt::HistogramBinningConfig binning{8U, -4.0, 4.0, 1.0};
    const std::vector<std::uint64_t> bins{4U, 0U, 2U, 6U, 6U, 3U, 0U, 5U};
    const std::optional<hbt::StatisticalRegion> region =
        hbt::select_delta_t_region(bins, 0U, binning);
    if (!region.has_value()) {
        return fail("non-empty delta-t histogram had no selected region");
    }
    if (region->first_bin != 0U || region->last_bin != 7U ||
        region->selected_count != 26U) {
        return fail("delta-t region did not retain the full histogram range");
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
 * @return true when mean, mean error, sigma, and sigma error are exact.
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
    const double mean_error = sigma / std::sqrt(17.0);
    const double sigma_error = sigma / std::sqrt(32.0);
    if (result.status != hbt::DeltaTStatisticsStatus::Valid ||
        !result.mean.has_value() || !result.mean_error.has_value() ||
        !result.sigma.has_value() || !result.sigma_error.has_value() ||
        !close(result.mean.value(), mean) ||
        !close(result.mean_error.value(), mean_error) ||
        !close(result.sigma.value(), sigma) ||
        !close(result.sigma_error.value(), sigma_error)) {
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
        result.mean_error.has_value() || result.sigma.has_value() ||
        result.sigma_error.has_value()) {
        return fail("invalid delta-t variance was hidden or clamped");
    }
    return true;
}

/**
 * @brief Verify N <= 1 is reported instead of inventing uncertainties.
 * @return true when the explicit insufficient-count state is preserved.
 */
bool verify_insufficient_delta_t_count() {
    const hbt::HistogramBinningConfig binning{3U, -1.5, 1.5, 1.0};
    const std::vector<std::uint64_t> bins{0U, 1U, 0U};
    const hbt::StatisticalRegion region{1U, 1U, 1U};
    const hbt::DeltaTHistogramResult result =
        hbt::calculate_delta_t_statistics(bins, 0U, binning, region);
    if (result.status != hbt::DeltaTStatisticsStatus::InsufficientCount ||
        !result.mean.has_value() || result.mean_error.has_value() ||
        !result.sigma.has_value() || result.sigma_error.has_value()) {
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
        !verify_radial_gaussian_core_region() ||
        !verify_osl_gaussian_core_region() ||
        !verify_osl_half_maximum_radius_seed() ||
        !verify_radial_half_maximum_radius_seed() ||
        !verify_osl_half_maximum_pava_recrossing() ||
        !verify_radial_half_maximum_rejects_endpoint_spike() ||
        !verify_gaussian_core_fallback() ||
        !verify_mixed_basin_grouping() ||
        !verify_largest_mixed_basin_selection() ||
        !verify_degenerate_fraction_basin_filter() ||
        !verify_no_nondegenerate_fraction_basin() ||
        !verify_origin_specific_fraction_basin_policy() ||
        !verify_equal_size_basin_q_tie_break() ||
        !verify_empty_regions() || !verify_delta_t_region() ||
        !verify_normalization() || !verify_delta_t_statistics() ||
        !verify_invalid_delta_t_variance() ||
        !verify_insufficient_delta_t_count()) {
        return EXIT_FAILURE;
    }
    return EXIT_SUCCESS;
}
