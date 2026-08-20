/**
 * @file minuit2_fitter_test.cpp
 * @brief Focused MIGRAD, multistart, and MINOS tests for post-sample analysis.
 */

#include "hbt/fits/binned_models.h"
#include "hbt/fits/minuit2_fitter.h"

#include <cmath>
#include <cstddef>
#include <cstdint>
#include <cstdlib>
#include <iostream>
#include <limits>
#include <vector>

namespace {

/**
 * @brief Report one failed Minuit2 fit condition.
 * @param message Description of the violated contract.
 * @return Always false.
 */
bool fail(const char* message) {
    std::cerr << "minuit2_fitter_test: " << message << ".\n";
    return false;
}

/**
 * @brief Convert model probabilities to deterministic positive raw counts.
 * @param probabilities Unit-normalized exact-bin probabilities.
 * @param total Nominal total count used for deterministic rounding.
 * @return Positive integer count in every selected bin.
 */
std::vector<std::uint64_t> make_counts(
    const std::vector<double>& probabilities,
    std::uint64_t total
) {
    std::vector<std::uint64_t> counts;
    counts.reserve(probabilities.size());
    for (const double probability : probabilities) {
        const double expected = probability * static_cast<double>(total);
        const auto rounded = static_cast<std::uint64_t>(
            std::llround(expected)
        );
        counts.push_back(rounded == 0U ? 1U : rounded);
    }
    return counts;
}

/**
 * @brief Return the exact uint64_t sum of one compact test count vector.
 * @param counts Raw test counts.
 * @return Exact sum; test magnitudes are intentionally far from overflow.
 */
std::uint64_t count_sum(const std::vector<std::uint64_t>& counts) {
    std::uint64_t sum = 0U;
    for (const std::uint64_t count : counts) {
        sum += count;
    }
    return sum;
}

/**
 * @brief Verify the pure Gaussian fit requires valid MIGRAD and MINOS.
 * @return true when a deterministic exact-model sample is fully valid.
 */
bool verify_gaussian_migrad_and_minos() {
    const hbt::HistogramBinningConfig binning{20U, 0.0, 10.0, 2.0};
    const hbt::StatisticalRegion model_region{0U, 19U, 1U};
    const std::vector<double> probabilities =
        hbt::gaussian_bin_probabilities(
            hbt::FitObservableFamily::OSL,
            binning,
            model_region,
            1.35
        );
    const std::vector<std::uint64_t> counts =
        make_counts(probabilities, 500000U);
    const hbt::StatisticalRegion region{
        0U,
        19U,
        count_sum(counts)
    };
    const hbt::GaussianFitResult fit = hbt::fit_gaussian_model(
        hbt::FitObservableFamily::OSL,
        counts,
        0U,
        binning,
        region
    );
    if (!fit.fully_valid || !fit.migrad.attempted ||
        !fit.migrad.valid || !fit.minos_radius.attempted ||
        !fit.minos_radius.lower_valid ||
        !fit.minos_radius.upper_valid || !fit.radius.has_value() ||
        !fit.q_min.has_value() || fit.fitted_pdf.size() != 20U) {
        return fail("exact Gaussian sample did not pass MIGRAD and MINOS");
    }
    if (fit.radius->value <= 0.0 || fit.radius->lower_error < 0.0 ||
        fit.radius->upper_error < 0.0) {
        return fail("Gaussian physical radius interval is invalid");
    }
    return true;
}

/**
 * @brief Verify Gaussian-core-anchored mixed consensus and selected MINOS.
 * @return true when five deterministic core starts reach a publishable basin.
 */
bool verify_mixed_multistart_and_minos() {
    const hbt::HistogramBinningConfig binning{20U, 0.0, 10.0, 2.0};
    const hbt::StatisticalRegion model_region{0U, 19U, 1U};
    const std::vector<double> probabilities = hbt::mixed_bin_probabilities(
        hbt::FitObservableFamily::OSL,
        binning,
        model_region,
        0.85,
        2.4,
        0.68
    );
    const std::vector<std::uint64_t> counts =
        make_counts(probabilities, 800000U);
    const hbt::StatisticalRegion region{
        0U,
        19U,
        count_sum(counts)
    };
    const hbt::GaussianFitResult gaussian = hbt::fit_gaussian_model(
        hbt::FitObservableFamily::OSL,
        counts,
        0U,
        binning,
        region
    );
    if (!gaussian.fully_valid) {
        return fail("mixed synthetic sample did not provide a valid Gaussian anchor");
    }

    const hbt::MixedFitResult mixed = hbt::fit_mixed_model(
        hbt::FitObservableFamily::OSL,
        counts,
        0U,
        binning,
        region,
        gaussian
    );
    if (mixed.starts_attempted != hbt::MixedFitResult::kCoreStartCount ||
        mixed.valid_starts < 4U || mixed.consensus_size < 4U ||
        !mixed.selected_core_start.has_value()) {
        return fail("mixed five-start core-basin consensus was not established");
    }
    for (const hbt::MigradDiagnostic& start : mixed.starts) {
        if (!start.attempted) {
            return fail("mixed contract did not attempt all five core starts");
        }
    }
    if (!mixed.fully_valid || !mixed.core_radius.has_value() ||
        !mixed.tail_radius.has_value() ||
        !mixed.core_fraction.has_value() ||
        !mixed.minos_core_radius.attempted ||
        !mixed.minos_tail_radius.attempted ||
        !mixed.minos_core_fraction.attempted ||
        !mixed.minos_core_radius.lower_valid ||
        !mixed.minos_core_radius.upper_valid ||
        !mixed.minos_tail_radius.lower_valid ||
        !mixed.minos_tail_radius.upper_valid ||
        !mixed.minos_core_fraction.lower_valid ||
        !mixed.minos_core_fraction.upper_valid) {
        return fail("selected mixed minimum did not pass required MINOS");
    }
    if (mixed.core_fraction->value <= 0.0 ||
        mixed.core_fraction->value >= 1.0) {
        return fail("mixed fit published a degenerate core fraction");
    }
    return true;
}


/**
 * @brief Verify deterministic classification of MIGRAD invalid states.
 * @return true when every required diagnostic maps to its explicit cause.
 */
bool verify_migrad_failure_classification() {
    const hbt::MigradDiagnostic valid{
        true,
        true,
        true,
        false,
        false,
        false,
        12,
        1.0
    };
    if (hbt::fit_failure_from_migrad(valid) !=
        hbt::FitFailureReason::None) {
        return fail("valid MIGRAD diagnostic was classified as invalid");
    }

    hbt::MigradDiagnostic diagnostic = valid;
    diagnostic.objective_failure = true;
    if (hbt::fit_failure_from_migrad(diagnostic) !=
        hbt::FitFailureReason::ObjectiveEvaluation) {
        return fail("objective failure was not classified explicitly");
    }
    diagnostic = valid;
    diagnostic.reached_call_limit = true;
    diagnostic.valid = false;
    if (hbt::fit_failure_from_migrad(diagnostic) !=
        hbt::FitFailureReason::MigradCallLimit) {
        return fail("MIGRAD call limit was not classified explicitly");
    }
    diagnostic = valid;
    diagnostic.above_max_edm = true;
    diagnostic.valid = false;
    if (hbt::fit_failure_from_migrad(diagnostic) !=
        hbt::FitFailureReason::MigradAboveMaxEdm) {
        return fail("MIGRAD EDM failure was not classified explicitly");
    }
    diagnostic = valid;
    diagnostic.q_min = std::nullopt;
    diagnostic.valid = false;
    if (hbt::fit_failure_from_migrad(diagnostic) !=
        hbt::FitFailureReason::NonFiniteMinimum) {
        return fail("non-finite MIGRAD minimum was not classified explicitly");
    }
    diagnostic = valid;
    diagnostic.valid_covariance = false;
    if (hbt::fit_failure_from_migrad(diagnostic) !=
        hbt::FitFailureReason::MigradInvalid) {
        return fail("invalid MIGRAD covariance was not rejected");
    }
    diagnostic = valid;
    diagnostic.valid = false;
    if (hbt::fit_failure_from_migrad(diagnostic) !=
        hbt::FitFailureReason::MigradInvalid) {
        return fail("generic invalid MIGRAD minimum lost its explicit state");
    }
    return true;
}

/**
 * @brief Verify deterministic classification of all required MINOS states.
 * @return true when each side failure and bounded-limit state is preserved.
 */
bool verify_minos_failure_classification() {
    const hbt::MinosDiagnostic valid{
        true,
        true,
        true,
        false,
        false,
        false,
        false,
        false,
        false
    };
    if (hbt::fit_failure_from_minos(valid, false) !=
        hbt::FitFailureReason::None) {
        return fail("valid MINOS diagnostic was classified as invalid");
    }

    struct Case {
        hbt::MinosDiagnostic diagnostic;
        bool reject_limits;
        hbt::FitFailureReason expected;
    };
    const Case cases[] = {
        {{true, true, true, false, false, false, false, true, false},
         false, hbt::FitFailureReason::MinosLowerNewMinimum},
        {{true, true, true, false, false, false, false, false, true},
         false, hbt::FitFailureReason::MinosUpperNewMinimum},
        {{true, true, true, false, false, true, false, false, false},
         false, hbt::FitFailureReason::MinosLowerCallLimit},
        {{true, true, true, false, false, false, true, false, false},
         false, hbt::FitFailureReason::MinosUpperCallLimit},
        {{true, true, true, true, false, false, false, false, false},
         true, hbt::FitFailureReason::MinosLowerLimit},
        {{true, true, true, false, true, false, false, false, false},
         true, hbt::FitFailureReason::MinosUpperLimit},
        {{true, false, true, false, false, false, false, false, false},
         false, hbt::FitFailureReason::MinosLowerInvalid},
        {{true, true, false, false, false, false, false, false, false},
         false, hbt::FitFailureReason::MinosUpperInvalid}
    };
    for (const Case& item : cases) {
        if (hbt::fit_failure_from_minos(
                item.diagnostic,
                item.reject_limits
            ) != item.expected) {
            return fail("MINOS invalid state lost its explicit cause");
        }
    }

    hbt::MinosDiagnostic radius_limit = valid;
    radius_limit.at_lower_limit = true;
    if (hbt::fit_failure_from_minos(radius_limit, false) !=
        hbt::FitFailureReason::None) {
        return fail("unbounded log-radius MINOS was given a physical limit");
    }
    return true;
}

/**
 * @brief Verify exact mixed-fraction endpoints are explicit degeneracies.
 * @return true when endpoint and out-of-domain states remain distinguishable.
 */
bool verify_core_fraction_classification() {
    if (hbt::mixed_core_fraction_failure(0.0) !=
            hbt::FitFailureReason::DegenerateCoreFraction ||
        hbt::mixed_core_fraction_failure(1.0) !=
            hbt::FitFailureReason::DegenerateCoreFraction ||
        hbt::mixed_core_fraction_failure(0.5) !=
            hbt::FitFailureReason::None ||
        hbt::mixed_core_fraction_failure(-0.1) !=
            hbt::FitFailureReason::NonFiniteMinimum ||
        hbt::mixed_core_fraction_failure(
            std::numeric_limits<double>::quiet_NaN()
        ) != hbt::FitFailureReason::NonFiniteMinimum) {
        return fail("mixed core-fraction states are not classified correctly");
    }
    return true;
}

/**
 * @brief Verify K >= P+1 is the deterministic pre-fit requirement.
 * @return true when underspecified fits are explicitly invalid.
 */
bool verify_insufficient_bins_are_reported() {
    const hbt::HistogramBinningConfig binning{3U, 0.0, 3.0, 1.0};
    const std::vector<std::uint64_t> counts{10U, 8U, 6U};
    const hbt::StatisticalRegion one_bin{0U, 0U, 10U};
    const hbt::GaussianFitResult gaussian = hbt::fit_gaussian_model(
        hbt::FitObservableFamily::OSL,
        counts,
        0U,
        binning,
        one_bin
    );
    if (gaussian.fully_valid ||
        gaussian.failure_reason != hbt::FitFailureReason::InsufficientBins ||
        gaussian.migrad.attempted || gaussian.minos_radius.attempted) {
        return fail("Gaussian K < P+1 was not explicitly rejected");
    }

    const hbt::StatisticalRegion three_bins{0U, 2U, 24U};
    const hbt::MixedFitResult mixed = hbt::fit_mixed_model(
        hbt::FitObservableFamily::OSL,
        counts,
        0U,
        binning,
        three_bins,
        gaussian
    );
    if (mixed.fully_valid ||
        mixed.failure_reason != hbt::FitFailureReason::InsufficientBins ||
        mixed.starts_attempted != 0U) {
        return fail("mixed K < P+1 was not explicitly rejected");
    }
    return true;
}

}  // namespace

/**
 * @brief Run focused MIGRAD, multistart, MINOS, and diagnostic-state checks.
 * @return EXIT_SUCCESS when every Minuit2 contract holds, otherwise failure.
 */
int main() {
    if (!verify_gaussian_migrad_and_minos() ||
        !verify_mixed_multistart_and_minos() ||
        !verify_migrad_failure_classification() ||
        !verify_minos_failure_classification() ||
        !verify_core_fraction_classification() ||
        !verify_insufficient_bins_are_reported()) {
        return EXIT_FAILURE;
    }
    return EXIT_SUCCESS;
}
