/**
 * @file binned_models_test.cpp
 * @brief Focused tests for exact-bin post-sample models and likelihood.
 */

#include "hbt/fits/binned_models.h"

#include <cmath>
#include <cstddef>
#include <cstdint>
#include <cstdlib>
#include <iostream>
#include <stdexcept>
#include <string>
#include <vector>

namespace {

/**
 * @brief Report one failed exact-bin model condition.
 * @param message Description of the violated contract.
 * @return Always false.
 */
bool fail(const char* message) {
    std::cerr << "binned_models_test: " << message << ".\n";
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
 * @brief Verify OSL and radial physical forms are not conflated.
 * @return true when their exact integrals follow their distinct formulas.
 */
bool verify_distinct_physical_forms() {
    const double pi = std::acos(-1.0);
    const double osl = hbt::gaussian_component_integral(
        hbt::FitObservableFamily::OSL,
        0.0,
        2.0,
        1.0
    );
    const double expected_osl = std::sqrt(pi) * std::erf(1.0);
    const double radial = hbt::gaussian_component_integral(
        hbt::FitObservableFamily::Radial,
        0.0,
        2.0,
        1.0
    );
    const double expected_radial =
        2.0 * std::sqrt(pi) * std::erf(1.0) - 4.0 / std::exp(1.0);
    if (!close(osl, expected_osl) || !close(radial, expected_radial) ||
        close(osl, radial)) {
        return fail("OSL and radial Gaussian integrals are not distinct");
    }
    return true;
}

/**
 * @brief Verify OSL and radial exponential forms retain distinct geometry.
 * @return true when exact integrals match their separate analytic formulas.
 */
bool verify_distinct_exponential_forms() {
    const double osl = hbt::exponential_component_integral(
        hbt::FitObservableFamily::OSL,
        0.0,
        2.0,
        1.0
    );
    const double expected_osl = 1.0 - std::exp(-2.0);
    const double radial = hbt::exponential_component_integral(
        hbt::FitObservableFamily::Radial,
        0.0,
        2.0,
        1.0
    );
    const double expected_radial = 2.0 - 10.0 * std::exp(-2.0);
    if (!close(osl, expected_osl) || !close(radial, expected_radial) ||
        close(osl, radial)) {
        return fail("OSL and radial exponential integrals are not distinct");
    }
    return true;
}

/**
 * @brief Verify each component is normalized on the selected exact-bin range.
 * @return true when probabilities sum to one and use selected edges.
 */
bool verify_exact_bin_normalization() {
    const hbt::HistogramBinningConfig binning{6U, 0.0, 3.0, 2.0};
    const hbt::StatisticalRegion region{1U, 4U, 100U};
    const std::vector<double> gaussian = hbt::gaussian_bin_probabilities(
        hbt::FitObservableFamily::OSL,
        binning,
        region,
        0.8
    );
    const std::vector<double> mixed = hbt::mixed_bin_probabilities(
        hbt::FitObservableFamily::Radial,
        binning,
        region,
        0.7,
        1.4,
        0.35
    );
    double gaussian_sum = 0.0;
    double mixed_sum = 0.0;
    for (const double value : gaussian) {
        gaussian_sum += value;
    }
    for (const double value : mixed) {
        mixed_sum += value;
    }
    if (gaussian.size() != 4U || mixed.size() != 4U ||
        !close(gaussian_sum, 1.0) || !close(mixed_sum, 1.0)) {
        return fail("component probabilities are not normalized on region");
    }
    return true;
}

/**
 * @brief Verify mixed-model endpoint fractions remain mathematical endpoints.
 * @return true when f_core 0 and 1 produce the pure tail and core models.
 */
bool verify_mixed_endpoint_degeneracies() {
    const hbt::HistogramBinningConfig binning{5U, 0.0, 5.0, 1.0};
    const hbt::StatisticalRegion region{0U, 4U, 100U};
    const double core_radius = 0.9;
    const double tail_radius = 1.7;
    const std::vector<double> pure_core = hbt::gaussian_bin_probabilities(
        hbt::FitObservableFamily::OSL,
        binning,
        region,
        core_radius
    );
    const std::vector<double> endpoint_core = hbt::mixed_bin_probabilities(
        hbt::FitObservableFamily::OSL,
        binning,
        region,
        core_radius,
        tail_radius,
        1.0
    );
    const std::vector<double> endpoint_tail = hbt::mixed_bin_probabilities(
        hbt::FitObservableFamily::OSL,
        binning,
        region,
        core_radius,
        tail_radius,
        0.0
    );

    std::vector<double> pure_tail;
    pure_tail.reserve(5U);
    double tail_norm = 0.0;
    for (std::size_t bin = 0U; bin < 5U; ++bin) {
        const double integral = hbt::exponential_component_integral(
            hbt::FitObservableFamily::OSL,
            static_cast<double>(bin),
            static_cast<double>(bin + 1U),
            tail_radius
        );
        pure_tail.push_back(integral);
        tail_norm += integral;
    }
    for (double& value : pure_tail) {
        value /= tail_norm;
    }

    for (std::size_t bin = 0U; bin < 5U; ++bin) {
        if (!close(endpoint_core[bin], pure_core[bin]) ||
            !close(endpoint_tail[bin], pure_tail[bin])) {
            return fail("mixed endpoint fraction changed a pure component");
        }
    }
    if (std::string(hbt::fit_failure_reason_token(
            hbt::FitFailureReason::DegenerateCoreFraction)) !=
        "degenerate_core_fraction") {
        return fail("mixed endpoint degeneracy has no stable report token");
    }
    if (std::string(hbt::fit_estimator_token(hbt::FitEstimator::Poisson)) !=
            "poisson" ||
        std::string(hbt::fit_estimator_token(hbt::FitEstimator::Neyman)) !=
            "neyman" ||
        std::string(hbt::fit_estimator_token(hbt::FitEstimator::Pearson)) !=
            "pearson") {
        return fail("fit estimator has no stable serialization token");
    }
    return true;
}


/**
 * @brief Verify radial exact-bin evaluation remains stable for very large R.
 * @return true when unbounded positive radii retain finite normalized models.
 */
bool verify_large_radius_radial_stability() {
    const hbt::HistogramBinningConfig binning{20U, 0.0, 10.0, 2.0};
    const hbt::StatisticalRegion region{0U, 19U, 100U};
    constexpr double radius = 1.0e8;

    const std::vector<double> gaussian = hbt::gaussian_bin_probabilities(
        hbt::FitObservableFamily::Radial,
        binning,
        region,
        radius
    );
    const std::vector<double> mixed = hbt::mixed_bin_probabilities(
        hbt::FitObservableFamily::Radial,
        binning,
        region,
        radius,
        radius,
        0.5
    );

    double gaussian_sum = 0.0;
    double mixed_sum = 0.0;
    for (const double value : gaussian) {
        if (!std::isfinite(value) || value <= 0.0) {
            return fail("large-R radial Gaussian produced an invalid bin");
        }
        gaussian_sum += value;
    }
    for (const double value : mixed) {
        if (!std::isfinite(value) || value <= 0.0) {
            return fail("large-R radial mixture produced an invalid bin");
        }
        mixed_sum += value;
    }
    if (!close(gaussian_sum, 1.0) || !close(mixed_sum, 1.0)) {
        return fail("large-R radial probabilities are not normalized");
    }

    const double lower = 9.5;
    const double upper = 10.0;
    const double geometric_limit =
        (upper * upper * upper - lower * lower * lower) / 3.0;
    const double gaussian_integral = hbt::gaussian_component_integral(
        hbt::FitObservableFamily::Radial,
        lower,
        upper,
        radius
    );
    const double exponential_integral = hbt::exponential_component_integral(
        hbt::FitObservableFamily::Radial,
        lower,
        upper,
        radius
    );
    const double gaussian_relative = std::fabs(
        gaussian_integral - geometric_limit
    ) / geometric_limit;
    const double exponential_relative = std::fabs(
        exponential_integral - geometric_limit
    ) / geometric_limit;
    if (gaussian_relative > 1.0e-12 || exponential_relative > 1.0e-6) {
        return fail("large-R radial integral lost its geometric limit");
    }
    return true;
}

/**
 * @brief Verify zero observed bins remain in Poisson deviance.
 * @return true when n_i == 0 contributes exactly 2*mu_i.
 */
bool verify_zero_count_likelihood_term() {
    const std::vector<std::uint64_t> bins{0U, 10U};
    const hbt::StatisticalRegion region{0U, 1U, 10U};
    const std::vector<double> probabilities{0.25, 0.75};
    const double actual = hbt::binned_poisson_deviance(
        bins,
        0U,
        region,
        probabilities
    );
    const double expected = 5.0 + 2.0 * (
        7.5 - 10.0 + 10.0 * std::log(10.0 / 7.5)
    );
    if (!close(actual, expected)) {
        return fail("zero-count bin was omitted from Poisson deviance");
    }
    return true;
}

/**
 * @brief Verify Neyman omits zero-count bins and uses observed denominators.
 * @return true when the historical n_i > 0 weighting is reproduced exactly.
 */
bool verify_neyman_chi_square_terms() {
    const std::vector<std::uint64_t> bins{0U, 10U};
    const hbt::StatisticalRegion region{0U, 1U, 10U};
    const std::vector<double> probabilities{0.25, 0.75};
    const double actual = hbt::binned_neyman_chi_square(
        bins,
        0U,
        region,
        probabilities
    );
    const double expected = (10.0 - 7.5) * (10.0 - 7.5) / 10.0;
    if (!close(actual, expected)) {
        return fail("Neyman chi-square did not omit zero-count bins exactly");
    }
    return true;
}

/**
 * @brief Verify Pearson keeps zero-count bins with expected denominators.
 * @return true when every selected bin contributes the documented term.
 */
bool verify_pearson_chi_square_terms() {
    const std::vector<std::uint64_t> bins{0U, 10U};
    const hbt::StatisticalRegion region{0U, 1U, 10U};
    const std::vector<double> probabilities{0.25, 0.75};
    const double actual = hbt::binned_pearson_chi_square(
        bins,
        0U,
        region,
        probabilities
    );
    const double expected = 2.5 +
        (10.0 - 7.5) * (10.0 - 7.5) / 7.5;
    if (!close(actual, expected)) {
        return fail("Pearson chi-square did not retain zero-count bins");
    }
    return true;
}

/**
 * @brief Verify likelihood N is anchored to raw selected counts.
 * @return true when inconsistent selected_count is rejected explicitly.
 */
bool verify_raw_count_anchor() {
    const std::vector<std::uint64_t> bins{1U, 2U, 3U};
    const hbt::StatisticalRegion inconsistent{0U, 2U, 7U};
    try {
        (void)hbt::binned_poisson_deviance(
            bins,
            0U,
            inconsistent,
            {0.2, 0.3, 0.5}
        );
    } catch (const std::invalid_argument&) {
        return true;
    }
    return fail("likelihood accepted N different from raw selected counts");
}

/**
 * @brief Verify the likelihood rejects non-unit model probabilities.
 * @return true when a finite positive probability vector with sum != 1 fails.
 */
bool verify_likelihood_probability_normalization() {
    const std::vector<std::uint64_t> bins{1U, 2U, 3U};
    const hbt::StatisticalRegion region{0U, 2U, 6U};
    try {
        (void)hbt::binned_poisson_deviance(
            bins,
            0U,
            region,
            {0.2, 0.3, 0.6}
        );
    } catch (const std::invalid_argument&) {
        return true;
    }
    return fail("likelihood accepted probabilities whose sum differs from one");
}

/**
 * @brief Verify presentation curves reuse integrated bin probabilities.
 * @return true when density is probability divided only by exact bin width.
 */
bool verify_integrated_curve_density() {
    const hbt::HistogramBinningConfig binning{4U, 0.0, 2.0, 2.0};
    const std::vector<double> probabilities{0.1, 0.2, 0.3, 0.4};
    const std::vector<double> density =
        hbt::probabilities_to_pdf(probabilities, binning);
    if (density.size() != probabilities.size() ||
        !close(density[0U], 0.2) || !close(density[3U], 0.8)) {
        return fail("fit presentation density did not reuse bin integral");
    }
    return true;
}

}  // namespace

/**
 * @brief Run focused exact-bin model and likelihood regression checks.
 * @return EXIT_SUCCESS when every model contract holds, otherwise failure.
 */
int main() {
    if (!verify_distinct_physical_forms() ||
        !verify_distinct_exponential_forms() ||
        !verify_exact_bin_normalization() ||
        !verify_mixed_endpoint_degeneracies() ||
        !verify_large_radius_radial_stability() ||
        !verify_zero_count_likelihood_term() ||
        !verify_neyman_chi_square_terms() ||
        !verify_pearson_chi_square_terms() ||
        !verify_raw_count_anchor() ||
        !verify_likelihood_probability_normalization() ||
        !verify_integrated_curve_density()) {
        return EXIT_FAILURE;
    }
    return EXIT_SUCCESS;
}
