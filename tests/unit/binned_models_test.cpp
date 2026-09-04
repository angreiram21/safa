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
 * @brief Verify normalized Gaussian and unnormalized mixed exact-bin models.
 * @return true when each API retains its documented normalization semantics.
 */
bool verify_exact_bin_model_semantics() {
    const hbt::HistogramBinningConfig binning{6U, 0.0, 3.0, 2.0};
    const hbt::StatisticalRegion region{1U, 4U, 100U};
    const std::vector<double> gaussian = hbt::gaussian_bin_probabilities(
        hbt::FitObservableFamily::OSL,
        binning,
        region,
        0.8
    );
    const double core_radius = 0.7;
    const double tail_radius = 1.4;
    const double f = 0.35;
    const std::vector<double> mixed = hbt::mixed_bin_integrals(
        hbt::FitObservableFamily::Radial,
        binning,
        region,
        core_radius,
        tail_radius,
        f
    );
    double gaussian_sum = 0.0;
    for (const double value : gaussian) {
        gaussian_sum += value;
    }
    if (gaussian.size() != 4U || mixed.size() != 4U ||
        !close(gaussian_sum, 1.0)) {
        return fail("exact-bin model cardinality or Gaussian normalization failed");
    }
    const double lower = 0.5;
    const double upper = 1.0;
    const double expected =
        f * hbt::gaussian_component_integral(
            hbt::FitObservableFamily::Radial, lower, upper, core_radius
        ) +
        (1.0 - f) * hbt::exponential_component_integral(
            hbt::FitObservableFamily::Radial, lower, upper, tail_radius
        );
    if (!close(mixed.front(), expected)) {
        return fail("mixed model normalized components instead of mixing integrals");
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
    std::vector<double> pure_core;
    pure_core.reserve(5U);
    for (std::size_t bin = 0U; bin < 5U; ++bin) {
        pure_core.push_back(hbt::gaussian_component_integral(
            hbt::FitObservableFamily::OSL,
            static_cast<double>(bin),
            static_cast<double>(bin + 1U),
            core_radius
        ));
    }
    const std::vector<double> endpoint_core = hbt::mixed_bin_integrals(
        hbt::FitObservableFamily::OSL,
        binning,
        region,
        core_radius,
        tail_radius,
        1.0
    );
    const std::vector<double> endpoint_tail = hbt::mixed_bin_integrals(
        hbt::FitObservableFamily::OSL,
        binning,
        region,
        core_radius,
        tail_radius,
        0.0
    );

    std::vector<double> pure_tail;
    pure_tail.reserve(5U);
    for (std::size_t bin = 0U; bin < 5U; ++bin) {
        const double integral = hbt::exponential_component_integral(
            hbt::FitObservableFamily::OSL,
            static_cast<double>(bin),
            static_cast<double>(bin + 1U),
            tail_radius
        );
        pure_tail.push_back(integral);
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
    const std::vector<double> mixed = hbt::mixed_bin_integrals(
        hbt::FitObservableFamily::Radial,
        binning,
        region,
        radius,
        radius,
        0.5
    );

    double gaussian_sum = 0.0;
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
    }
    if (!close(gaussian_sum, 1.0)) {
        return fail("large-R radial Gaussian probabilities are not normalized");
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
 * @brief Verify analytic mixed A and the exact Neyman p_i == 0 limit.
 * @return true when the fixed-shape amplitude is exact and a populated zero
 *         model bin contributes n_i instead of invalidating the objective.
 */
bool verify_profiled_mixed_neyman_amplitude() {
    const std::vector<std::uint64_t> bins{10U, 5U, 2U};
    const hbt::StatisticalRegion region{0U, 2U, 17U};
    const std::vector<double> p{0.2, 0.1, 0.0};
    const double amplitude = hbt::neyman_optimal_mixed_amplitude(
        bins, 0U, region, p
    );
    const double expected_amplitude = 50.0 / 17.0;
    if (!close(amplitude, expected_amplitude)) {
        return fail("analytic mixed Neyman amplitude is incorrect");
    }
    const double chi_square = hbt::binned_neyman_chi_square_from_integrals(
        bins, 0U, region, p, amplitude
    );
    if (!close(chi_square, 2.0)) {
        return fail("populated p_i == 0 bin did not contribute exactly n_i");
    }
    return true;
}

/**
 * @brief Verify mixed presentation density includes fitted A and permits zero.
 */
bool verify_mixed_fitted_density() {
    const hbt::HistogramBinningConfig binning{3U, 0.0, 1.5, 2.0};
    const std::vector<double> p{0.2, 0.1, 0.0};
    const std::vector<double> density = hbt::mixed_integrals_to_pdf(
        p, 2.0, binning
    );
    if (density.size() != 3U || !close(density[0U], 0.8) ||
        !close(density[1U], 0.4) || !close(density[2U], 0.0)) {
        return fail("mixed fitted density lost A scaling or zero bins");
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
 * @brief Verify the optional model normalization rescales expected counts.
 * @return true when free-amplitude objective scaling is applied exactly.
 */
bool verify_free_normalization_multiplier() {
    const std::vector<std::uint64_t> bins{2U, 4U};
    const hbt::StatisticalRegion region{0U, 1U, 6U};
    const std::vector<double> probabilities{1.0 / 3.0, 2.0 / 3.0};
    const double neyman = hbt::binned_neyman_chi_square(
        bins, 0U, region, probabilities, 0.5
    );
    const double pearson = hbt::binned_pearson_chi_square(
        bins, 0U, region, probabilities, 0.5
    );
    if (!close(neyman, 1.5) || !close(pearson, 3.0)) {
        return fail("free normalization did not rescale expected counts");
    }
    try {
        (void)hbt::binned_poisson_deviance(
            bins, 0U, region, probabilities, 0.0
        );
    } catch (const std::invalid_argument&) {
        return true;
    }
    return fail("non-positive free normalization was accepted");
}

/**
 * @brief Verify log-expected objectives reproduce ordinary finite objectives.
 * @return true when only the numerical representation changes.
 */
bool verify_log_expected_objective_consistency() {
    const std::vector<std::uint64_t> bins{2U, 4U};
    const hbt::StatisticalRegion region{0U, 1U, 6U};
    const std::vector<double> probabilities{1.0 / 3.0, 2.0 / 3.0};
    constexpr double normalization = 0.5;
    const std::vector<double> log_expected{
        std::log(6.0 * normalization * probabilities[0U]),
        std::log(6.0 * normalization * probabilities[1U])
    };

    const double poisson = hbt::binned_poisson_deviance(
        bins, 0U, region, probabilities, normalization
    );
    const double neyman = hbt::binned_neyman_chi_square(
        bins, 0U, region, probabilities, normalization
    );
    const double pearson = hbt::binned_pearson_chi_square(
        bins, 0U, region, probabilities, normalization
    );
    if (!close(
            hbt::binned_poisson_deviance_from_log_expected(
                bins, 0U, region, log_expected
            ),
            poisson
        ) ||
        !close(
            hbt::binned_neyman_chi_square_from_log_expected(
                bins, 0U, region, log_expected
            ),
            neyman
        ) ||
        !close(
            hbt::binned_pearson_chi_square_from_log_expected(
                bins, 0U, region, log_expected
            ),
            pearson
        )) {
        return fail("log-expected objectives changed finite objective values");
    }
    return true;
}

/**
 * @brief Verify logarithmic Gaussian integrals remain finite in extreme tails.
 * @return true when direct free-amplitude objectives avoid double underflow.
 */
bool verify_log_gaussian_tail_objective() {
    const hbt::HistogramBinningConfig binning{3U, 0.0, 300.0, 0.01};
    const hbt::StatisticalRegion region{0U, 2U, 111U};
    const std::vector<std::uint64_t> bins{100U, 10U, 1U};

    for (const hbt::FitObservableFamily family : {
             hbt::FitObservableFamily::OSL,
             hbt::FitObservableFamily::Radial
         }) {
        const double log_tail = hbt::gaussian_component_log_integral(
            family, 200.0, 300.0, 2.0
        );
        if (!std::isfinite(log_tail) || !(log_tail < -1000.0)) {
            return fail("far-tail Gaussian log integral is not stable");
        }
        const std::vector<double> log_expected =
            hbt::gaussian_bin_log_expected_counts(
                family, binning, region, 2.0, 1.0
            );
        if (log_expected.size() != 3U ||
            !std::isfinite(log_expected[2U])) {
            return fail("far-tail Gaussian expected counts are not logarithmic");
        }
        const double poisson =
            hbt::binned_poisson_deviance_from_log_expected(
                bins, 0U, region, log_expected
            );
        const double neyman =
            hbt::binned_neyman_chi_square_from_log_expected(
                bins, 0U, region, log_expected
            );
        if (!std::isfinite(poisson) ||
            poisson == std::numeric_limits<double>::max() ||
            !std::isfinite(neyman) ||
            neyman == std::numeric_limits<double>::max()) {
            return fail("free-amplitude Gaussian objective failed in far tail");
        }
    }
    return true;
}

/**
 * @brief Verify logarithmic and ordinary Gaussian integrals agree in range.
 * @return true when the log representation preserves the exact-bin model.
 */
bool verify_log_gaussian_integral_consistency() {
    for (const hbt::FitObservableFamily family : {
             hbt::FitObservableFamily::OSL,
             hbt::FitObservableFamily::Radial
         }) {
        const double integral = hbt::gaussian_component_integral(
            family, 1.2, 1.5, 2.3
        );
        const double log_integral = hbt::gaussian_component_log_integral(
            family, 1.2, 1.5, 2.3
        );
        if (!close(std::exp(log_integral), integral)) {
            return fail("log Gaussian integral changed the exact-bin model");
        }
    }
    return true;
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
        !verify_exact_bin_model_semantics() ||
        !verify_mixed_endpoint_degeneracies() ||
        !verify_large_radius_radial_stability() ||
        !verify_zero_count_likelihood_term() ||
        !verify_neyman_chi_square_terms() ||
        !verify_profiled_mixed_neyman_amplitude() ||
        !verify_mixed_fitted_density() ||
        !verify_pearson_chi_square_terms() ||
        !verify_raw_count_anchor() ||
        !verify_likelihood_probability_normalization() ||
        !verify_free_normalization_multiplier() ||
        !verify_log_expected_objective_consistency() ||
        !verify_log_gaussian_tail_objective() ||
        !verify_log_gaussian_integral_consistency() ||
        !verify_integrated_curve_density()) {
        return EXIT_FAILURE;
    }
    return EXIT_SUCCESS;
}
