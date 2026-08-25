/**
 * @file binned_models.cpp
 * @brief Exact-bin post-sample OSL and radial model calculations.
 */

#include "hbt/fits/binned_models.h"

#include "hbt/fits/statistical_analysis.h"

#include <cmath>
#include <limits>
#include <stdexcept>

namespace hbt {
namespace {

/** Maximum scaled edge for the cancellation-free power series branch. */
constexpr double kSeriesScaledEdge = 0.25;

/** Lower scaled edge above which complementary Gaussian tails are stable. */
constexpr double kGaussianTailEdge = 4.0;

/** Scaled edge above which asymptotic Gaussian survival logs are used. */
constexpr long double kGaussianLogTailEdge = 8.0L;

/**
 * @brief Scaled OSL Gaussian survival with exp(+u^2) factored out.
 * @param value Dimensionless u>0.
 * @return exp(u^2) * integral_u^inf exp(-t^2) dt.
 */
long double scaled_osl_gaussian_survival(long double value) {
    if (!(value > 0.0L) || !std::isfinite(value)) {
        throw std::invalid_argument(
            "HBT analysis model: invalid scaled Gaussian tail edge"
        );
    }
    const long double inverse_two_u2 = 0.5L / (value * value);
    long double term = 1.0L;
    long double sum = 1.0L;
    long double previous_abs = std::fabs(term);
    for (std::size_t order = 1U; order < 256U; ++order) {
        term *= -static_cast<long double>(2U * order - 1U) *
            inverse_two_u2;
        const long double current_abs = std::fabs(term);
        if (current_abs > previous_abs) {
            break;
        }
        const long double next = sum + term;
        if (next == sum) {
            break;
        }
        sum = next;
        previous_abs = current_abs;
    }
    const long double scaled = 0.5L * sum / value;
    if (!std::isfinite(scaled) || !(scaled > 0.0L)) {
        throw std::invalid_argument(
            "HBT analysis model: invalid scaled Gaussian survival"
        );
    }
    return scaled;
}

/** @brief Natural log of integral_u^inf exp(-t^2) dt. */
long double log_osl_gaussian_survival(long double value) {
    if (!std::isfinite(value) || value < 0.0L) {
        throw std::invalid_argument(
            "HBT analysis model: invalid Gaussian survival edge"
        );
    }
    if (value < kGaussianLogTailEdge) {
        const long double pi = std::acos(-1.0L);
        const long double survival = 0.5L * std::sqrt(pi) * std::erfc(value);
        if (!std::isfinite(survival) || !(survival > 0.0L)) {
            throw std::invalid_argument(
                "HBT analysis model: invalid Gaussian survival"
            );
        }
        return std::log(survival);
    }
    return -value * value +
        std::log(scaled_osl_gaussian_survival(value));
}

/** @brief Natural log of integral_u^inf t^2 exp(-t^2) dt. */
long double log_radial_gaussian_survival(long double value) {
    if (!std::isfinite(value) || value < 0.0L) {
        throw std::invalid_argument(
            "HBT analysis model: invalid radial Gaussian survival edge"
        );
    }
    if (value < kGaussianLogTailEdge) {
        const long double pi = std::acos(-1.0L);
        const long double survival =
            0.25L * std::sqrt(pi) * std::erfc(value) +
            0.5L * value * std::exp(-value * value);
        if (!std::isfinite(survival) || !(survival > 0.0L)) {
            throw std::invalid_argument(
                "HBT analysis model: invalid radial Gaussian survival"
            );
        }
        return std::log(survival);
    }
    const long double scaled =
        0.5L * value +
        0.5L * scaled_osl_gaussian_survival(value);
    if (!std::isfinite(scaled) || !(scaled > 0.0L)) {
        throw std::invalid_argument(
            "HBT analysis model: invalid scaled radial Gaussian survival"
        );
    }
    return -value * value + std::log(scaled);
}

/** @brief Stable log(exp(log_a)-exp(log_b)) for log_a>log_b. */
long double log_positive_difference(long double log_a, long double log_b) {
    if (!std::isfinite(log_a) || !std::isfinite(log_b) || !(log_a > log_b)) {
        throw std::invalid_argument(
            "HBT analysis model: invalid logarithmic positive difference"
        );
    }
    const long double delta = log_b - log_a;
    const long double correction = -std::expm1(delta);
    if (!std::isfinite(correction) || !(correction > 0.0L)) {
        throw std::invalid_argument(
            "HBT analysis model: unresolved logarithmic bin integral"
        );
    }
    return log_a + std::log(correction);
}

/**
 * @brief Require valid physical integration inputs.
 * @param lower Non-negative lower edge.
 * @param upper Upper edge greater than lower.
 * @param radius Strictly positive radius.
 * @throws std::invalid_argument If an input is non-finite or out of domain.
 */
void require_integral_inputs(double lower, double upper, double radius) {
    if (!std::isfinite(lower) || !std::isfinite(upper) ||
        !std::isfinite(radius) || lower < 0.0 || upper <= lower ||
        radius <= 0.0) {
        throw std::invalid_argument(
            "HBT analysis model: invalid integration input"
        );
    }
}

/**
 * @brief Validate selected region against configured binning.
 * @param binning Validated uniform histogram binning.
 * @param region Selected contiguous statistical region.
 * @throws std::out_of_range If region indices exceed binning.
 * @throws std::invalid_argument If region order or selected count is invalid.
 */
void require_region(
    const HistogramBinningConfig& binning,
    const StatisticalRegion& region
) {
    if (region.first_bin > region.last_bin ||
        region.last_bin >= binning.nbins) {
        throw std::out_of_range(
            "HBT analysis model: selected region is out of range"
        );
    }
    if (region.selected_count == 0U) {
        throw std::invalid_argument(
            "HBT analysis model: selected count is zero"
        );
    }
}

/**
 * @brief Return x^power for a finite non-negative physical edge.
 * @param value Non-negative finite edge.
 * @param power Non-negative integer power.
 * @return Finite power when representable.
 * @throws std::invalid_argument If the power is not finite.
 */
double finite_power(double value, std::size_t power) {
    const double result = std::pow(value, static_cast<double>(power));
    if (!std::isfinite(result)) {
        throw std::invalid_argument(
            "HBT analysis model: integration power overflow"
        );
    }
    return result;
}

/**
 * @brief Integrate x^power times exp(-x^2/(4R^2)) by its stable series.
 * @param lower Non-negative lower edge.
 * @param upper Strictly larger upper edge.
 * @param radius Strictly positive radius.
 * @param power Non-negative integer x power.
 * @return Positive finite exact-edge integral at machine precision.
 * @throws std::invalid_argument If the series becomes non-finite or invalid.
 *
 * This branch is used only when upper/(2R) <= 0.25, where the alternating
 * Taylor series decreases rapidly. Iteration stops only when the next term no
 * longer changes the representable floating-point sum; no scientific epsilon
 * or physical cutoff is introduced.
 */
double gaussian_large_radius_series(
    double lower,
    double upper,
    double radius,
    std::size_t power
) {
    const double radius2 = radius * radius;
    if (!std::isfinite(radius2)) {
        const double leading = (
            finite_power(upper, power + 1U) -
            finite_power(lower, power + 1U)
        ) / static_cast<double>(power + 1U);
        if (!std::isfinite(leading) || leading <= 0.0) {
            throw std::invalid_argument(
                "HBT analysis model: invalid Gaussian series integral"
            );
        }
        return leading;
    }

    double factorial = 1.0;
    double scale_power = 1.0;
    double sum = 0.0;
    for (std::size_t term_index = 0U;; ++term_index) {
        if (term_index > 0U) {
            factorial *= static_cast<double>(term_index);
            scale_power *= 4.0 * radius2;
        }
        const std::size_t exponent = power + 2U * term_index + 1U;
        const double edge_difference =
            finite_power(upper, exponent) - finite_power(lower, exponent);
        const double denominator = factorial * scale_power *
            static_cast<double>(exponent);
        if (!std::isfinite(denominator) || denominator <= 0.0) {
            break;
        }
        double term = edge_difference / denominator;
        if ((term_index % 2U) != 0U) {
            term = -term;
        }
        if (!std::isfinite(term)) {
            throw std::invalid_argument(
                "HBT analysis model: non-finite Gaussian series term"
            );
        }
        const double next = sum + term;
        if (next == sum) {
            break;
        }
        sum = next;
    }

    if (!std::isfinite(sum) || sum <= 0.0) {
        throw std::invalid_argument(
            "HBT analysis model: invalid Gaussian series integral"
        );
    }
    return sum;
}

/**
 * @brief Integrate x^power times exp(-x/R) by its stable large-R series.
 * @param lower Non-negative lower edge.
 * @param upper Strictly larger upper edge.
 * @param radius Strictly positive radius.
 * @param power Non-negative integer x power.
 * @return Positive finite exact-edge integral at machine precision.
 * @throws std::invalid_argument If the series becomes non-finite or invalid.
 *
 * This branch is used only when upper/R <= 0.25. Terms are summed until the
 * next term is below representable resolution of the accumulated result.
 */
double exponential_large_radius_series(
    double lower,
    double upper,
    double radius,
    std::size_t power
) {
    double factorial = 1.0;
    double radius_power = 1.0;
    double sum = 0.0;
    for (std::size_t term_index = 0U;; ++term_index) {
        if (term_index > 0U) {
            factorial *= static_cast<double>(term_index);
            radius_power *= radius;
        }
        const std::size_t exponent = power + term_index + 1U;
        const double edge_difference =
            finite_power(upper, exponent) - finite_power(lower, exponent);
        const double denominator = factorial * radius_power *
            static_cast<double>(exponent);
        if (!std::isfinite(denominator) || denominator <= 0.0) {
            break;
        }
        double term = edge_difference / denominator;
        if ((term_index % 2U) != 0U) {
            term = -term;
        }
        if (!std::isfinite(term)) {
            throw std::invalid_argument(
                "HBT analysis model: non-finite exponential series term"
            );
        }
        const double next = sum + term;
        if (next == sum) {
            break;
        }
        sum = next;
    }

    if (!std::isfinite(sum) || sum <= 0.0) {
        throw std::invalid_argument(
            "HBT analysis model: invalid exponential series integral"
        );
    }
    return sum;
}

/**
 * @brief Return the dimensionless OSL Gaussian exact-bin weight.
 * @param lower Non-negative lower edge.
 * @param upper Strictly larger upper edge.
 * @param radius Strictly positive radius.
 * @return Non-negative Gaussian weight with the common scale removed.
 */
double osl_gaussian_weight(double lower, double upper, double radius) {
    const double lower_scaled = lower / (2.0 * radius);
    const double upper_scaled = upper / (2.0 * radius);
    if (lower_scaled >= kGaussianTailEdge) {
        return std::erfc(lower_scaled) - std::erfc(upper_scaled);
    }
    return std::erf(upper_scaled) - std::erf(lower_scaled);
}

/**
 * @brief Return the radial Gaussian scaled antiderivative F(u).
 * @param value Dimensionless non-negative u = r/(2R).
 * @return F(u) with F'(u) = u^2 exp(-u^2).
 */
double radial_gaussian_cdf(double value) {
    const double pi = std::acos(-1.0);
    return 0.25 * std::sqrt(pi) * std::erf(value) -
        0.5 * value * std::exp(-value * value);
}

/**
 * @brief Return the radial Gaussian scaled survival integral T(u).
 * @param value Dimensionless non-negative u = r/(2R).
 * @return Integral from u to infinity of t^2 exp(-t^2) dt.
 */
double radial_gaussian_survival(double value) {
    const double pi = std::acos(-1.0);
    return 0.25 * std::sqrt(pi) * std::erfc(value) +
        0.5 * value * std::exp(-value * value);
}

/**
 * @brief Return a dimensionless radial Gaussian exact-bin weight.
 * @param lower Non-negative lower edge.
 * @param upper Strictly larger upper edge.
 * @param radius Strictly positive radius.
 * @return Non-negative weight with the common 8R^3 scale removed.
 */
double radial_gaussian_weight(
    double lower,
    double upper,
    double radius
) {
    const double lower_scaled = lower / (2.0 * radius);
    const double upper_scaled = upper / (2.0 * radius);
    if (upper_scaled <= kSeriesScaledEdge) {
        const double integral = gaussian_large_radius_series(
            lower,
            upper,
            radius,
            2U
        );
        const double radius3 = radius * radius * radius;
        if (!std::isfinite(radius3) || radius3 <= 0.0) {
            return integral;
        }
        return integral / (8.0 * radius3);
    }
    if (lower_scaled >= kGaussianTailEdge) {
        return radial_gaussian_survival(lower_scaled) -
            radial_gaussian_survival(upper_scaled);
    }
    return radial_gaussian_cdf(upper_scaled) -
        radial_gaussian_cdf(lower_scaled);
}

/**
 * @brief Return a dimensionless OSL exponential exact-bin weight.
 * @param lower Non-negative lower edge.
 * @param upper Strictly larger upper edge.
 * @param radius Strictly positive radius.
 * @return Non-negative weight with the common R scale removed.
 */
double osl_exponential_weight(double lower, double upper, double radius) {
    const double lower_scaled = lower / radius;
    const double width_scaled = (upper - lower) / radius;
    return std::exp(-lower_scaled) * (-std::expm1(-width_scaled));
}

/**
 * @brief Return the radial exponential scaled survival integral.
 * @param value Dimensionless non-negative u = r/R.
 * @return Integral from u to infinity of t^2 exp(-t) dt.
 */
double radial_exponential_survival(double value) {
    return std::exp(-value) *
        (value * value + 2.0 * value + 2.0);
}

/**
 * @brief Return a dimensionless radial exponential exact-bin weight.
 * @param lower Non-negative lower edge.
 * @param upper Strictly larger upper edge.
 * @param radius Strictly positive radius.
 * @return Non-negative weight with the common R^3 scale removed.
 */
double radial_exponential_weight(
    double lower,
    double upper,
    double radius
) {
    const double upper_scaled = upper / radius;
    if (upper_scaled <= kSeriesScaledEdge) {
        const double integral = exponential_large_radius_series(
            lower,
            upper,
            radius,
            2U
        );
        const double radius3 = radius * radius * radius;
        if (!std::isfinite(radius3) || radius3 <= 0.0) {
            return integral;
        }
        return integral / radius3;
    }
    const double lower_scaled = lower / radius;
    return radial_exponential_survival(lower_scaled) -
        radial_exponential_survival(upper_scaled);
}

/**
 * @brief Calculate one unnormalized exact-bin shape weight.
 * @param family OSL or radial physical family.
 * @param lower Non-negative lower edge.
 * @param upper Strictly larger upper edge.
 * @param radius Strictly positive radius.
 * @param gaussian true for Gaussian, false for exponential.
 * @return Non-negative finite weight with any common R scale removed.
 */
double component_weight(
    FitObservableFamily family,
    double lower,
    double upper,
    double radius,
    bool gaussian
) {
    require_integral_inputs(lower, upper, radius);
    double weight = 0.0;
    if (gaussian) {
        weight = family == FitObservableFamily::OSL
            ? osl_gaussian_weight(lower, upper, radius)
            : radial_gaussian_weight(lower, upper, radius);
    } else {
        weight = family == FitObservableFamily::OSL
            ? osl_exponential_weight(lower, upper, radius)
            : radial_exponential_weight(lower, upper, radius);
    }
    if (!std::isfinite(weight) || weight < 0.0) {
        throw std::invalid_argument(
            "HBT analysis model: invalid component bin weight"
        );
    }
    return weight;
}

/**
 * @brief Normalize non-negative exact-bin component weights.
 * @param weights Finite non-negative per-bin component weights.
 * @param require_positive Whether every normalized bin must be positive.
 * @return Unit-sum probabilities.
 * @throws std::invalid_argument If normalization or a required bin is invalid.
 */
std::vector<double> normalize_weights(
    const std::vector<double>& weights,
    bool require_positive
) {
    double normalization = 0.0;
    for (const double weight : weights) {
        if (!std::isfinite(weight) || weight < 0.0) {
            throw std::invalid_argument(
                "HBT analysis model: invalid component bin weight"
            );
        }
        normalization += weight;
    }
    if (!std::isfinite(normalization) || normalization <= 0.0) {
        throw std::invalid_argument(
            "HBT analysis model: invalid component normalization"
        );
    }

    std::vector<double> probabilities;
    probabilities.reserve(weights.size());
    for (const double weight : weights) {
        const double probability = weight / normalization;
        if (!std::isfinite(probability) || probability < 0.0 ||
            (require_positive && probability <= 0.0)) {
            throw std::invalid_argument(
                "HBT analysis model: invalid normalized bin probability"
            );
        }
        probabilities.push_back(probability);
    }
    return probabilities;
}

/**
 * @brief Calculate one component's exact selected-bin probabilities.
 * @param family OSL or radial physical family.
 * @param binning Validated histogram binning.
 * @param region Selected contiguous region.
 * @param radius Strictly positive component radius.
 * @param gaussian true for Gaussian, false for exponential.
 * @param require_positive Whether all returned probabilities must be positive.
 * @return Unit-sum exact-bin probabilities.
 */
std::vector<double> component_probabilities(
    FitObservableFamily family,
    const HistogramBinningConfig& binning,
    const StatisticalRegion& region,
    double radius,
    bool gaussian,
    bool require_positive
) {
    require_region(binning, region);
    std::vector<double> weights;
    weights.reserve(region.last_bin - region.first_bin + 1U);

    for (std::size_t bin = region.first_bin;
         bin <= region.last_bin;
         ++bin) {
        const double lower = histogram_bin_lower_edge(binning, bin);
        const double upper = histogram_bin_upper_edge(binning, bin);
        weights.push_back(
            component_weight(family, lower, upper, radius, gaussian)
        );
    }
    return normalize_weights(weights, require_positive);
}

}  // namespace

double gaussian_component_log_integral(
    FitObservableFamily family,
    double lower,
    double upper,
    double radius
) {
    require_integral_inputs(lower, upper, radius);
    const long double radius_ld = static_cast<long double>(radius);
    const long double lower_scaled =
        static_cast<long double>(lower) / (2.0L * radius_ld);
    const long double upper_scaled =
        static_cast<long double>(upper) / (2.0L * radius_ld);

    long double log_integral = 0.0L;
    if (upper_scaled <= static_cast<long double>(kSeriesScaledEdge)) {
        const double integral = gaussian_large_radius_series(
            lower,
            upper,
            radius,
            family == FitObservableFamily::OSL ? 0U : 2U
        );
        log_integral = std::log(static_cast<long double>(integral));
    } else {
        const long double log_lower_survival =
            family == FitObservableFamily::OSL
                ? log_osl_gaussian_survival(lower_scaled)
                : log_radial_gaussian_survival(lower_scaled);
        const long double log_upper_survival =
            family == FitObservableFamily::OSL
                ? log_osl_gaussian_survival(upper_scaled)
                : log_radial_gaussian_survival(upper_scaled);
        const long double log_dimensionless = log_positive_difference(
            log_lower_survival, log_upper_survival
        );
        if (family == FitObservableFamily::OSL) {
            log_integral = std::log(2.0L * radius_ld) + log_dimensionless;
        } else {
            log_integral = std::log(8.0L * radius_ld * radius_ld * radius_ld) +
                log_dimensionless;
        }
    }

    const double result = static_cast<double>(log_integral);
    if (!std::isfinite(result)) {
        throw std::invalid_argument(
            "HBT analysis model: invalid logarithmic Gaussian integral"
        );
    }
    return result;
}

std::vector<double> gaussian_bin_log_expected_counts(
    FitObservableFamily family,
    const HistogramBinningConfig& binning,
    const StatisticalRegion& region,
    double radius,
    double amplitude
) {
    require_region(binning, region);
    if (!std::isfinite(radius) || radius <= 0.0 ||
        !std::isfinite(amplitude) || amplitude <= 0.0) {
        throw std::invalid_argument(
            "HBT Gaussian fit: invalid radius or amplitude"
        );
    }
    const double log_scale =
        std::log(static_cast<double>(region.selected_count)) +
        std::log(amplitude);
    if (!std::isfinite(log_scale)) {
        throw std::invalid_argument(
            "HBT Gaussian fit: invalid expected-count scale"
        );
    }

    std::vector<double> log_expected;
    log_expected.reserve(region.last_bin - region.first_bin + 1U);
    for (std::size_t bin = region.first_bin; bin <= region.last_bin; ++bin) {
        const double log_integral = gaussian_component_log_integral(
            family,
            histogram_bin_lower_edge(binning, bin),
            histogram_bin_upper_edge(binning, bin),
            radius
        );
        const double value = log_scale + log_integral;
        if (!std::isfinite(value)) {
            throw std::invalid_argument(
                "HBT Gaussian fit: invalid logarithmic expected count"
            );
        }
        log_expected.push_back(value);
    }
    return log_expected;
}

double gaussian_component_integral(
    FitObservableFamily family,
    double lower,
    double upper,
    double radius
) {
    require_integral_inputs(lower, upper, radius);

    double result = 0.0;
    switch (family) {
        case FitObservableFamily::OSL:
            if (upper / (2.0 * radius) <= kSeriesScaledEdge) {
                result = gaussian_large_radius_series(
                    lower,
                    upper,
                    radius,
                    0U
                );
            } else {
                const double pi = std::acos(-1.0);
                result = radius * std::sqrt(pi) *
                    osl_gaussian_weight(lower, upper, radius);
            }
            break;
        case FitObservableFamily::Radial:
            if (upper / (2.0 * radius) <= kSeriesScaledEdge) {
                result = gaussian_large_radius_series(
                    lower,
                    upper,
                    radius,
                    2U
                );
            } else {
                const double radius3 = radius * radius * radius;
                result = 8.0 * radius3 *
                    radial_gaussian_weight(lower, upper, radius);
            }
            break;
    }

    if (!std::isfinite(result) || result <= 0.0) {
        throw std::invalid_argument(
            "HBT analysis model: invalid Gaussian component integral"
        );
    }
    return result;
}

double exponential_component_integral(
    FitObservableFamily family,
    double lower,
    double upper,
    double radius
) {
    require_integral_inputs(lower, upper, radius);

    double result = 0.0;
    switch (family) {
        case FitObservableFamily::OSL:
            if (upper / radius <= kSeriesScaledEdge) {
                result = exponential_large_radius_series(
                    lower,
                    upper,
                    radius,
                    0U
                );
            } else {
                result = radius *
                    osl_exponential_weight(lower, upper, radius);
            }
            break;
        case FitObservableFamily::Radial:
            if (upper / radius <= kSeriesScaledEdge) {
                result = exponential_large_radius_series(
                    lower,
                    upper,
                    radius,
                    2U
                );
            } else {
                const double radius3 = radius * radius * radius;
                result = radius3 *
                    radial_exponential_weight(lower, upper, radius);
            }
            break;
    }

    if (!std::isfinite(result) || result <= 0.0) {
        throw std::invalid_argument(
            "HBT analysis model: invalid exponential component integral"
        );
    }
    return result;
}

std::vector<double> gaussian_bin_probabilities(
    FitObservableFamily family,
    const HistogramBinningConfig& binning,
    const StatisticalRegion& region,
    double radius
) {
    return component_probabilities(
        family,
        binning,
        region,
        radius,
        true,
        true
    );
}

std::vector<double> mixed_bin_probabilities(
    FitObservableFamily family,
    const HistogramBinningConfig& binning,
    const StatisticalRegion& region,
    double core_radius,
    double tail_radius,
    double core_fraction
) {
    if (!std::isfinite(core_fraction) || core_fraction < 0.0 ||
        core_fraction > 1.0) {
        throw std::invalid_argument(
            "HBT analysis model: core fraction is outside [0,1]"
        );
    }

    if (core_fraction == 1.0) {
        return component_probabilities(
            family,
            binning,
            region,
            core_radius,
            true,
            true
        );
    }
    if (core_fraction == 0.0) {
        return component_probabilities(
            family,
            binning,
            region,
            tail_radius,
            false,
            true
        );
    }

    const std::vector<double> gaussian = component_probabilities(
        family,
        binning,
        region,
        core_radius,
        true,
        false
    );
    const std::vector<double> exponential = component_probabilities(
        family,
        binning,
        region,
        tail_radius,
        false,
        false
    );

    std::vector<double> probabilities;
    probabilities.reserve(gaussian.size());
    for (std::size_t index = 0U; index < gaussian.size(); ++index) {
        const double probability = core_fraction * gaussian[index] +
            (1.0 - core_fraction) * exponential[index];
        if (!std::isfinite(probability) || probability <= 0.0) {
            throw std::invalid_argument(
                "HBT analysis model: invalid mixed-model bin probability"
            );
        }
        probabilities.push_back(probability);
    }
    return probabilities;
}

double binned_poisson_deviance(
    const std::vector<std::uint64_t>& bins,
    std::size_t offset,
    const StatisticalRegion& region,
    const std::vector<double>& probabilities,
    double normalization
) {
    if (!std::isfinite(normalization) || normalization <= 0.0) {
        throw std::invalid_argument(
            "HBT analysis objective: invalid model normalization"
        );
    }
    const std::size_t selected_bins =
        region.last_bin - region.first_bin + 1U;
    if (probabilities.size() != selected_bins) {
        throw std::invalid_argument(
            "HBT analysis objective: probability cardinality mismatch"
        );
    }
    if (offset > bins.size() ||
        region.first_bin > region.last_bin ||
        region.first_bin > bins.size() - offset ||
        selected_bins > bins.size() - offset - region.first_bin) {
        throw std::out_of_range(
            "HBT analysis objective: selected raw histogram range is "
            "unavailable"
        );
    }
    if (region.selected_count == 0U) {
        throw std::invalid_argument(
            "HBT analysis objective: selected count is zero"
        );
    }

    std::uint64_t raw_sum = 0U;
    long double probability_sum = 0.0L;
    for (std::size_t index = 0U; index < selected_bins; ++index) {
        const std::uint64_t raw = bins[offset + region.first_bin + index];
        if (raw > std::numeric_limits<std::uint64_t>::max() - raw_sum) {
            throw std::overflow_error(
                "HBT analysis objective: selected raw count overflow"
            );
        }
        raw_sum += raw;
        const double probability = probabilities[index];
        if (!std::isfinite(probability) || probability <= 0.0) {
            throw std::invalid_argument(
                "HBT analysis objective: invalid model probability"
            );
        }
        probability_sum += static_cast<long double>(probability);
    }
    if (raw_sum != region.selected_count) {
        throw std::invalid_argument(
            "HBT analysis objective: selected count differs from raw counts"
        );
    }

    const long double normalization_tolerance =
        static_cast<long double>(std::numeric_limits<double>::epsilon()) *
        static_cast<long double>(selected_bins);
    if (!std::isfinite(probability_sum) ||
        std::fabs(probability_sum - 1.0L) > normalization_tolerance) {
        throw std::invalid_argument(
            "HBT analysis objective: model probabilities are not normalized"
        );
    }

    double deviance = 0.0;
    for (std::size_t index = 0U; index < selected_bins; ++index) {
        const double expected = static_cast<double>(region.selected_count) *
            normalization * probabilities[index];
        if (!std::isfinite(expected) || expected <= 0.0) {
            throw std::invalid_argument(
                "HBT analysis objective: invalid expected count"
            );
        }
        const std::uint64_t raw = bins[offset + region.first_bin + index];
        if (raw == 0U) {
            deviance += 2.0 * expected;
        } else {
            const double observed = static_cast<double>(raw);
            deviance += 2.0 * (
                expected - observed +
                observed * std::log(observed / expected)
            );
        }
    }
    if (!std::isfinite(deviance)) {
        throw std::invalid_argument(
            "HBT analysis objective: non-finite Poisson deviance"
        );
    }
    return deviance;
}

namespace {

/**
 * @brief Validate common binned-estimator bookkeeping and normalization.
 * @param bins Raw histogram counts.
 * @param offset First raw counter belonging to the logical histogram.
 * @param region Selected contiguous statistical region.
 * @param probabilities Model probabilities aligned with selected bins.
 * @return Number of selected bins after all checks pass.
 * @throws std::invalid_argument For invalid probabilities or count bookkeeping.
 * @throws std::out_of_range If the selected raw range is unavailable.
 */
std::size_t validate_chi_square_inputs(
    const std::vector<std::uint64_t>& bins,
    std::size_t offset,
    const StatisticalRegion& region,
    const std::vector<double>& probabilities
) {
    const std::size_t selected_bins =
        region.last_bin - region.first_bin + 1U;
    if (probabilities.size() != selected_bins) {
        throw std::invalid_argument(
            "HBT analysis chi-square: probability cardinality mismatch"
        );
    }
    if (offset > bins.size() ||
        region.first_bin > region.last_bin ||
        region.first_bin > bins.size() - offset ||
        selected_bins > bins.size() - offset - region.first_bin) {
        throw std::out_of_range(
            "HBT analysis chi-square: selected raw histogram range is unavailable"
        );
    }
    if (region.selected_count == 0U) {
        throw std::invalid_argument(
            "HBT analysis chi-square: selected count is zero"
        );
    }

    std::uint64_t raw_sum = 0U;
    long double probability_sum = 0.0L;
    for (std::size_t index = 0U; index < selected_bins; ++index) {
        const std::uint64_t raw = bins[offset + region.first_bin + index];
        if (raw > std::numeric_limits<std::uint64_t>::max() - raw_sum) {
            throw std::overflow_error(
                "HBT analysis chi-square: selected raw count overflow"
            );
        }
        raw_sum += raw;
        const double probability = probabilities[index];
        if (!std::isfinite(probability) || probability <= 0.0) {
            throw std::invalid_argument(
                "HBT analysis chi-square: invalid model probability"
            );
        }
        probability_sum += static_cast<long double>(probability);
    }
    if (raw_sum != region.selected_count) {
        throw std::invalid_argument(
            "HBT analysis chi-square: selected count differs from raw counts"
        );
    }
    const long double normalization_tolerance =
        static_cast<long double>(std::numeric_limits<double>::epsilon()) *
        static_cast<long double>(selected_bins);
    if (!std::isfinite(probability_sum) ||
        std::fabs(probability_sum - 1.0L) > normalization_tolerance) {
        throw std::invalid_argument(
            "HBT analysis chi-square: model probabilities are not normalized"
        );
    }
    return selected_bins;
}

}  // namespace

double binned_neyman_chi_square(
    const std::vector<std::uint64_t>& bins,
    std::size_t offset,
    const StatisticalRegion& region,
    const std::vector<double>& probabilities,
    double normalization
) {
    if (!std::isfinite(normalization) || normalization <= 0.0) {
        throw std::invalid_argument(
            "HBT analysis Neyman chi-square: invalid model normalization"
        );
    }
    const std::size_t selected_bins = validate_chi_square_inputs(
        bins, offset, region, probabilities
    );
    double chi_square = 0.0;
    for (std::size_t index = 0U; index < selected_bins; ++index) {
        const std::uint64_t raw = bins[offset + region.first_bin + index];
        if (raw == 0U) {
            continue;
        }
        const double observed = static_cast<double>(raw);
        const double expected = static_cast<double>(region.selected_count) *
            normalization * probabilities[index];
        if (!std::isfinite(expected) || expected <= 0.0) {
            throw std::invalid_argument(
                "HBT analysis Neyman chi-square: invalid expected count"
            );
        }
        const double residual = observed - expected;
        chi_square += residual * residual / observed;
    }
    if (!std::isfinite(chi_square)) {
        throw std::invalid_argument(
            "HBT analysis Neyman chi-square: non-finite objective evaluation"
        );
    }
    return chi_square;
}

double binned_pearson_chi_square(
    const std::vector<std::uint64_t>& bins,
    std::size_t offset,
    const StatisticalRegion& region,
    const std::vector<double>& probabilities,
    double normalization
) {
    if (!std::isfinite(normalization) || normalization <= 0.0) {
        throw std::invalid_argument(
            "HBT analysis Pearson chi-square: invalid model normalization"
        );
    }
    const std::size_t selected_bins = validate_chi_square_inputs(
        bins, offset, region, probabilities
    );
    double chi_square = 0.0;
    for (std::size_t index = 0U; index < selected_bins; ++index) {
        const double observed = static_cast<double>(
            bins[offset + region.first_bin + index]
        );
        const double expected = static_cast<double>(region.selected_count) *
            normalization * probabilities[index];
        if (!std::isfinite(expected) || expected <= 0.0) {
            throw std::invalid_argument(
                "HBT analysis Pearson chi-square: invalid expected count"
            );
        }
        const double residual = observed - expected;
        chi_square += residual * residual / expected;
    }
    if (!std::isfinite(chi_square)) {
        throw std::invalid_argument(
            "HBT analysis Pearson chi-square: non-finite objective evaluation"
        );
    }
    return chi_square;
}

namespace {

/** @brief Validate raw-count bookkeeping for logarithmic expected counts. */
std::size_t validate_log_expected_inputs(
    const std::vector<std::uint64_t>& bins,
    std::size_t offset,
    const StatisticalRegion& region,
    const std::vector<double>& log_expected
) {
    const std::size_t selected_bins =
        region.last_bin - region.first_bin + 1U;
    if (log_expected.size() != selected_bins) {
        throw std::invalid_argument(
            "HBT analysis objective: log-expected cardinality mismatch"
        );
    }
    if (offset > bins.size() ||
        region.first_bin > region.last_bin ||
        region.first_bin > bins.size() - offset ||
        selected_bins > bins.size() - offset - region.first_bin) {
        throw std::out_of_range(
            "HBT analysis objective: selected raw histogram range is unavailable"
        );
    }
    if (region.selected_count == 0U) {
        throw std::invalid_argument(
            "HBT analysis objective: selected count is zero"
        );
    }

    std::uint64_t raw_sum = 0U;
    for (std::size_t index = 0U; index < selected_bins; ++index) {
        const std::uint64_t raw = bins[offset + region.first_bin + index];
        if (raw > std::numeric_limits<std::uint64_t>::max() - raw_sum) {
            throw std::overflow_error(
                "HBT analysis objective: selected raw count overflow"
            );
        }
        raw_sum += raw;
        if (!std::isfinite(log_expected[index])) {
            throw std::invalid_argument(
                "HBT analysis objective: invalid logarithmic expected count"
            );
        }
    }
    if (raw_sum != region.selected_count) {
        throw std::invalid_argument(
            "HBT analysis objective: selected count differs from raw counts"
        );
    }
    return selected_bins;
}

/** @brief Exponentiate one finite log count using long-double dynamic range. */
long double expected_from_log(double log_expected) {
    const long double log_value = static_cast<long double>(log_expected);
    const long double log_max = std::log(
        std::numeric_limits<long double>::max()
    );
    if (log_value > log_max) {
        return std::numeric_limits<long double>::infinity();
    }
    return std::exp(log_value);
}

/** @brief Convert a non-negative long-double objective to Minuit double. */
double finite_objective_value(long double value) {
    if (!std::isfinite(value) || value < 0.0L ||
        value >= static_cast<long double>(
            std::numeric_limits<double>::max()
        )) {
        return std::numeric_limits<double>::max();
    }
    return static_cast<double>(value);
}

}  // namespace

double binned_poisson_deviance_from_log_expected(
    const std::vector<std::uint64_t>& bins,
    std::size_t offset,
    const StatisticalRegion& region,
    const std::vector<double>& log_expected
) {
    const std::size_t selected_bins = validate_log_expected_inputs(
        bins, offset, region, log_expected
    );
    long double deviance = 0.0L;
    for (std::size_t index = 0U; index < selected_bins; ++index) {
        const std::uint64_t raw = bins[offset + region.first_bin + index];
        const long double expected = expected_from_log(log_expected[index]);
        if (!std::isfinite(expected)) {
            return std::numeric_limits<double>::max();
        }
        if (raw == 0U) {
            deviance += 2.0L * expected;
        } else {
            const long double observed = static_cast<long double>(raw);
            deviance += 2.0L * (
                expected - observed + observed * (
                    std::log(observed) -
                    static_cast<long double>(log_expected[index])
                )
            );
        }
        if (!std::isfinite(deviance)) {
            return std::numeric_limits<double>::max();
        }
    }
    if (deviance < 0.0L &&
        std::fabs(deviance) <=
            64.0L * std::numeric_limits<long double>::epsilon()) {
        deviance = 0.0L;
    }
    return finite_objective_value(deviance);
}

double binned_neyman_chi_square_from_log_expected(
    const std::vector<std::uint64_t>& bins,
    std::size_t offset,
    const StatisticalRegion& region,
    const std::vector<double>& log_expected
) {
    const std::size_t selected_bins = validate_log_expected_inputs(
        bins, offset, region, log_expected
    );
    long double chi_square = 0.0L;
    for (std::size_t index = 0U; index < selected_bins; ++index) {
        const std::uint64_t raw = bins[offset + region.first_bin + index];
        if (raw == 0U) {
            continue;
        }
        const long double expected = expected_from_log(log_expected[index]);
        if (!std::isfinite(expected)) {
            return std::numeric_limits<double>::max();
        }
        const long double observed = static_cast<long double>(raw);
        const long double residual = observed - expected;
        chi_square += residual * residual / observed;
        if (!std::isfinite(chi_square)) {
            return std::numeric_limits<double>::max();
        }
    }
    return finite_objective_value(chi_square);
}

double binned_pearson_chi_square_from_log_expected(
    const std::vector<std::uint64_t>& bins,
    std::size_t offset,
    const StatisticalRegion& region,
    const std::vector<double>& log_expected
) {
    const std::size_t selected_bins = validate_log_expected_inputs(
        bins, offset, region, log_expected
    );
    long double chi_square = 0.0L;
    for (std::size_t index = 0U; index < selected_bins; ++index) {
        const long double expected = expected_from_log(log_expected[index]);
        if (!std::isfinite(expected)) {
            return std::numeric_limits<double>::max();
        }
        const long double observed = static_cast<long double>(
            bins[offset + region.first_bin + index]
        );
        if (expected == 0.0L) {
            if (observed == 0.0L) {
                continue;
            }
            return std::numeric_limits<double>::max();
        }
        const long double residual = observed - expected;
        chi_square += residual * residual / expected;
        if (!std::isfinite(chi_square)) {
            return std::numeric_limits<double>::max();
        }
    }
    return finite_objective_value(chi_square);
}

std::vector<double> probabilities_to_pdf(
    const std::vector<double>& probabilities,
    const HistogramBinningConfig& binning
) {
    const double width = 1.0 / binning.inverse_bin_width;
    std::vector<double> densities;
    densities.reserve(probabilities.size());
    for (const double probability : probabilities) {
        if (!std::isfinite(probability) || probability < 0.0) {
            throw std::invalid_argument(
                "HBT analysis model: invalid probability for presentation PDF"
            );
        }
        const double density = probability / width;
        if (!std::isfinite(density)) {
            throw std::invalid_argument(
                "HBT analysis model: non-finite presentation PDF"
            );
        }
        densities.push_back(density);
    }
    return densities;
}

}  // namespace hbt
