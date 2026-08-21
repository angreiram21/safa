/**
 * @file binned_models.h
 * @brief Exact-bin post-sample OSL and radial model calculations.
 */

#ifndef HBT_FITS_BINNED_MODELS_H
#define HBT_FITS_BINNED_MODELS_H

#include "hbt/config/hbt_config.h"
#include "hbt/fits/fit_results.h"

#include <cstddef>
#include <cstdint>
#include <vector>

namespace hbt {

/**
 * @brief Integrate the unnormalized Gaussian component over exact edges.
 * @param family OSL or radial physical model family.
 * @param lower Non-negative inclusive lower edge.
 * @param upper Finite upper edge strictly greater than lower.
 * @param radius Strictly positive finite Gaussian radius R.
 * @return Exact analytic component integral over [lower, upper].
 * @throws std::invalid_argument If any precondition is violated or the
 *         resulting integral is not finite and strictly positive.
 */
[[nodiscard]] double gaussian_component_integral(
    FitObservableFamily family,
    double lower,
    double upper,
    double radius
);

/**
 * @brief Integrate the unnormalized exponential component over exact edges.
 * @param family OSL or radial physical model family.
 * @param lower Non-negative inclusive lower edge.
 * @param upper Finite upper edge strictly greater than lower.
 * @param radius Strictly positive finite exponential tail radius.
 * @return Exact analytic component integral over [lower, upper].
 * @throws std::invalid_argument If any precondition is violated or the
 *         resulting integral is not finite and strictly positive.
 */
[[nodiscard]] double exponential_component_integral(
    FitObservableFamily family,
    double lower,
    double upper,
    double radius
);

/**
 * @brief Return exact-bin probabilities for a normalized Gaussian model.
 * @param family OSL or radial physical model family.
 * @param binning Validated uniform histogram binning starting at zero.
 * @param region Selected contiguous statistical region.
 * @param radius Strictly positive finite Gaussian radius R.
 * @return One normalized probability per selected bin.
 * @throws std::invalid_argument If model evaluation is invalid.
 * @throws std::out_of_range If the region exceeds configured binning.
 *
 * The component is normalized to unit integral over exactly the selected
 * region. No free amplitude or bin-center approximation is used.
 */
[[nodiscard]] std::vector<double> gaussian_bin_probabilities(
    FitObservableFamily family,
    const HistogramBinningConfig& binning,
    const StatisticalRegion& region,
    double radius
);

/**
 * @brief Return exact-bin probabilities for the normalized mixed model.
 * @param family OSL or radial physical model family.
 * @param binning Validated uniform histogram binning starting at zero.
 * @param region Selected contiguous statistical region.
 * @param core_radius Strictly positive finite Gaussian core radius.
 * @param tail_radius Strictly positive finite exponential tail radius.
 * @param core_fraction Physical mixture fraction in [0,1].
 * @return One normalized mixture probability per selected bin.
 * @throws std::invalid_argument If model evaluation is invalid.
 * @throws std::out_of_range If the region exceeds configured binning.
 *
 * Gaussian and exponential components are normalized independently over the
 * same selected region before mixing.
 */
[[nodiscard]] std::vector<double> mixed_bin_probabilities(
    FitObservableFamily family,
    const HistogramBinningConfig& binning,
    const StatisticalRegion& region,
    double core_radius,
    double tail_radius,
    double core_fraction
);

/**
 * @brief Evaluate binned Poisson deviance from counts and probabilities.
 * @param bins Slot-major raw histogram count storage.
 * @param offset First raw counter belonging to the logical histogram.
 * @param region Selected contiguous statistical region.
 * @param probabilities Model probabilities aligned with selected bins.
 * @return D = -2 log L up to shape-independent constants.
 * @throws std::invalid_argument If probabilities are non-positive, non-finite,
 *         not unit-normalized within double roundoff, cardinality differs, or
 *         selected_count differs from raw counts.
 * @throws std::out_of_range If the selected raw range is unavailable.
 *
 * Observed n_i == 0 bins remain in the sum and contribute exactly 2*mu_i.
 * The unit-sum allowance is derived from machine epsilon and selected-bin
 * count only. No scientific epsilon, renormalization, Gaussian count errors,
 * or empty-bin removal is used.
 */
[[nodiscard]] double binned_poisson_deviance(
    const std::vector<std::uint64_t>& bins,
    std::size_t offset,
    const StatisticalRegion& region,
    const std::vector<double>& probabilities
);

/**
 * @brief Evaluate Neyman chi-square from raw counts and model probabilities.
 * @param bins Slot-major raw histogram count storage.
 * @param offset First raw counter belonging to the logical histogram.
 * @param region Selected contiguous statistical region.
 * @param probabilities Unit-normalized model probabilities for selected bins.
 * @return chi2_N = sum_{n_i>0} (n_i-mu_i)^2 / n_i.
 * @throws std::invalid_argument If model probabilities or selected-count
 *         bookkeeping are invalid.
 * @throws std::out_of_range If the selected raw range is unavailable.
 *
 * Expected counts are mu_i = N_fit p_i with N_fit equal to the exact selected
 * raw count. Bins with n_i == 0 are omitted, matching the historical SAFA
 * Neyman weighting. The model remains normalized over the complete selected
 * region; no free amplitude is introduced.
 */
[[nodiscard]] double binned_neyman_chi_square(
    const std::vector<std::uint64_t>& bins,
    std::size_t offset,
    const StatisticalRegion& region,
    const std::vector<double>& probabilities
);

/**
 * @brief Evaluate Pearson chi-square from raw counts and probabilities.
 * @param bins Slot-major raw histogram count storage.
 * @param offset First raw counter belonging to the logical histogram.
 * @param region Selected contiguous statistical region.
 * @param probabilities Unit-normalized model probabilities for selected bins.
 * @return chi2_P = sum_i (n_i-mu_i)^2 / mu_i.
 * @throws std::invalid_argument If model probabilities or selected-count
 *         bookkeeping are invalid.
 * @throws std::out_of_range If the selected raw range is unavailable.
 *
 * Zero-count observed bins remain in Pearson chi-square because mu_i is the
 * denominator. The normalization and expected-count convention are identical
 * to the Poisson and Neyman estimators.
 */
[[nodiscard]] double binned_pearson_chi_square(
    const std::vector<std::uint64_t>& bins,
    std::size_t offset,
    const StatisticalRegion& region,
    const std::vector<double>& probabilities
);

/**
 * @brief Convert exact-bin model probabilities into density values.
 * @param probabilities Probabilities over selected uniform bins.
 * @param binning Owning validated uniform binning.
 * @return Probability divided by exact bin width for every selected bin.
 * @throws std::invalid_argument If a probability is non-finite or negative.
 *
 * This is presentation state only. It uses the same integrated probabilities
 * as the likelihood and never evaluates a model at bin centers.
 */
[[nodiscard]] std::vector<double> probabilities_to_pdf(
    const std::vector<double>& probabilities,
    const HistogramBinningConfig& binning
);

}  // namespace hbt

#endif  // HBT_FITS_BINNED_MODELS_H
