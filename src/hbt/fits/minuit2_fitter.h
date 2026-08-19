/**
 * @file minuit2_fitter.h
 * @brief Minuit2 MIGRAD/MINOS fits for post-sample HBT shape histograms.
 */

#ifndef HBT_FITS_MINUIT2_FITTER_H
#define HBT_FITS_MINUIT2_FITTER_H

#include "hbt/config/hbt_config.h"
#include "hbt/fits/fit_results.h"

#include <cstddef>
#include <cstdint>
#include <optional>
#include <vector>

namespace hbt {

/**
 * @brief Fit the normalized pure-Gaussian shape model with MIGRAD and MINOS.
 * @param family OSL or radial physical model family.
 * @param bins Slot-major raw uint64_t histogram count storage.
 * @param offset First raw counter belonging to the logical histogram.
 * @param binning Validated uniform binning for this histogram family.
 * @param region Selected contiguous statistical region.
 * @return Complete fit result including explicit failure diagnostics.
 * @throws std::out_of_range If the selected raw slot is unavailable.
 *
 * The only free parameter is log(R), which makes R strictly positive without
 * arbitrary physical bounds. The objective uses raw counts and exact-bin
 * probabilities. MINOS is run only after a valid MIGRAD minimum.
 */
[[nodiscard]] GaussianFitResult fit_gaussian_model(
    FitObservableFamily family,
    const std::vector<std::uint64_t>& bins,
    std::size_t offset,
    const HistogramBinningConfig& binning,
    const StatisticalRegion& region
);

/**
 * @brief Fit the normalized Gaussian-plus-exponential model.
 * @param family OSL or radial physical model family.
 * @param bins Slot-major raw uint64_t histogram count storage.
 * @param offset First raw counter belonging to the logical histogram.
 * @param binning Validated uniform binning for this histogram family.
 * @param region Selected contiguous statistical region.
 * @param gaussian_result Independent pure-Gaussian result for optional start B.
 * @return Complete multistart fit result and explicit diagnostics.
 * @throws std::out_of_range If the selected raw slot is unavailable.
 *
 * Start A derives both radius scales from raw moments. Start B exists only
 * when @p gaussian_result is fully valid; it reuses its fitted R for the core
 * and retains the data-derived tail seed. MIGRAD runs independently for each
 * start. MINOS runs only on the valid minimum with the smallest objective.
 */
[[nodiscard]] MixedFitResult fit_mixed_model(
    FitObservableFamily family,
    const std::vector<std::uint64_t>& bins,
    std::size_t offset,
    const HistogramBinningConfig& binning,
    const StatisticalRegion& region,
    const GaussianFitResult& gaussian_result
);

}  // namespace hbt

#endif  // HBT_FITS_MINUIT2_FITTER_H
