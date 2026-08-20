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
 * @param gaussian_result Valid truncated Gaussian-core result from the same
 *        histogram and observable geometry. Its radius seeds every core start.
 * @return Complete consensus-multistart fit result and explicit diagnostics.
 * @throws std::out_of_range If the selected raw slot is unavailable.
 *
 * Five deterministic starts share R_core(0)=R_G^core from
 * @p gaussian_result. Their R_tail and f_core seeds are (R_tail,mom,0.50),
 * (0.5 R_tail,mom,0.50), (2 R_tail,mom,0.50), (R_tail,mom,0.25), and
 * (R_tail,mom,0.75). A result is publishable only when at least four starts
 * converge to the same numerical basin. Q selects the best realization only
 * inside that consensus basin; no ordering of R_core and R_tail is imposed.
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
