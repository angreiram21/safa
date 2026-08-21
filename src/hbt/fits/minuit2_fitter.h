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
 * @brief Fit one normalized pure-Gaussian estimator with MIGRAD and MINOS.
 * @param family OSL or radial physical model family.
 * @param bins Slot-major raw uint64_t histogram count storage.
 * @param offset First raw counter belonging to the logical histogram.
 * @param binning Validated uniform binning for this histogram family.
 * @param region Validated 10%-core statistical region used by the Gaussian.
 * @param estimator Statistical objective minimized by this independent fit.
 * @param half_maximum_seed Gaussian R seed obtained from the half-maximum width
 *        of the full selected histogram and converted to model R units.
 * @return Complete fit result including both start diagnostics and MINOS.
 * @throws std::out_of_range If the selected raw slot is unavailable.
 *
 * The two deterministic starts are the moment-derived Gaussian radius and
 * @p half_maximum_seed. Both minimize the same estimator independently; the
 * valid minimum with the smallest q is selected. The only free parameter is
 * log(R), so R remains strictly positive without arbitrary physical bounds.
 */
[[nodiscard]] GaussianFitResult fit_gaussian_model(
    FitObservableFamily family,
    const std::vector<std::uint64_t>& bins,
    std::size_t offset,
    const HistogramBinningConfig& binning,
    const StatisticalRegion& region,
    FitEstimator estimator,
    double half_maximum_seed
);

/**
 * @brief Fit one normalized Gaussian-plus-exponential estimator.
 * @param family OSL or radial physical model family.
 * @param bins Slot-major raw uint64_t histogram count storage.
 * @param offset First raw counter belonging to the logical histogram.
 * @param binning Validated uniform binning for this histogram family.
 * @param region Full selected statistical region used by the mixed model.
 * @param estimator Statistical objective minimized by this independent fit.
 * @param gaussian_result Valid 10%-core Gaussian result from the same estimator;
 *        its fitted radius provides the R_G member of the core-seed set.
 * @param half_maximum_seed Gaussian R seed converted from histogram FWHM.
 * @return Complete 36-start fit result and explicit diagnostics.
 * @throws std::out_of_range If the selected raw slot is unavailable.
 *
 * The deterministic Cartesian product is
 * R_core={R_G,0.5R_HM,R_HM,2R_HM},
 * R_tail={0.5,1,2}R_tail,mom, and f_core={0.25,0.50,0.75}. Valid minima are
 * grouped into numerical basins. The basin whose geometric-mean R_core is
 * closest to R_HM in logarithmic relative scale is identified as the physical
 * Gaussian-core basin; the smallest q inside that basin is selected for MINOS.
 * Basin multiplicity is diagnostic only, and no ordering of R_core and R_tail
 * is imposed.
 */
[[nodiscard]] MixedFitResult fit_mixed_model(
    FitObservableFamily family,
    const std::vector<std::uint64_t>& bins,
    std::size_t offset,
    const HistogramBinningConfig& binning,
    const StatisticalRegion& region,
    FitEstimator estimator,
    const GaussianFitResult& gaussian_result,
    double half_maximum_seed
);

}  // namespace hbt

#endif  // HBT_FITS_MINUIT2_FITTER_H
