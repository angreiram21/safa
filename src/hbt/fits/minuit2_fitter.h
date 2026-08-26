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
 * @brief Fit one free-amplitude pure-Gaussian estimator with MIGRAD and MINOS.
 * @param family OSL or radial physical model family.
 * @param bins Slot-major raw uint64_t histogram count storage.
 * @param offset First raw counter belonging to the logical histogram.
 * @param binning Validated uniform binning for this histogram family.
 * @param region Full contiguous shape region used by the Gaussian.
 * @param estimator Statistical objective minimized by this independent fit.
 * @param half_maximum_seed Gaussian R seed obtained from the half-maximum width
 *        of the full selected histogram and converted to model R units.
 * @return Complete fit result including both start diagnostics and MINOS.
 * @throws std::out_of_range If the selected raw slot is unavailable.
 *
 * The two deterministic starts are the moment-derived Gaussian radius and
 * @p half_maximum_seed. Both minimize the same estimator independently; the
 * valid minimum with the smallest q is selected. The free parameters are
 * log(R) and log(A_G), so both radius and Gaussian amplitude remain strictly
 * positive without arbitrary physical bounds. Expected counts are evaluated
 * directly as N_selected*A_G times the exact Gaussian bin integral; the pure
 * Gaussian is not normalized to unit probability over the fit region. The
 * initial A_G at each radius start is the inverse full-region Gaussian
 * integral, reproducing the former normalized shape only as an initial state.
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
 * @param gaussian_result Valid full-range free-amplitude Gaussian result from the same estimator;
 *        its fitted radius provides the R_G member of the core-seed set.
 * @param half_maximum_seed Gaussian R seed converted from histogram FWHM.
 * @param core_fraction_policy Origin-dependent physical basin-fraction policy.
 * @return Complete 36-start fit result, terminal start endpoints, and explicit
 *         diagnostics.
 * @throws std::out_of_range If the selected raw slot is unavailable.
 *
 * The deterministic Cartesian product is
 * R_core={R_G,0.5R_HM,R_HM,2R_HM},
 * R_tail={0.5,1,2}R_tail,mom, and f_core={0.25,0.50,0.75}. Valid minima are
 * grouped into numerical basins. Basin fractions are filtered by the supplied
 * origin policy: PRD requires 0.1 < mean(f_core) < 0.9, while P and PR require
 * 0.1 < mean(f_core) < 0.99 to reject the near-pure-Gaussian degeneracy. Among
 * the remaining basins, every origin selects the basin reached by the largest
 * number of converged deterministic starts. Equal-size basins are ranked by
 * their smallest q and then by lowest start index. Once the basin is fixed,
 * its endpoints are ordered by increasing terminal q. Each endpoint supplies
 * coordinates for a fresh post-selection MIGRAD pass followed by MINOS. The
 * first fully publishable result is accepted; a MIGRAD or MINOS failure retries
 * the next endpoint only within that same basin. If all endpoints fail, the
 * lowest-q endpoint remains the primary failure diagnostic. These fallback
 * passes do not participate in basin selection. R_HM remains in the
 * deterministic core-seed set but does not rank final basins. If every basin is degenerate, the mixed
 * fit is invalidated with DegenerateCoreFraction.
 * The terminal physical R_core, R_tail and f_core coordinates of every start
 * are retained for post-run basin inspection, including finite endpoints from
 * starts that fail the acceptance contract. These endpoint diagnostics do not
 * participate in basin selection or fit validity. Basin multiplicity is
 * diagnostic only after basin selection, and no ordering of R_core and R_tail
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
    double half_maximum_seed,
    MixedCoreFractionPolicy core_fraction_policy
);

}  // namespace hbt

#endif  // HBT_FITS_MINUIT2_FITTER_H
