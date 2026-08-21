/**
 * @file statistical_analysis.h
 * @brief Pure post-sample region, moment, and normalization operations.
 */

#ifndef HBT_FITS_STATISTICAL_ANALYSIS_H
#define HBT_FITS_STATISTICAL_ANALYSIS_H

#include "hbt/config/hbt_config.h"
#include "hbt/fits/fit_results.h"

#include <cstddef>
#include <cstdint>
#include <optional>
#include <vector>

namespace hbt {

/**
 * @brief Select the OSL/radial region from zero through the right-tail cut.
 * @param bins Slot-major raw histogram count storage.
 * @param offset First raw counter belonging to the logical histogram.
 * @param binning Validated uniform binning for the logical histogram.
 * @return Selected region, or std::nullopt when the histogram is empty.
 * @throws std::out_of_range If the requested logical histogram is absent.
 * @throws std::overflow_error If the selected uint64_t count sum overflows.
 *
 * Leading empty bins are retained. The first empty bin strictly to the right
 * of the right edge of the first contiguous global-maximum plateau is excluded
 * and fixes the tail cut. No later island is reincorporated.
 */
[[nodiscard]] std::optional<StatisticalRegion> select_shape_region(
    const std::vector<std::uint64_t>& bins,
    std::size_t offset,
    const HistogramBinningConfig& binning
);

/**
 * @brief Select the compact Gaussian-core fit region inside a shape region.
 * @param family OSL or radial physical model family.
 * @param bins Slot-major raw histogram count storage.
 * @param offset First raw counter belonging to the logical histogram.
 * @param binning Validated uniform binning for the logical histogram.
 * @param full_region Existing full statistical region retained by the mixed fit.
 * @param threshold_fraction Positive fraction of the monotonic core reference
 *        level at which the first excluded bin is identified. Production uses
 *        0.10; 0.05 is retained for systematic studies.
 * @return Core region starting at bin zero, or std::nullopt when a core cannot
 *         be defined from the supplied counts.
 * @throws std::invalid_argument If @p threshold_fraction is outside (0,1).
 * @throws std::out_of_range If the logical histogram or region is unavailable.
 * @throws std::overflow_error If the selected uint64_t count sum overflows.
 *
 * The selector never changes fit counts. It applies a non-increasing PAVA
 * estimate only to decide the upper edge of the Gaussian fit. For radial
 * histograms the estimate starts at the right edge of the first contiguous
 * global-maximum plateau. For OSL histograms it starts at bin zero because the
 * absolute-coordinate Gaussian has its physical maximum at x=0. The first
 * threshold-crossing bin is excluded. The full region remains unchanged and
 * is still used by the mixed model.
 */
[[nodiscard]] std::optional<StatisticalRegion> select_gaussian_core_region(
    FitObservableFamily family,
    const std::vector<std::uint64_t>& bins,
    std::size_t offset,
    const HistogramBinningConfig& binning,
    const StatisticalRegion& full_region,
    double threshold_fraction = 0.10
);


/**
 * @brief Derive a Gaussian-radius seed from the histogram half-maximum width.
 * @param family OSL or radial physical model family.
 * @param bins Slot-major raw histogram count storage.
 * @param offset First raw counter belonging to the logical histogram.
 * @param binning Validated uniform binning for the logical histogram.
 * @param full_region Existing full selected shape region.
 * @return Radius seed in the model's R parameterization, or std::nullopt when
 *         the required half-maximum crossing(s) cannot be identified.
 * @throws std::out_of_range If the logical histogram or region is unavailable.
 *
 * The half-height locations are linearly interpolated between neighboring bin
 * centers. OSL histograms contain |x|, so the observed right half-maximum
 * location is mirrored about zero to obtain FWHM and
 * R = FWHM/(4 sqrt(ln 2)). Radial histograms use both measured crossings around
 * the non-zero mode and R = FWHM/2.3098847205021675, the exact numerical width
 * ratio of r^2 exp[-r^2/(4R^2)]. Raw counts and N_selected are never modified.
 */
[[nodiscard]] std::optional<double> half_maximum_radius_seed(
    FitObservableFamily family,
    const std::vector<std::uint64_t>& bins,
    std::size_t offset,
    const HistogramBinningConfig& binning,
    const StatisticalRegion& full_region
);

/**
 * @brief Select the signed delta-t region around its modal plateau.
 * @param bins Slot-major raw histogram count storage.
 * @param offset First raw counter belonging to the logical histogram.
 * @param binning Validated signed delta-t binning.
 * @return Selected region, or std::nullopt when the histogram is empty.
 * @throws std::out_of_range If the requested logical histogram is absent.
 * @throws std::overflow_error If the selected uint64_t count sum overflows.
 *
 * The first empty bin on each side of the modal plateau is excluded and ends
 * that side of the selected region.
 */
[[nodiscard]] std::optional<StatisticalRegion> select_delta_t_region(
    const std::vector<std::uint64_t>& bins,
    std::size_t offset,
    const HistogramBinningConfig& binning
);

/**
 * @brief Build the final normalized distribution for one selected region.
 * @param bins Slot-major raw histogram count storage.
 * @param offset First raw counter belonging to the logical histogram.
 * @param binning Validated uniform histogram binning.
 * @param region Selected contiguous statistical region.
 * @return One normalized result per selected raw bin.
 * @throws std::invalid_argument If selected_count is zero or inconsistent.
 * @throws std::out_of_range If the requested region is outside raw storage.
 * @throws std::overflow_error If re-summing selected counts overflows.
 *
 * This operation is for presentation state only and must be called after fits
 * or delta-t statistics. It does not modify or retain the raw counts.
 */
[[nodiscard]] std::vector<NormalizedHistogramBin> normalize_region(
    const std::vector<std::uint64_t>& bins,
    std::size_t offset,
    const HistogramBinningConfig& binning,
    const StatisticalRegion& region
);

/**
 * @brief Compute signed delta-t moments from raw selected counts.
 * @param bins Slot-major raw histogram count storage.
 * @param offset First raw counter belonging to the logical histogram.
 * @param binning Validated signed delta-t binning.
 * @param region Selected peak-centered region.
 * @return Delta-t result without its normalized presentation bins.
 * @throws std::invalid_argument If region.selected_count is zero or differs
 *         from the selected raw count sum.
 * @throws std::out_of_range If the region lies outside raw storage.
 * @throws std::overflow_error If re-summing selected counts overflows.
 *
 * The population variance is evaluated exactly from weighted bin centers. A
 * negative or non-finite variance is reported and is never clamped to zero.
 */
[[nodiscard]] DeltaTHistogramResult calculate_delta_t_statistics(
    const std::vector<std::uint64_t>& bins,
    std::size_t offset,
    const HistogramBinningConfig& binning,
    const StatisticalRegion& region
);

/**
 * @brief Return the center of one validated uniform histogram bin.
 * @param binning Validated uniform histogram binning.
 * @param bin_index Zero-based bin index.
 * @return Arithmetic bin center in the family's physical units.
 * @throws std::out_of_range If @p bin_index is outside the binning.
 */
[[nodiscard]] double histogram_bin_center(
    const HistogramBinningConfig& binning,
    std::size_t bin_index
);

/**
 * @brief Return the lower edge of one validated uniform histogram bin.
 * @param binning Validated uniform histogram binning.
 * @param bin_index Zero-based bin index.
 * @return Exact lower bin edge in the family's physical units.
 * @throws std::out_of_range If @p bin_index is outside the binning.
 */
[[nodiscard]] double histogram_bin_lower_edge(
    const HistogramBinningConfig& binning,
    std::size_t bin_index
);

/**
 * @brief Return the upper edge of one validated uniform histogram bin.
 * @param binning Validated uniform histogram binning.
 * @param bin_index Zero-based bin index.
 * @return Exact upper bin edge in the family's physical units.
 * @throws std::out_of_range If @p bin_index is outside the binning.
 */
[[nodiscard]] double histogram_bin_upper_edge(
    const HistogramBinningConfig& binning,
    std::size_t bin_index
);

}  // namespace hbt

#endif  // HBT_FITS_STATISTICAL_ANALYSIS_H
