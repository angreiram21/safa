/**
 * @file statistical_analysis.cpp
 * @brief Pure post-sample region, moment, and normalization operations.
 */

#include "hbt/fits/statistical_analysis.h"

#include <cmath>
#include <limits>
#include <stdexcept>

namespace hbt {
namespace {

/**
 * @brief Validate one logical slot lies inside flattened raw storage.
 * @param bins Flattened raw count storage.
 * @param offset First requested counter.
 * @param nbins Number of requested counters.
 * @throws std::out_of_range If the requested range is not fully present.
 */
void require_slot_range(
    const std::vector<std::uint64_t>& bins,
    std::size_t offset,
    std::size_t nbins
) {
    if (offset > bins.size() || nbins > bins.size() - offset) {
        throw std::out_of_range(
            "HBT statistical analysis: raw histogram slot is out of range"
        );
    }
}

/**
 * @brief Add one raw count with explicit uint64_t overflow detection.
 * @param sum Running selected count.
 * @param value Count to add.
 * @return Exact updated count.
 * @throws std::overflow_error If the addition is not representable.
 */
std::uint64_t checked_add(std::uint64_t sum, std::uint64_t value) {
    if (value > std::numeric_limits<std::uint64_t>::max() - sum) {
        throw std::overflow_error(
            "HBT statistical analysis: selected count overflow"
        );
    }
    return sum + value;
}

/**
 * @brief Sum one inclusive raw-histogram interval.
 * @param bins Flattened raw count storage.
 * @param offset First logical histogram counter.
 * @param first Inclusive logical bin index.
 * @param last Inclusive logical bin index.
 * @return Exact uint64_t selected count.
 * @throws std::overflow_error If the count sum is not representable.
 */
std::uint64_t sum_region(
    const std::vector<std::uint64_t>& bins,
    std::size_t offset,
    std::size_t first,
    std::size_t last
) {
    std::uint64_t sum = 0U;
    for (std::size_t bin = first; bin <= last; ++bin) {
        sum = checked_add(sum, bins[offset + bin]);
    }
    return sum;
}

/**
 * @brief Locate the first contiguous plateau of the global maximum.
 * @param bins Flattened raw count storage.
 * @param offset First logical histogram counter.
 * @param nbins Logical histogram bin count.
 * @return Inclusive left and right modal-plateau indices.
 * @throws std::logic_error If called for an all-zero histogram.
 */
std::pair<std::size_t, std::size_t> modal_plateau(
    const std::vector<std::uint64_t>& bins,
    std::size_t offset,
    std::size_t nbins
) {
    std::uint64_t maximum = 0U;
    for (std::size_t bin = 0U; bin < nbins; ++bin) {
        if (bins[offset + bin] > maximum) {
            maximum = bins[offset + bin];
        }
    }
    if (maximum == 0U) {
        throw std::logic_error(
            "HBT statistical analysis: modal plateau requested for empty "
            "histogram"
        );
    }

    std::size_t left = 0U;
    while (bins[offset + left] != maximum) {
        ++left;
    }

    std::size_t right = left;
    while (right + 1U < nbins &&
           bins[offset + right + 1U] == maximum) {
        ++right;
    }
    return {left, right};
}

/**
 * @brief Return monotonic PAVA levels for one contiguous interval.
 * @param bins Flattened raw count storage.
 * @param offset First logical histogram counter.
 * @param first Inclusive first bin to estimate.
 * @param last Inclusive last bin to estimate.
 * @param nondecreasing true for a non-decreasing estimate, false for a
 *        non-increasing estimate.
 * @return One monotonic long-double level per input bin.
 *
 * Adjacent blocks are pooled whenever their fitted means violate the requested
 * ordering. Raw counts are never modified; the returned levels are a
 * shape-constrained auxiliary estimate only.
 */
std::vector<long double> monotonic_pava(
    const std::vector<std::uint64_t>& bins,
    std::size_t offset,
    std::size_t first,
    std::size_t last,
    bool nondecreasing
) {
    struct Block {
        long double sum;
        std::size_t count;
    };

    std::vector<Block> blocks;
    blocks.reserve(last - first + 1U);
    for (std::size_t bin = first; bin <= last; ++bin) {
        blocks.push_back({
            static_cast<long double>(bins[offset + bin]),
            1U
        });
        while (blocks.size() >= 2U) {
            const Block& left = blocks[blocks.size() - 2U];
            const Block& right = blocks.back();
            const long double left_mean =
                left.sum / static_cast<long double>(left.count);
            const long double right_mean =
                right.sum / static_cast<long double>(right.count);
            const bool ordered = nondecreasing
                ? left_mean <= right_mean
                : left_mean >= right_mean;
            if (ordered) {
                break;
            }
            Block merged{
                left.sum + right.sum,
                left.count + right.count
            };
            blocks.pop_back();
            blocks.back() = merged;
        }
    }

    std::vector<long double> levels;
    levels.reserve(last - first + 1U);
    for (const Block& block : blocks) {
        const long double mean =
            block.sum / static_cast<long double>(block.count);
        for (std::size_t index = 0U; index < block.count; ++index) {
            levels.push_back(mean);
        }
    }
    return levels;
}

/**
 * @brief Return non-increasing PAVA levels for one contiguous interval.
 * @param bins Flattened raw count storage.
 * @param offset First logical histogram counter.
 * @param first Inclusive first bin to estimate.
 * @param last Inclusive last bin to estimate.
 * @return One non-increasing long-double level per input bin.
 */
std::vector<long double> nonincreasing_pava(
    const std::vector<std::uint64_t>& bins,
    std::size_t offset,
    std::size_t first,
    std::size_t last
) {
    return monotonic_pava(bins, offset, first, last, false);
}

/**
 * @brief Return non-decreasing PAVA levels for one contiguous interval.
 * @param bins Flattened raw count storage.
 * @param offset First logical histogram counter.
 * @param first Inclusive first bin to estimate.
 * @param last Inclusive last bin to estimate.
 * @return One non-decreasing long-double level per input bin.
 */
std::vector<long double> nondecreasing_pava(
    const std::vector<std::uint64_t>& bins,
    std::size_t offset,
    std::size_t first,
    std::size_t last
) {
    return monotonic_pava(bins, offset, first, last, true);
}

/**
 * @brief Squared-error contribution of one PAVA block.
 * @param sum Sum of raw counts in the block.
 * @param sum_squares Sum of squared raw counts in the block.
 * @param count Number of bins in the block.
 * @return Least-squares error around the block mean.
 */
long double pava_block_squared_error(
    long double sum,
    long double sum_squares,
    std::size_t count
) {
    const long double size = static_cast<long double>(count);
    const long double error = sum_squares - (sum * sum) / size;
    return error > 0.0L ? error : 0.0L;
}

/**
 * @brief Prefix errors of non-decreasing least-squares isotonic regressions.
 * @param bins Flattened raw count storage.
 * @param offset First logical histogram counter.
 * @param first Inclusive first bin.
 * @param last Inclusive last bin.
 * @return For each local index i, the PAVA squared error on [first,first+i].
 *
 * The PAVA stack is updated incrementally, so all prefix errors are obtained in
 * linear amortized time rather than by refitting every prefix independently.
 */
std::vector<long double> prefix_nondecreasing_pava_errors(
    const std::vector<std::uint64_t>& bins,
    std::size_t offset,
    std::size_t first,
    std::size_t last
) {
    struct Block {
        long double sum;
        long double sum_squares;
        std::size_t count;
    };

    std::vector<Block> blocks;
    blocks.reserve(last - first + 1U);
    std::vector<long double> errors;
    errors.reserve(last - first + 1U);
    long double total_error = 0.0L;

    for (std::size_t bin = first; bin <= last; ++bin) {
        const long double value = static_cast<long double>(bins[offset + bin]);
        blocks.push_back({value, value * value, 1U});
        while (blocks.size() >= 2U) {
            const Block& left = blocks[blocks.size() - 2U];
            const Block& right = blocks.back();
            const long double left_mean =
                left.sum / static_cast<long double>(left.count);
            const long double right_mean =
                right.sum / static_cast<long double>(right.count);
            if (left_mean <= right_mean) {
                break;
            }

            total_error -= pava_block_squared_error(
                left.sum, left.sum_squares, left.count
            );
            total_error -= pava_block_squared_error(
                right.sum, right.sum_squares, right.count
            );
            Block merged{
                left.sum + right.sum,
                left.sum_squares + right.sum_squares,
                left.count + right.count
            };
            blocks.pop_back();
            blocks.back() = merged;
            total_error += pava_block_squared_error(
                merged.sum, merged.sum_squares, merged.count
            );
        }
        errors.push_back(total_error);
    }
    return errors;
}

/**
 * @brief Suffix errors of non-increasing least-squares isotonic regressions.
 * @param bins Flattened raw count storage.
 * @param offset First logical histogram counter.
 * @param first Inclusive first bin.
 * @param last Inclusive last bin.
 * @return For each local index i, the PAVA squared error on [first+i,last].
 *
 * Reading the histogram from right to left converts a non-increasing suffix in
 * the original order into a non-decreasing prefix in reverse order.
 */
std::vector<long double> suffix_nonincreasing_pava_errors(
    const std::vector<std::uint64_t>& bins,
    std::size_t offset,
    std::size_t first,
    std::size_t last
) {
    struct Block {
        long double sum;
        long double sum_squares;
        std::size_t count;
    };

    const std::size_t size = last - first + 1U;
    std::vector<Block> blocks;
    blocks.reserve(size);
    std::vector<long double> errors(size, 0.0L);
    long double total_error = 0.0L;

    for (std::size_t reverse = 0U; reverse < size; ++reverse) {
        const std::size_t bin = last - reverse;
        const long double value = static_cast<long double>(bins[offset + bin]);
        blocks.push_back({value, value * value, 1U});
        while (blocks.size() >= 2U) {
            const Block& left = blocks[blocks.size() - 2U];
            const Block& right = blocks.back();
            const long double left_mean =
                left.sum / static_cast<long double>(left.count);
            const long double right_mean =
                right.sum / static_cast<long double>(right.count);
            if (left_mean <= right_mean) {
                break;
            }

            total_error -= pava_block_squared_error(
                left.sum, left.sum_squares, left.count
            );
            total_error -= pava_block_squared_error(
                right.sum, right.sum_squares, right.count
            );
            Block merged{
                left.sum + right.sum,
                left.sum_squares + right.sum_squares,
                left.count + right.count
            };
            blocks.pop_back();
            blocks.back() = merged;
            total_error += pava_block_squared_error(
                merged.sum, merged.sum_squares, merged.count
            );
        }
        errors[bin - first] = total_error;
    }
    return errors;
}

/**
 * @brief Build the least-squares unimodal PAVA estimate of a radial histogram.
 * @param bins Flattened raw count storage.
 * @param offset First logical histogram counter.
 * @param first Inclusive first bin.
 * @param last Inclusive last bin.
 * @return Shape-constrained levels and one optimal modal-bin index.
 *
 * Every candidate mode is ranked by the sum of the non-decreasing prefix PAVA
 * error and the non-increasing suffix PAVA error. Prefix and suffix errors are
 * accumulated in linear amortized time. The selected mode is then reconstructed
 * with one PAVA pass on each side. This suppresses isolated endpoint spikes
 * without smoothing or changing the raw histogram used by the fits.
 */
std::pair<std::vector<long double>, std::size_t> unimodal_pava(
    const std::vector<std::uint64_t>& bins,
    std::size_t offset,
    std::size_t first,
    std::size_t last
) {
    const std::vector<long double> prefix_errors =
        prefix_nondecreasing_pava_errors(bins, offset, first, last);
    const std::vector<long double> suffix_errors =
        suffix_nonincreasing_pava_errors(bins, offset, first, last);

    std::size_t mode = first;
    long double best_error = prefix_errors.front() + suffix_errors.front();
    for (std::size_t local = 1U; local < prefix_errors.size(); ++local) {
        const long double error = prefix_errors[local] + suffix_errors[local];
        if (error < best_error) {
            best_error = error;
            mode = first + local;
        }
    }

    const std::vector<long double> left =
        nondecreasing_pava(bins, offset, first, mode);
    const std::vector<long double> right =
        nonincreasing_pava(bins, offset, mode, last);
    const long double raw_mode = static_cast<long double>(bins[offset + mode]);

    std::vector<long double> levels;
    levels.reserve(last - first + 1U);
    for (std::size_t local = 0U; local + 1U < left.size(); ++local) {
        levels.push_back(left[local]);
    }
    levels.push_back(raw_mode);
    for (std::size_t local = 1U; local < right.size(); ++local) {
        levels.push_back(right[local]);
    }
    return {levels, mode};
}

/**
 * @brief Linearly interpolate one half-height crossing between bin centers.
 * @param x0 Center of the first bracketing bin.
 * @param y0 Count in the first bracketing bin.
 * @param x1 Center of the second bracketing bin.
 * @param y1 Count in the second bracketing bin.
 * @param half Target half-maximum count level.
 * @return Finite crossing coordinate, or std::nullopt for a degenerate pair.
 */
std::optional<double> interpolate_half_crossing(
    double x0,
    long double y0,
    double x1,
    long double y1,
    long double half
) {
    const long double denominator = y1 - y0;
    if (denominator == 0.0L) {
        return std::nullopt;
    }
    const long double fraction = (half - y0) / denominator;
    if (fraction < 0.0L || fraction > 1.0L) {
        return std::nullopt;
    }
    const long double crossing = static_cast<long double>(x0) + fraction *
        (static_cast<long double>(x1) - static_cast<long double>(x0));
    const double value = static_cast<double>(crossing);
    if (!std::isfinite(value)) {
        return std::nullopt;
    }
    return value;
}

}  // namespace

std::optional<StatisticalRegion> select_shape_region(
    const std::vector<std::uint64_t>& bins,
    std::size_t offset,
    const HistogramBinningConfig& binning
) {
    require_slot_range(bins, offset, binning.nbins);

    bool any = false;
    for (std::size_t bin = 0U; bin < binning.nbins; ++bin) {
        if (bins[offset + bin] != 0U) {
            any = true;
            break;
        }
    }
    if (!any) {
        return std::nullopt;
    }

    const std::size_t last = binning.nbins - 1U;
    return StatisticalRegion{0U, last, sum_region(bins, offset, 0U, last)};
}

std::optional<StatisticalRegion> select_gaussian_core_region(
    FitObservableFamily family,
    const std::vector<std::uint64_t>& bins,
    std::size_t offset,
    const HistogramBinningConfig& binning,
    const StatisticalRegion& full_region,
    double threshold_fraction
) {
    if (!std::isfinite(threshold_fraction) ||
        threshold_fraction <= 0.0 || threshold_fraction >= 1.0) {
        throw std::invalid_argument(
            "HBT Gaussian core selection: threshold fraction must be in (0,1)"
        );
    }
    require_slot_range(bins, offset, binning.nbins);
    if (full_region.first_bin != 0U ||
        full_region.last_bin >= binning.nbins ||
        full_region.selected_count == 0U) {
        throw std::out_of_range(
            "HBT Gaussian core selection: full shape region is invalid"
        );
    }

    std::size_t branch_first = 0U;
    std::size_t branch_last = full_region.last_bin;
    long double reference = 0.0L;

    switch (family) {
        case FitObservableFamily::Radial: {
            const auto mode = modal_plateau(bins, offset, binning.nbins);
            if (mode.second > branch_last) {
                return std::nullopt;
            }
            branch_first = mode.second;
            reference = static_cast<long double>(bins[offset + mode.second]);
            break;
        }
        case FitObservableFamily::OSL: {
            for (std::size_t bin = 0U; bin <= full_region.last_bin; ++bin) {
                if (bins[offset + bin] == 0U) {
                    if (bin == 0U) {
                        return std::nullopt;
                    }
                    branch_last = bin - 1U;
                    break;
                }
            }
            if (branch_last < branch_first) {
                return std::nullopt;
            }
            break;
        }
    }

    const std::vector<long double> levels = nonincreasing_pava(
        bins,
        offset,
        branch_first,
        branch_last
    );
    if (levels.empty()) {
        return std::nullopt;
    }
    if (family == FitObservableFamily::OSL) {
        reference = levels.front();
    }
    if (!(reference > 0.0L)) {
        return std::nullopt;
    }

    const long double threshold =
        static_cast<long double>(threshold_fraction) * reference;
    std::size_t last = branch_last;
    bool crossed = false;
    for (std::size_t local = 0U; local < levels.size(); ++local) {
        const std::size_t bin = branch_first + local;
        if (levels[local] <= threshold) {
            if (bin == 0U) {
                return std::nullopt;
            }
            last = bin - 1U;
            crossed = true;
            break;
        }
    }
    if (!crossed) {
        last = branch_last;
    }
    if (last >= binning.nbins) {
        return std::nullopt;
    }
    return StatisticalRegion{0U, last, sum_region(bins, offset, 0U, last)};
}

std::optional<double> half_maximum_radius_seed(
    FitObservableFamily family,
    const std::vector<std::uint64_t>& bins,
    std::size_t offset,
    const HistogramBinningConfig& binning,
    const StatisticalRegion& full_region
) {
    require_slot_range(bins, offset, binning.nbins);
    if (full_region.first_bin != 0U ||
        full_region.last_bin >= binning.nbins ||
        full_region.selected_count == 0U) {
        throw std::out_of_range(
            "HBT half-maximum seed: full shape region is invalid"
        );
    }

    constexpr double kSqrtLn2 = 0.83255461115769775635;
    constexpr double kRadialFwhmOverRadius = 2.3098847205021675;

    switch (family) {
        case FitObservableFamily::OSL: {
            const std::vector<long double> levels = nonincreasing_pava(
                bins,
                offset,
                full_region.first_bin,
                full_region.last_bin
            );
            if (levels.empty() || !(levels.front() > 0.0L)) {
                return std::nullopt;
            }
            const long double half = 0.5L * levels.front();
            std::optional<double> right_crossing;
            for (std::size_t local = 1U; local < levels.size(); ++local) {
                if (levels[local] <= half) {
                    const std::size_t previous = local - 1U;
                    right_crossing = interpolate_half_crossing(
                        histogram_bin_center(
                            binning,
                            full_region.first_bin + previous
                        ),
                        levels[previous],
                        histogram_bin_center(
                            binning,
                            full_region.first_bin + local
                        ),
                        levels[local],
                        half
                    );
                    break;
                }
            }
            if (!right_crossing.has_value()) {
                return std::nullopt;
            }

            // The stored observable is |x|. Mirroring the positive crossing
            // gives FWHM = 2*x_half for a zero-centered Gaussian.
            const double fwhm = 2.0 * right_crossing.value();
            const double radius = fwhm / (4.0 * kSqrtLn2);
            return std::isfinite(radius) && radius > 0.0
                ? std::optional<double>{radius}
                : std::nullopt;
        }
        case FitObservableFamily::Radial: {
            const auto estimate = unimodal_pava(
                bins,
                offset,
                full_region.first_bin,
                full_region.last_bin
            );
            const std::vector<long double>& levels = estimate.first;
            if (levels.empty()) {
                return std::nullopt;
            }
            const std::size_t mode_local =
                estimate.second - full_region.first_bin;
            if (mode_local >= levels.size() || !(levels[mode_local] > 0.0L)) {
                return std::nullopt;
            }

            const long double maximum = levels[mode_local];
            const long double half = 0.5L * maximum;
            std::size_t mode_first = mode_local;
            while (mode_first > 0U &&
                   levels[mode_first - 1U] == maximum) {
                --mode_first;
            }
            std::size_t mode_last = mode_local;
            while (mode_last + 1U < levels.size() &&
                   levels[mode_last + 1U] == maximum) {
                ++mode_last;
            }
            if (mode_first == 0U) {
                return std::nullopt;
            }

            std::optional<double> left_crossing;
            for (std::size_t local = mode_first; local > 0U; --local) {
                const std::size_t left = local - 1U;
                if (levels[left] <= half) {
                    left_crossing = interpolate_half_crossing(
                        histogram_bin_center(
                            binning,
                            full_region.first_bin + left
                        ),
                        levels[left],
                        histogram_bin_center(
                            binning,
                            full_region.first_bin + local
                        ),
                        levels[local],
                        half
                    );
                    break;
                }
            }

            std::optional<double> right_crossing;
            for (std::size_t local = mode_last + 1U;
                 local < levels.size();
                 ++local) {
                if (levels[local] <= half) {
                    const std::size_t previous = local - 1U;
                    right_crossing = interpolate_half_crossing(
                        histogram_bin_center(
                            binning,
                            full_region.first_bin + previous
                        ),
                        levels[previous],
                        histogram_bin_center(
                            binning,
                            full_region.first_bin + local
                        ),
                        levels[local],
                        half
                    );
                    break;
                }
            }
            if (!left_crossing.has_value() ||
                !right_crossing.has_value() ||
                right_crossing.value() <= left_crossing.value()) {
                return std::nullopt;
            }

            const double fwhm =
                right_crossing.value() - left_crossing.value();
            const double radius = fwhm / kRadialFwhmOverRadius;
            return std::isfinite(radius) && radius > 0.0
                ? std::optional<double>{radius}
                : std::nullopt;
        }
    }
    return std::nullopt;
}

std::optional<StatisticalRegion> select_delta_t_region(
    const std::vector<std::uint64_t>& bins,
    std::size_t offset,
    const HistogramBinningConfig& binning
) {
    require_slot_range(bins, offset, binning.nbins);

    bool any = false;
    for (std::size_t bin = 0U; bin < binning.nbins; ++bin) {
        if (bins[offset + bin] != 0U) {
            any = true;
            break;
        }
    }
    if (!any) {
        return std::nullopt;
    }

    const std::size_t first = 0U;
    const std::size_t last = binning.nbins - 1U;
    return StatisticalRegion{
        first,
        last,
        sum_region(bins, offset, first, last)
    };
}

double histogram_bin_lower_edge(
    const HistogramBinningConfig& binning,
    std::size_t bin_index
) {
    if (bin_index >= binning.nbins) {
        throw std::out_of_range(
            "HBT statistical analysis: bin index is out of range"
        );
    }
    const double width = 1.0 / binning.inverse_bin_width;
    return binning.minimum + static_cast<double>(bin_index) * width;
}

double histogram_bin_upper_edge(
    const HistogramBinningConfig& binning,
    std::size_t bin_index
) {
    if (bin_index >= binning.nbins) {
        throw std::out_of_range(
            "HBT statistical analysis: bin index is out of range"
        );
    }
    const double width = 1.0 / binning.inverse_bin_width;
    return binning.minimum + static_cast<double>(bin_index + 1U) * width;
}

double histogram_bin_center(
    const HistogramBinningConfig& binning,
    std::size_t bin_index
) {
    const double lower = histogram_bin_lower_edge(binning, bin_index);
    const double upper = histogram_bin_upper_edge(binning, bin_index);
    return 0.5 * (lower + upper);
}

std::vector<NormalizedHistogramBin> normalize_region(
    const std::vector<std::uint64_t>& bins,
    std::size_t offset,
    const HistogramBinningConfig& binning,
    const StatisticalRegion& region
) {
    require_slot_range(bins, offset, binning.nbins);
    if (region.first_bin > region.last_bin ||
        region.last_bin >= binning.nbins ||
        region.selected_count == 0U) {
        throw std::invalid_argument(
            "HBT analysis normalization: invalid selected region"
        );
    }
    const std::uint64_t count = sum_region(
        bins,
        offset,
        region.first_bin,
        region.last_bin
    );
    if (count != region.selected_count) {
        throw std::invalid_argument(
            "HBT analysis normalization: selected count differs from raw counts"
        );
    }

    const double width = 1.0 / binning.inverse_bin_width;
    const double denominator =
        static_cast<double>(region.selected_count) * width;
    std::vector<NormalizedHistogramBin> result;
    result.reserve(region.last_bin - region.first_bin + 1U);

    for (std::size_t bin = region.first_bin;
         bin <= region.last_bin;
         ++bin) {
        const std::uint64_t raw = bins[offset + bin];
        const double lower = histogram_bin_lower_edge(binning, bin);
        const double upper = histogram_bin_upper_edge(binning, bin);
        result.push_back({
            bin,
            lower,
            upper,
            0.5 * (lower + upper),
            static_cast<double>(raw) / denominator,
            std::sqrt(static_cast<double>(raw)) / denominator
        });
    }
    return result;
}

DeltaTHistogramResult calculate_delta_t_statistics(
    const std::vector<std::uint64_t>& bins,
    std::size_t offset,
    const HistogramBinningConfig& binning,
    const StatisticalRegion& region
) {
    require_slot_range(bins, offset, binning.nbins);
    if (region.first_bin > region.last_bin ||
        region.last_bin >= binning.nbins ||
        region.selected_count == 0U) {
        throw std::invalid_argument(
            "HBT analysis delta_t: invalid selected region"
        );
    }
    const std::uint64_t count = sum_region(
        bins,
        offset,
        region.first_bin,
        region.last_bin
    );
    if (count != region.selected_count) {
        throw std::invalid_argument(
            "HBT analysis delta_t: selected count differs from raw counts"
        );
    }

    double weighted_sum = 0.0;
    double weighted_square_sum = 0.0;
    for (std::size_t bin = region.first_bin;
         bin <= region.last_bin;
         ++bin) {
        const double center = histogram_bin_center(binning, bin);
        const double weight = static_cast<double>(bins[offset + bin]);
        weighted_sum += weight * center;
        weighted_square_sum += weight * center * center;
    }

    const double n = static_cast<double>(region.selected_count);
    const double mean = weighted_sum / n;
    const double variance = weighted_square_sum / n - mean * mean;

    DeltaTHistogramResult result{
        region,
        region.selected_count,
        {},
        DeltaTStatisticsStatus::Valid,
        mean,
        std::nullopt,
        std::nullopt,
        std::nullopt
    };

    if (!std::isfinite(variance) || variance < 0.0) {
        result.status = DeltaTStatisticsStatus::InvalidVariance;
        result.mean = std::isfinite(mean)
            ? std::optional<double>{mean}
            : std::nullopt;
        return result;
    }

    const double sigma = std::sqrt(variance);
    result.sigma = sigma;
    if (region.selected_count <= 1U) {
        result.status = DeltaTStatisticsStatus::InsufficientCount;
        return result;
    }

    const double mean_error = sigma /
        std::sqrt(static_cast<double>(region.selected_count));
    const double sigma_error = sigma /
        std::sqrt(2.0 * static_cast<double>(region.selected_count - 1U));
    if (!std::isfinite(mean) || !std::isfinite(mean_error) ||
        !std::isfinite(sigma) || !std::isfinite(sigma_error)) {
        result.status = DeltaTStatisticsStatus::InvalidVariance;
        result.mean = std::isfinite(mean)
            ? std::optional<double>{mean}
            : std::nullopt;
        result.sigma = std::isfinite(sigma)
            ? std::optional<double>{sigma}
            : std::nullopt;
        return result;
    }

    result.mean_error = mean_error;
    result.sigma_error = sigma_error;
    return result;
}

}  // namespace hbt
