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

    const auto mode = modal_plateau(bins, offset, binning.nbins);
    std::size_t last = binning.nbins - 1U;
    for (std::size_t bin = mode.second + 1U;
         bin < binning.nbins;
         ++bin) {
        if (bins[offset + bin] == 0U) {
            last = bin - 1U;
            break;
        }
    }

    return StatisticalRegion{0U, last, sum_region(bins, offset, 0U, last)};
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

    const auto mode = modal_plateau(bins, offset, binning.nbins);
    std::size_t first = 0U;
    std::size_t last = binning.nbins - 1U;

    for (std::size_t bin = mode.first; bin > 0U; --bin) {
        if (bins[offset + bin - 1U] == 0U) {
            first = bin;
            break;
        }
    }
    for (std::size_t bin = mode.second + 1U;
         bin < binning.nbins;
         ++bin) {
        if (bins[offset + bin] == 0U) {
            last = bin - 1U;
            break;
        }
    }

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

    const double sigma_error = sigma /
        std::sqrt(2.0 * static_cast<double>(region.selected_count - 1U));
    if (!std::isfinite(mean) || !std::isfinite(sigma) ||
        !std::isfinite(sigma_error)) {
        result.status = DeltaTStatisticsStatus::InvalidVariance;
        result.mean = std::isfinite(mean)
            ? std::optional<double>{mean}
            : std::nullopt;
        result.sigma = std::isfinite(sigma)
            ? std::optional<double>{sigma}
            : std::nullopt;
        return result;
    }

    result.sigma_error = sigma_error;
    return result;
}

}  // namespace hbt
