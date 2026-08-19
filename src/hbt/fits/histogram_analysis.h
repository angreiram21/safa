/**
 * @file histogram_analysis.h
 * @brief Complete post-sample HBT statistical orchestration.
 */

#ifndef HBT_FITS_HISTOGRAM_ANALYSIS_H
#define HBT_FITS_HISTOGRAM_ANALYSIS_H

#include "hbt/config/hbt_config.h"
#include "hbt/fits/fit_results.h"
#include "hbt/histograms/raw_histograms.h"

namespace hbt {

/**
 * @brief Build post-sample HBT results from completed raw counts.
 * @param config Validated HBT configuration that owns histogram metadata.
 * @param raw Complete-sample immutable raw histogram state.
 * @return Separate derived state aligned with product/origin/global/slice
 *         indices already present in @p raw and @p config.
 * @throws std::logic_error If raw state dimensions disagree with config.
 * @throws std::invalid_argument If validated metadata is internally invalid.
 * @throws std::overflow_error If selected raw count accumulation overflows.
 *
 * This function is synchronous. It does not retain references to @p config or
 * @p raw and never modifies raw counts. It performs no event, subevent, pair,
 * frame, kinematic, routing, or histogram accumulation work. All fits and
 * delta-t statistics are completed before final presentation normalization is
 * materialized in the returned state.
 */
[[nodiscard]] HistogramAnalysisState analyze_histograms(
    const HBTConfig& config,
    const RawHistogramState& raw
);

/**
 * @brief Require derived state to match the configured raw-histogram layout.
 * @param config Validated HBT configuration defining product/origin/slices.
 * @param raw Complete-sample raw histogram state.
 * @param derived post-sample state to validate against both inputs.
 * @throws std::logic_error If any product, origin, slice, region, selected
 *         count, normalized-bin cardinality, or fit-series cardinality is
 *         structurally inconsistent.
 * @throws std::invalid_argument If configured origin/slice structure is
 *         invalid.
 *
 * The check is read-only and performs no fit, normalization, or serialization.
 * No references or pointers are retained after the synchronous call.
 */
void require_histogram_analysis_layout(
    const HBTConfig& config,
    const RawHistogramState& raw,
    const HistogramAnalysisState& derived
);

}  // namespace hbt

#endif  // HBT_FITS_HISTOGRAM_ANALYSIS_H
