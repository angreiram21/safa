/**
 * @file pair_processing_summary.h
 * @brief Data summaries produced by HBT pair processing.
 */

#ifndef SMASH_AFTERBURNER_ANALYSIS_HBT_PAIR_PAIR_PROCESSING_SUMMARY_H
#define SMASH_AFTERBURNER_ANALYSIS_HBT_PAIR_PAIR_PROCESSING_SUMMARY_H

#include "hbt/pair/pair_count_summary.h"
#include "hbt/pair/pair_origin_route_count_summary.h"
#include "hbt/pair/pair_slice_count_summary.h"
#include "hbt/reporting/rejected_pair_report.h"

#include <cstddef>
#include <vector>

namespace hbt {

/**
 * @brief Pair-processing summary for one independent HBT subevent.
 *
 * pair_counts stores every physical pair formed by the iterator. The valid and
 * numerical-rejection summaries partition those formed pairs exactly by
 * primitive channel. A pair is numerically valid only after every numerical
 * gate required for its path succeeds, including frame observables when that
 * gate is active. origin_route_counts records requested routes committed for
 * valid pairs. pair_slice_counts records only valid routed pairs that reach one
 * configured kinetic slice.
 */
struct HBTPairSubeventSummary {
    std::size_t outer_event_number;  ///< One-based outer-event number.
    int subevent_id;                 ///< Subevent identifier from Afterburner.
    PairCountSummary pair_counts;    ///< All physical pairs formed by channel.
    /// Numerically valid pairs by primitive channel.
    PairCountSummary valid_pair_counts;
    /// Numerically rejected formed pairs by channel.
    PairCountSummary numerical_rejection_counts;
    /// Valid pairs routed to requested nested origin selections.
    PairOriginRouteCountSummary origin_route_counts;
    /// Valid routed pairs counted by kinetic slice, origin, and channel.
    PairSliceCountSummary pair_slice_counts;
};

/**
 * @brief Pair-processing summary for one completed HBT run.
 *
 * total_pair_counts preserves the frozen formed-pair counting contract.
 * total_valid_pair_counts and total_numerical_rejection_counts partition those
 * formed pairs exactly by primitive channel. total_origin_route_counts records
 * how numerically valid pairs were routed to requested nested selections.
 * total_pair_slice_counts records the routed subset that reaches configured
 * kinetic slices. numerical_rejections retains every complete recoverable
 * pair-rejection record.
 */
struct HBTPairProcessingSummary {
    PairCountSummary total_pair_counts;  ///< All formed pairs by channel.
    PairCountSummary total_valid_pair_counts;  ///< Numerically valid pairs.
    /// Numerically rejected formed pairs by channel.
    PairCountSummary total_numerical_rejection_counts;
    /// Run-total valid pairs routed by requested origin selection.
    PairOriginRouteCountSummary total_origin_route_counts;
    /// Run-total routed pairs counted by slice, origin, and channel.
    PairSliceCountSummary total_pair_slice_counts;
    RejectedPairReport numerical_rejections;  ///< Complete rejection records.
    std::vector<HBTPairSubeventSummary> subevents;  ///< Local summaries.
};

}  // namespace hbt

#endif  // SMASH_AFTERBURNER_ANALYSIS_HBT_PAIR_PAIR_PROCESSING_SUMMARY_H
