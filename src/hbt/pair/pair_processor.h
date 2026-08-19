/**
 * @file pair_processor.h
 * @brief Pair-level kinematics, validation, origin routing, and accounting.
 */

#ifndef SMASH_AFTERBURNER_ANALYSIS_HBT_PAIR_PAIR_PROCESSOR_H
#define SMASH_AFTERBURNER_ANALYSIS_HBT_PAIR_PAIR_PROCESSOR_H

#include "hbt/channels/primitive_channel.h"
#include "hbt/config/hbt_config.h"
#include "hbt/event/event_buffers.h"
#include "hbt/pair/pair_frame_consumer.h"
#include "hbt/pair/pair_processing_summary.h"
#include "hbt/reporting/rejected_pair_report.h"

#include <cstddef>
#include <vector>

namespace hbt {

/**
 * @brief Complete result of processing all required pairs in one subevent.
 */
struct PairSubeventProcessingResult {
    HBTPairSubeventSummary summary;  ///< Counts, routes, and local identity.
    RejectedPairReport numerical_rejections;  ///< Complete local rejections.
};

/**
 * @brief Process all required physical pairs in one current subevent.
 * @param outer_event_number One-based outer-event number for diagnostics.
 * @param subevent_id Current Afterburner subevent identifier.
 * @param buffers Accepted particles from exactly this subevent.
 * @param required_channels Unique required primitive channels in startup order.
 * @param origin_mode Configured nested origin-routing mode.
 * @param pair_slicing Validated kT/mT slicing configuration.
 * @param frame_consumer Immediate downstream frame-observable consumer.
 * @return Local pair counts, origin/slice routes, and rejection records.
 * @throws std::invalid_argument If channels or origin_mode are invalid.
 * @throws std::overflow_error If a local policy counter would overflow.
 * @throws std::logic_error If pair accounting or origin routing is
 *         inconsistent.
 *
 * Every physical pair is formed exactly once by for_each_pair(). PairKinematics
 * is calculated exactly once for that pair. Non-finite kT or mT retains the
 * pair in formed accounting and records one recoverable rejection. Otherwise
 * origin routing is resolved once and kinetic slicing at most once. Frame
 * observables are calculated exactly once when the frame gate is active and
 * are validated immediately. A non-finite frame result is rejected before
 * route counts are committed or any consumer is called. Every remaining pair
 * is counted as numerically valid; in-domain frame results are delivered once
 * by const reference with the already resolved P/PR/PRD routes and slice. No
 * histogramming, fitting, or output is performed here.
 */
[[nodiscard]] PairSubeventProcessingResult process_subevent_pairs(
    std::size_t outer_event_number,
    int subevent_id,
    const EventBuffers& buffers,
    const std::vector<PrimitiveChannelId>& required_channels,
    OriginMode origin_mode,
    const PairSlicingConfig& pair_slicing,
    PairFrameConsumer& frame_consumer
);

/**
 * @brief Accumulate one completed local pair-processing result into run totals.
 * @param total Existing run-total summary to update.
 * @param local Completed local result to consume.
 * @throws std::invalid_argument If channel structures or origin modes differ.
 * @throws std::overflow_error If any run-total count would overflow.
 * @throws std::logic_error If pair accounting or origin routing is
 *         inconsistent.
 *
 * Formed, valid, rejected, origin-route, and slice-route counts are
 * accumulated into temporary copies and checked before run-total counts are
 * committed. Local rejection records are then appended and the local subevent
 * summary is preserved in processing order.
 */
void accumulate_pair_processing_result(
    HBTPairProcessingSummary& total,
    PairSubeventProcessingResult local
);

/**
 * @brief Reduce one ordered event-local pair summary into run totals.
 * @param total Existing run-total pair-processing summary.
 * @param local Completed event-local summary to consume.
 * @throws std::invalid_argument If summary structures or origin modes differ.
 * @throws std::overflow_error If any uint64_t total would overflow.
 * @throws std::logic_error If pair, routing, or rejection invariants fail.
 *
 * Count reductions are integer-only and checked before commit. Rejection
 * records and subevent summaries are appended in the order stored by @p local.
 * Callers therefore obtain deterministic event-major ordering by invoking this
 * function only after all workers join and by reducing event results in
 * ascending order. The function is synchronous and not internally synchronized;
 * @p total requires exclusive mutable ownership for the call.
 */
void accumulate_pair_processing_summary(
    HBTPairProcessingSummary& total,
    HBTPairProcessingSummary local
);

}  // namespace hbt

#endif  // SMASH_AFTERBURNER_ANALYSIS_HBT_PAIR_PAIR_PROCESSOR_H
