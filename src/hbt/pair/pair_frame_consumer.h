/**
 * @file pair_frame_consumer.h
 * @brief Consumer boundary for already calculated pair-frame observables.
 */

#ifndef SMASH_AFTERBURNER_ANALYSIS_HBT_PAIR_PAIR_FRAME_CONSUMER_H
#define SMASH_AFTERBURNER_ANALYSIS_HBT_PAIR_PAIR_FRAME_CONSUMER_H

#include "hbt/channels/primitive_channel.h"
#include "hbt/pair/pair_frame_observables.h"
#include "hbt/pair/pair_kinematics.h"
#include "hbt/pair/pair_origin_routing.h"
#include "hbt/pair/pair_slice_routing.h"

#include <cstddef>

namespace hbt {

/**
 * @brief Routing context accompanying one frame-observable consumption call.
 *
 * channel_index is produced by primitive-channel traversal and must be reused
 * directly downstream. origin_routes contains all already resolved compatible
 * P/PR/PRD destinations for this one physical pair. pair_slice_route points to
 * the already resolved slice when slicing is active and is nullptr when
 * slicing is disabled. The pointer is valid only for the duration of consume().
 */
struct PairFrameRouteContext {
    std::size_t channel_index;              ///< Reused required-channel index.
    PrimitiveChannelId channel;             ///< Primitive channel identifier.
    PairOriginRoutes origin_routes;          ///< All resolved origin routes.
    const PairSliceRoute* pair_slice_route;  ///< Reused slice, or nullptr.
};

/**
 * @brief Immediate consumer of one already calculated pair-frame result.
 *
 * The processor invokes consume() exactly once for every physical pair that
 * passes the frame-observable gate. Implementations receive existing
 * PairKinematics and PairFrameObservables by const reference together with all
 * already resolved routing destinations. Consumers must consume synchronously
 * and must not retain references or pointers beyond the call.
 */
class PairFrameConsumer {
public:
    /**
     * @brief Destroy the consumer through its interface.
     */
    virtual ~PairFrameConsumer() = default;

    /**
     * @brief Consume one already calculated pair-frame result.
     * @param context Reused channel, origin, and slice-routing context.
     * @param kinematics Existing pair kinematics calculated once upstream.
     * @param observables Existing finite frame observables calculated once.
     *
     * Implementations may accumulate downstream products but must not perform
     * pair formation, pair-kinematics calculation, frame transformations, or
     * route reconstruction.
     */
    virtual void consume(
        const PairFrameRouteContext& context,
        const PairKinematics& kinematics,
        const PairFrameObservables& observables
    ) = 0;
};

}  // namespace hbt

#endif  // SMASH_AFTERBURNER_ANALYSIS_HBT_PAIR_PAIR_FRAME_CONSUMER_H
