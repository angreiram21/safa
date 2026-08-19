/**
 * @file hbt_startup_state.h
 * @brief Resolved HBT state required before event processing begins.
 */

#ifndef HBT_STARTUP_HBT_STARTUP_STATE_H
#define HBT_STARTUP_HBT_STARTUP_STATE_H

#include "hbt/channels/primitive_channel.h"
#include "hbt/selection/hbt_selection.h"
#include "hbt/species/species.h"

#include <vector>

namespace hbt {

/**
 * @brief Resolved HBT selection and its derived processing requirements.
 *
 * This structure contains the scientific HBT selection together with the
 * unique primitive channels and particle species required to execute that
 * selection.
 *
 * The required primitive channels and species are resolved before event
 * processing begins so that later event-preparation code can classify
 * particles without repeatedly deriving these requirements.
 *
 * The order of required_primitive_channels follows stable first occurrence in
 * the stored HBTSelection. The order of required_species follows stable first
 * occurrence while traversing those required primitive channels in canonical
 * species-A then species-B order.
 *
 * This structure does not read events, identify particles, construct pairs,
 * calculate pair observables, or allocate analysis products.
 */
struct HBTStartupState {
    /// Resolved scientific HBT selection.
    HBTSelection selection;
    /// Unique primitive channels required by the resolved selection.
    std::vector<PrimitiveChannelId> required_primitive_channels;
    /// Unique particle species required by the primitive channels.
    std::vector<SpeciesId> required_species;
};

}  // namespace hbt

#endif
