/**
 * @file required_channels.h
 * @brief Derivation of primitive HBT channels required by a selection.
 *
 * This file declares the operation that derives the unique primitive channels
 * required to produce all final products contained in an HBTSelection.
 *
 * The operation resolves product composition only. It does not validate the
 * selection, construct particle pairs, calculate observables, allocate
 * histograms, or perform any analysis processing.
 */

#ifndef HBT_SELECTION_REQUIRED_CHANNELS_H
#define HBT_SELECTION_REQUIRED_CHANNELS_H

#include "hbt/selection/hbt_selection.h"

#include <vector>

namespace hbt {

    /**
     * @brief Derive the unique primitive channels required by an HBT selection.
     *
     * The products are inspected in their stored order, and the primitive
     * channels inside each product are inspected in their stored order.
     *
     * Each PrimitiveChannelId appears at most once in the returned vector.
     * When a primitive channel occurs in multiple products, its first
     * occurrence determines its position in the result.
     *
     * Preserving first-occurrence order provides deterministic channel ordering
     * without introducing an independent ordering convention.
     *
     * This function does not validate the selection. An empty selection, or a
     * selection containing no primitive channels, therefore produces an empty
     * result.
     *
     * @param selection Requested HBT analysis selection.
     *
     * @return Unique primitive channels required to produce all selected
     *         products, ordered by first occurrence.
     */
    std::vector<PrimitiveChannelId> required_primitive_channels(
        const HBTSelection& selection
    );

}  // namespace hbt

#endif  // HBT_SELECTION_REQUIRED_CHANNELS_H
