/**
 * @file species_requirement.h
 * @brief Query of particle-species requirements for an HBT selection.
 *
 * This file declares the operation used to determine whether a particle
 * species is required by the currently selected primitive HBT channels.
 *
 * The required-species collection must have been derived before event
 * processing begins. This operation performs membership testing only; particle
 * identification, event reading, particle storage, pair construction, and
 * observable calculation are outside its responsibility.
 */

#ifndef HBT_SELECTION_SPECIES_REQUIREMENT_H
#define HBT_SELECTION_SPECIES_REQUIREMENT_H

#include "hbt/species/species.h"

#include <vector>

namespace hbt {

    /**
     * @brief Determine whether a particle species is required by HBT.
     *
     * The function tests whether @p species occurs in the supplied collection
     * of species required by the selected primitive HBT channels.
     *
     * The supplied collection is not modified.
     *
     * @param species Canonical particle species to test.
     * @param required_species Species required by the current HBT selection.
     *
     * @return true if @p species is required by HBT, otherwise false.
     */
    bool is_species_required(
        SpeciesId species,
        const std::vector<SpeciesId>& required_species
    ) noexcept;

}  // namespace hbt

#endif  // HBT_SELECTION_SPECIES_REQUIREMENT_H
