/**
 * @file species_requirement.cpp
 * @brief Query of particle-species requirements for an HBT selection.
 */

#include "hbt/selection/species_requirement.h"

#include <algorithm>

namespace hbt {

    bool is_species_required(
        SpeciesId species,
        const std::vector<SpeciesId>& required_species
    ) noexcept {
        return std::find(
            required_species.begin(),
            required_species.end(),
            species
        ) != required_species.end();
    }

}  // namespace hbt
