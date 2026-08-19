/**
 * @file required_species.cpp
 * @brief Derivation of particle species required by primitive HBT channels.
 */

#include "hbt/selection/required_species.h"

#include "hbt/channels/channel_catalog.h"

#include <algorithm>

namespace hbt {

    std::vector<SpeciesId> required_species(
        const std::vector<PrimitiveChannelId>& channels
    ) {
        std::vector<SpeciesId> species;

        for (const PrimitiveChannelId channel : channels) {
            const PrimitiveChannel& definition =
                primitive_channel_definition(channel);

            const auto existing_species_a = std::find(
                species.begin(),
                species.end(),
                definition.species_a
            );

            if (existing_species_a == species.end()) {
                species.push_back(definition.species_a);
            }

            const auto existing_species_b = std::find(
                species.begin(),
                species.end(),
                definition.species_b
            );

            if (existing_species_b == species.end()) {
                species.push_back(definition.species_b);
            }
        }

        return species;
    }

}  // namespace hbt
