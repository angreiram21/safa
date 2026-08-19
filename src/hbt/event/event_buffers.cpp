/**
 * @file event_buffers.cpp
 * @brief Implementation of per-subevent HBT species buffers.
 */

#include "hbt/event/event_buffers.h"

#include <stdexcept>
#include <utility>

namespace hbt {

std::size_t EventBuffers::index_for_species(SpeciesId species) {
    const auto index = static_cast<std::size_t>(species);

    if (index >= kSpeciesBufferCount) {
        throw std::invalid_argument(
            "EventBuffers: invalid SpeciesId"
        );
    }

    return index;
}

void EventBuffers::add(Particle particle) {
    buffers_[index_for_species(particle.species)].push_back(
        std::move(particle)
    );
}

const std::vector<Particle>& EventBuffers::get(SpeciesId species) const {
    return buffers_[index_for_species(species)];
}

void EventBuffers::clear() noexcept {
    for (std::vector<Particle>& buffer : buffers_) {
        buffer.clear();
    }
}

}  // namespace hbt
