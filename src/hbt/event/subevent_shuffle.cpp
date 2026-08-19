/**
 * @file subevent_shuffle.cpp
 * @brief Deterministic accepted-particle shuffling for independent subevents.
 */

#include "hbt/event/subevent_shuffle.h"

#include <algorithm>
#include <cstdint>
#include <random>

namespace hbt {
namespace {

/**
 * @brief Build the deterministic random seed for one independent subevent.
 * @param outer_event_number One-based outer-event number.
 * @param subevent_id Afterburner subevent identifier.
 * @return Stable 64-bit seed depending only on event and subevent identity.
 */
std::uint64_t subevent_shuffle_seed(
    std::size_t outer_event_number,
    int subevent_id
) noexcept {
    constexpr std::uint64_t kSeed = UINT64_C(0x9e3779b97f4a7c15);
    constexpr std::uint64_t kEventMix = UINT64_C(0xbf58476d1ce4e5b9);
    constexpr std::uint64_t kSubeventMix = UINT64_C(0x94d049bb133111eb);

    std::uint64_t value = kSeed;
    value ^= static_cast<std::uint64_t>(outer_event_number) +
             kEventMix + (value << 6U) + (value >> 2U);
    value ^= static_cast<std::uint64_t>(subevent_id) +
             kSubeventMix + (value << 6U) + (value >> 2U);
    return value;
}

}  // namespace

void shuffle_subevent_particles(
    std::vector<Particle>& particles,
    std::size_t outer_event_number,
    int subevent_id
) {
    std::mt19937_64 generator(
        subevent_shuffle_seed(outer_event_number, subevent_id)
    );
    std::shuffle(particles.begin(), particles.end(), generator);
}

}  // namespace hbt
