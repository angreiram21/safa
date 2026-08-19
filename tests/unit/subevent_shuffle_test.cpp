/**
 * @file subevent_shuffle_test.cpp
 * @brief Unit tests for deterministic per-subevent particle shuffling.
 */

#include "hbt/event/subevent_shuffle.h"

#include "common/four_vector.h"
#include "hbt/event/origin_flags.h"
#include "hbt/species/species.h"

#include <algorithm>
#include <cstddef>
#include <iostream>
#include <vector>

namespace {

hbt::Particle make_particle(int identity) {
    return {
        hbt::SpeciesId::PiPlus,
        {static_cast<double>(identity), 0.0, 0.0, 0.0},
        {1.0, static_cast<double>(identity), 0.0, 0.0},
        0.139 + static_cast<double>(identity) * 1.0e-6,
        {true, identity % 2 == 0, identity % 3 == 0},
        identity,
        identity % 3 - 1
    };
}

std::vector<hbt::Particle> make_particles() {
    std::vector<hbt::Particle> particles;
    particles.reserve(16U);
    for (int identity = 1; identity <= 16; ++identity) {
        particles.push_back(make_particle(identity));
    }
    return particles;
}

std::vector<int> identities(const std::vector<hbt::Particle>& particles) {
    std::vector<int> result;
    result.reserve(particles.size());
    for (const hbt::Particle& particle : particles) {
        result.push_back(particle.raw_pdg);
    }
    return result;
}

bool particle_fields_match_identity(const hbt::Particle& particle) {
    const int identity = particle.raw_pdg;
    const double value = static_cast<double>(identity);
    return
        particle.species == hbt::SpeciesId::PiPlus &&
        particle.position.x0 == value &&
        particle.momentum.x1 == value &&
        particle.invariant_mass_gev == 0.139 + value * 1.0e-6 &&
        particle.origin.primordial &&
        particle.origin.primordial_rescattering == (identity % 2 == 0) &&
        particle.origin.primordial_rescattering_decay == (identity % 3 == 0) &&
        particle.raw_charge == identity % 3 - 1;
}

bool test_deterministic_particle_permutation() {
    const std::vector<hbt::Particle> original = make_particles();
    std::vector<hbt::Particle> first = original;
    std::vector<hbt::Particle> second = original;

    hbt::shuffle_subevent_particles(first, 3U, 17);
    hbt::shuffle_subevent_particles(second, 3U, 17);

    return
        identities(first) == identities(second) &&
        identities(first) != identities(original);
}

bool test_event_and_subevent_identity_affect_permutation() {
    std::vector<hbt::Particle> reference = make_particles();
    std::vector<hbt::Particle> other_subevent = make_particles();
    std::vector<hbt::Particle> other_event = make_particles();

    hbt::shuffle_subevent_particles(reference, 3U, 17);
    hbt::shuffle_subevent_particles(other_subevent, 3U, 18);
    hbt::shuffle_subevent_particles(other_event, 4U, 17);

    return
        identities(reference) != identities(other_subevent) &&
        identities(reference) != identities(other_event);
}

bool test_particle_identity_is_preserved() {
    const std::vector<hbt::Particle> original = make_particles();
    std::vector<hbt::Particle> shuffled = original;

    hbt::shuffle_subevent_particles(shuffled, 4U, 9);

    std::vector<int> original_ids = identities(original);
    std::vector<int> shuffled_ids = identities(shuffled);
    std::sort(original_ids.begin(), original_ids.end());
    std::sort(shuffled_ids.begin(), shuffled_ids.end());
    if (original_ids != shuffled_ids) {
        return false;
    }

    return std::all_of(
        shuffled.begin(),
        shuffled.end(),
        particle_fields_match_identity
    );
}

bool test_empty_and_single_particle_inputs() {
    std::vector<hbt::Particle> empty;
    hbt::shuffle_subevent_particles(empty, 1U, 0);
    if (!empty.empty()) {
        return false;
    }

    std::vector<hbt::Particle> single{make_particle(7)};
    hbt::shuffle_subevent_particles(single, 1U, 0);
    return
        single.size() == 1U &&
        single.front().raw_pdg == 7 &&
        particle_fields_match_identity(single.front());
}

}  // namespace

/**
 * @brief Run the deterministic per-subevent shuffle unit tests.
 * @return EXIT_SUCCESS when every shuffle contract holds, otherwise
 *         EXIT_FAILURE.
 */
int main() {
    if (!test_deterministic_particle_permutation()) {
        std::cerr << "subevent_shuffle_test: determinism failed.\n";
        return 1;
    }
    if (!test_event_and_subevent_identity_affect_permutation()) {
        std::cerr << "subevent_shuffle_test: subevent independence failed.\n";
        return 1;
    }
    if (!test_particle_identity_is_preserved()) {
        std::cerr << "subevent_shuffle_test: particle preservation failed.\n";
        return 1;
    }
    if (!test_empty_and_single_particle_inputs()) {
        std::cerr << "subevent_shuffle_test: trivial-input handling failed.\n";
        return 1;
    }
    return 0;
}
