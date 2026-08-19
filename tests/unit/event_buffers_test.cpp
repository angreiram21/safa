/**
 * @file event_buffers_test.cpp
 * @brief Unit tests for per-subevent HBT particle buffers.
 */

#include "hbt/event/event_buffers.h"

#include <cstdlib>
#include <exception>
#include <iostream>
#include <stdexcept>
#include <utility>

namespace {

/**
 * @brief Construct one accepted particle with recognizable test values.
 * @param species Canonical species stored in the particle.
 * @param marker Numeric marker copied into position and momentum.
 * @return Complete Particle value for buffer tests.
 */
hbt::Particle make_particle(hbt::SpeciesId species, double marker) {
    return {
        species,
        {marker, marker + 1.0, marker + 2.0, marker + 3.0},
        {marker + 4.0, marker + 5.0, marker + 6.0, marker + 7.0},
        marker + 8.0,
        {true, true, true},
        211,
        1
    };
}

/**
 * @brief Verify species separation and insertion-order preservation.
 * @return `true` when independent species buffers contain the expected values.
 */
bool verify_species_grouping() {
    hbt::EventBuffers buffers;

    buffers.add(make_particle(hbt::SpeciesId::PiPlus, 1.0));
    buffers.add(make_particle(hbt::SpeciesId::KPlus, 10.0));
    buffers.add(make_particle(hbt::SpeciesId::PiPlus, 20.0));

    const auto& pions = buffers.get(hbt::SpeciesId::PiPlus);
    const auto& kaons = buffers.get(hbt::SpeciesId::KPlus);
    const auto& protons = buffers.get(hbt::SpeciesId::Proton);

    if (pions.size() != 2U || kaons.size() != 1U || !protons.empty()) {
        std::cerr << "event_buffers_test: species grouping is incorrect.\n";
        return false;
    }

    if (pions[0].position.x0 != 1.0 || pions[1].position.x0 != 20.0 ||
        kaons[0].position.x0 != 10.0) {
        std::cerr << "event_buffers_test: insertion order was not preserved.\n";
        return false;
    }

    return true;
}

/**
 * @brief Verify clear removes every current-subevent particle.
 * @return `true` when all previously populated buffers become empty.
 */
bool verify_clear() {
    hbt::EventBuffers buffers;
    buffers.add(make_particle(hbt::SpeciesId::PiMinus, 2.0));
    buffers.add(make_particle(hbt::SpeciesId::LambdaBar, 3.0));

    buffers.clear();

    return buffers.get(hbt::SpeciesId::PiMinus).empty() &&
           buffers.get(hbt::SpeciesId::LambdaBar).empty();
}

/**
 * @brief Verify that the SpeciesId::Count sentinel is rejected.
 * @return `true` when add and get both throw std::invalid_argument for Count.
 */
bool verify_count_sentinel_rejected() {
    hbt::EventBuffers buffers;
    bool add_threw = false;
    bool get_threw = false;

    try {
        buffers.add(make_particle(hbt::SpeciesId::Count, 1.0));
    } catch (const std::invalid_argument&) {
        add_threw = true;
    }

    try {
        static_cast<void>(buffers.get(hbt::SpeciesId::Count));
    } catch (const std::invalid_argument&) {
        get_threw = true;
    }

    if (!add_threw || !get_threw) {
        std::cerr
            << "event_buffers_test: SpeciesId::Count was not rejected.\n";
        return false;
    }

    return true;
}

/**
 * @brief Verify invalid SpeciesId values are rejected for add and get.
 * @return `true` when both invalid access paths throw std::invalid_argument.
 */
bool verify_invalid_species_rejected() {
    const auto invalid = static_cast<hbt::SpeciesId>(999);
    hbt::EventBuffers buffers;
    bool add_threw = false;
    bool get_threw = false;

    try {
        buffers.add(make_particle(invalid, 1.0));
    } catch (const std::invalid_argument&) {
        add_threw = true;
    }

    try {
        static_cast<void>(buffers.get(invalid));
    } catch (const std::invalid_argument&) {
        get_threw = true;
    }

    if (!add_threw || !get_threw) {
        std::cerr
            << "event_buffers_test: invalid SpeciesId was not rejected.\n";
        return false;
    }

    return true;
}

}  // namespace

/**
 * @brief Run the complete EventBuffers unit-test collection.
 * @return EXIT_SUCCESS when every test passes, otherwise EXIT_FAILURE.
 */
int main() {
    bool success = true;
    success = verify_species_grouping() && success;
    success = verify_clear() && success;
    success = verify_count_sentinel_rejected() && success;
    success = verify_invalid_species_rejected() && success;
    return success ? EXIT_SUCCESS : EXIT_FAILURE;
}
