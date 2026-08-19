/**
 * @file species_requirement_test.cpp
 * @brief Unit tests for HBT particle-species requirement queries.
 *
 * This test verifies that is_species_required() correctly reports whether a
 * particle species occurs in the supplied collection of species required by
 * the current HBT selection.
 */

#include "hbt/selection/species_requirement.h"

#include <cstdlib>
#include <iostream>
#include <vector>

namespace {

    /**
     * @brief Verify one particle-species requirement query.
     *
     * @param test_name Human-readable name of the test case.
     * @param species Particle species to query.
     * @param required_species Collection of species required by HBT.
     * @param expected Expected membership result.
     *
     * @return true when the query result matches @p expected, otherwise false.
     */
    bool verify_query(
        const char* test_name,
        hbt::SpeciesId species,
        const std::vector<hbt::SpeciesId>& required_species,
        bool expected
    ) {
        const bool actual =
            hbt::is_species_required(species, required_species);

        if (actual != expected) {
            std::cerr
                << "species_requirement_test: "
                << test_name
                << " returned "
                << (actual ? "true" : "false")
                << "; expected "
                << (expected ? "true" : "false")
                << ".\n";
            return false;
        }

        return true;
    }

    /**
     * @brief Verify that no species is required by an empty collection.
     *
     * @return true when the query correctly returns false, otherwise false.
     */
    bool verify_empty_collection() {
        const std::vector<hbt::SpeciesId> required_species{};

        return verify_query(
            "empty collection",
            hbt::SpeciesId::PiPlus,
            required_species,
            false
        );
    }

    /**
     * @brief Verify a species present in a single-element collection.
     *
     * @return true when the query correctly returns true, otherwise false.
     */
    bool verify_single_species_present() {
        const std::vector<hbt::SpeciesId> required_species{
            hbt::SpeciesId::Proton
        };

        return verify_query(
            "single species present",
            hbt::SpeciesId::Proton,
            required_species,
            true
        );
    }

    /**
     * @brief Verify a species absent from a single-element collection.
     *
     * @return true when the query correctly returns false, otherwise false.
     */
    bool verify_single_species_absent() {
        const std::vector<hbt::SpeciesId> required_species{
            hbt::SpeciesId::Proton
        };

        return verify_query(
            "single species absent",
            hbt::SpeciesId::PiPlus,
            required_species,
            false
        );
    }

    /**
     * @brief Verify detection of a species at the beginning of a collection.
     *
     * @return true when the query correctly returns true, otherwise false.
     */
    bool verify_species_at_beginning() {
        const std::vector<hbt::SpeciesId> required_species{
            hbt::SpeciesId::KMinus,
            hbt::SpeciesId::Proton,
            hbt::SpeciesId::PiPlus,
            hbt::SpeciesId::Lambda
        };

        return verify_query(
            "species at beginning",
            hbt::SpeciesId::KMinus,
            required_species,
            true
        );
    }

    /**
     * @brief Verify detection of a species in the middle of a collection.
     *
     * @return true when the query correctly returns true, otherwise false.
     */
    bool verify_species_in_middle() {
        const std::vector<hbt::SpeciesId> required_species{
            hbt::SpeciesId::KMinus,
            hbt::SpeciesId::Proton,
            hbt::SpeciesId::PiPlus,
            hbt::SpeciesId::Lambda
        };

        return verify_query(
            "species in middle",
            hbt::SpeciesId::PiPlus,
            required_species,
            true
        );
    }

    /**
     * @brief Verify detection of a species at the end of a collection.
     *
     * @return true when the query correctly returns true, otherwise false.
     */
    bool verify_species_at_end() {
        const std::vector<hbt::SpeciesId> required_species{
            hbt::SpeciesId::KMinus,
            hbt::SpeciesId::Proton,
            hbt::SpeciesId::PiPlus,
            hbt::SpeciesId::Lambda
        };

        return verify_query(
            "species at end",
            hbt::SpeciesId::Lambda,
            required_species,
            true
        );
    }

    /**
     * @brief Verify rejection of a species absent from a multi-species
     * collection.
     *
     * @return true when the query correctly returns false, otherwise false.
     */
    bool verify_species_absent() {
        const std::vector<hbt::SpeciesId> required_species{
            hbt::SpeciesId::KMinus,
            hbt::SpeciesId::Proton,
            hbt::SpeciesId::PiPlus,
            hbt::SpeciesId::Lambda
        };

        return verify_query(
            "species absent",
            hbt::SpeciesId::SigmaZero,
            required_species,
            false
        );
    }

}  // namespace

/**
 * @brief Run the particle-species requirement unit tests.
 *
 * @return EXIT_SUCCESS when every test passes, otherwise EXIT_FAILURE.
 */
int main() {
    bool success = true;

    if (!verify_empty_collection()) {
        success = false;
    }

    if (!verify_single_species_present()) {
        success = false;
    }

    if (!verify_single_species_absent()) {
        success = false;
    }

    if (!verify_species_at_beginning()) {
        success = false;
    }

    if (!verify_species_in_middle()) {
        success = false;
    }

    if (!verify_species_at_end()) {
        success = false;
    }

    if (!verify_species_absent()) {
        success = false;
    }

    return success ? EXIT_SUCCESS : EXIT_FAILURE;
}
