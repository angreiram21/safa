/**
 * @file required_species_test.cpp
 * @brief Unit tests for required particle-species derivation.
 *
 * This test verifies that required_species() derives the unique particle
 * species required by primitive HBT channels while preserving canonical
 * species-A/species-B order and first-occurrence order.
 */

#include "hbt/selection/required_species.h"

#include <cstddef>
#include <cstdlib>
#include <iostream>
#include <stdexcept>
#include <vector>

namespace {

    /**
     * @brief Verify the particle species derived for one channel collection.
     *
     * The result must contain exactly the expected species in exactly the
     * expected order.
     *
     * @param test_name Human-readable name of the test case.
     * @param channels Primitive HBT channels to inspect.
     * @param expected Expected unique particle species.
     *
     * @return true when the derived result matches the expectation, otherwise
     *         false.
     */
    bool verify_result(
        const char* test_name,
        const std::vector<hbt::PrimitiveChannelId>& channels,
        const std::vector<hbt::SpeciesId>& expected
    ) {
        const std::vector<hbt::SpeciesId> actual =
            hbt::required_species(channels);

        if (actual.size() != expected.size()) {
            std::cerr
                << "required_species_test: "
                << test_name
                << " returned "
                << actual.size()
                << " species; expected "
                << expected.size()
                << ".\n";
            return false;
        }

        bool success = true;

        for (std::size_t index = 0; index < expected.size(); ++index) {
            if (actual[index] != expected[index]) {
                std::cerr
                    << "required_species_test: "
                    << test_name
                    << " has an incorrect species at index "
                    << index
                    << ".\n";
                success = false;
            }
        }

        return success;
    }

    /**
     * @brief Verify that an empty channel collection requires no species.
     *
     * @return true when the derived result is empty, otherwise false.
     */
    bool verify_empty_channels() {
        const std::vector<hbt::PrimitiveChannelId> channels{};
        const std::vector<hbt::SpeciesId> expected{};

        return verify_result(
            "empty channels",
            channels,
            expected
        );
    }

    /**
     * @brief Verify species derivation for an identical-particle channel.
     *
     * Both canonical channel roles contain the same species, which must appear
     * only once in the result.
     *
     * @return true when the identical channel contributes one species,
     *         otherwise false.
     */
    bool verify_identical_channel() {
        const std::vector<hbt::PrimitiveChannelId> channels{
            hbt::PrimitiveChannelId::PiPlusPiPlus
        };

        const std::vector<hbt::SpeciesId> expected{
            hbt::SpeciesId::PiPlus
        };

        return verify_result(
            "identical channel",
            channels,
            expected
        );
    }

    /**
     * @brief Verify canonical species-A/species-B order for a cross channel.
     *
     * @return true when species A appears before species B, otherwise false.
     */
    bool verify_cross_channel() {
        const std::vector<hbt::PrimitiveChannelId> channels{
            hbt::PrimitiveChannelId::KMinusProton
        };

        const std::vector<hbt::SpeciesId> expected{
            hbt::SpeciesId::KMinus,
            hbt::SpeciesId::Proton
        };

        return verify_result(
            "cross channel",
            channels,
            expected
        );
    }

    /**
     * @brief Verify several channels that do not share particle species.
     *
     * @return true when all species are returned in first-occurrence order,
     *         otherwise false.
     */
    bool verify_channels_without_shared_species() {
        const std::vector<hbt::PrimitiveChannelId> channels{
            hbt::PrimitiveChannelId::PiPlusPiPlus,
            hbt::PrimitiveChannelId::KMinusProton,
            hbt::PrimitiveChannelId::LambdaLambda
        };

        const std::vector<hbt::SpeciesId> expected{
            hbt::SpeciesId::PiPlus,
            hbt::SpeciesId::KMinus,
            hbt::SpeciesId::Proton,
            hbt::SpeciesId::Lambda
        };

        return verify_result(
            "channels without shared species",
            channels,
            expected
        );
    }

    /**
     * @brief Verify deduplication and first-occurrence order for shared
     * species.
     *
     * Species shared by several primitive channels must appear only once.
     * Their position must be determined by the first canonical A/B occurrence
     * encountered while traversing the supplied channels.
     *
     * @return true when shared species are deduplicated without changing order,
     *         otherwise false.
     */
    bool verify_shared_species_and_order() {
        const std::vector<hbt::PrimitiveChannelId> channels{
            hbt::PrimitiveChannelId::KMinusProton,
            hbt::PrimitiveChannelId::PiPlusProton,
            hbt::PrimitiveChannelId::PiPlusSigmaBarMinus,
            hbt::PrimitiveChannelId::PiMinusSigmaPlus
        };

        const std::vector<hbt::SpeciesId> expected{
            hbt::SpeciesId::KMinus,
            hbt::SpeciesId::Proton,
            hbt::SpeciesId::PiPlus,
            hbt::SpeciesId::SigmaBarMinus,
            hbt::SpeciesId::PiMinus,
            hbt::SpeciesId::SigmaPlus
        };

        return verify_result(
            "shared species and order",
            channels,
            expected
        );
    }

    /**
     * @brief Verify propagation of an invalid primitive-channel identifier.
     *
     * @return true when required_species() throws std::invalid_argument,
     *         otherwise false.
     */
    bool verify_invalid_channel_id() {
        const std::vector<hbt::PrimitiveChannelId> channels{
            static_cast<hbt::PrimitiveChannelId>(999)
        };

        try {
            static_cast<void>(
                hbt::required_species(channels)
            );
        } catch (const std::invalid_argument&) {
            return true;
        } catch (...) {
            std::cerr
                << "required_species_test: invalid channel produced an "
                << "unexpected exception type.\n";
            return false;
        }

        std::cerr
            << "required_species_test: invalid channel was not rejected.\n";

        return false;
    }

}  // namespace

/**
 * @brief Run the required particle-species unit tests.
 *
 * @return EXIT_SUCCESS when every test passes, otherwise EXIT_FAILURE.
 */
int main() {
    bool success = true;

    if (!verify_empty_channels()) {
        success = false;
    }

    if (!verify_identical_channel()) {
        success = false;
    }

    if (!verify_cross_channel()) {
        success = false;
    }

    if (!verify_channels_without_shared_species()) {
        success = false;
    }

    if (!verify_shared_species_and_order()) {
        success = false;
    }

    if (!verify_invalid_channel_id()) {
        success = false;
    }

    return success ? EXIT_SUCCESS : EXIT_FAILURE;
}
