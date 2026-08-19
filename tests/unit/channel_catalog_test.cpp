/**
 * @file channel_catalog_test.cpp
 * @brief Unit tests for the canonical primitive HBT channel catalogue.
 */

#include "hbt/channels/channel_catalog.h"

#include <array>
#include <cstddef>
#include <cstdlib>
#include <iostream>
#include <stdexcept>
#include <string_view>

namespace {

/**
 * @brief Expected complete metadata for one primitive HBT channel.
 */
struct ExpectedChannel {
    hbt::PrimitiveChannelId id;  ///< Expected primitive-channel identifier.
    hbt::SpeciesId species_a;    ///< Expected canonical species in role A.
    hbt::SpeciesId species_b;    ///< Expected canonical species in role B.
    std::string_view name;       ///< Expected stable canonical ASCII name.
};

/**
 * @brief Complete expected primitive-channel catalogue.
 */
constexpr std::array<ExpectedChannel, 21> kExpectedChannels{{
    {
        hbt::PrimitiveChannelId::PiPlusPiPlus,
        hbt::SpeciesId::PiPlus,
        hbt::SpeciesId::PiPlus,
        "pi_plus_pi_plus"
    },
    {
        hbt::PrimitiveChannelId::PiMinusPiMinus,
        hbt::SpeciesId::PiMinus,
        hbt::SpeciesId::PiMinus,
        "pi_minus_pi_minus"
    },
    {
        hbt::PrimitiveChannelId::PiZeroPiZero,
        hbt::SpeciesId::PiZero,
        hbt::SpeciesId::PiZero,
        "pi_zero_pi_zero"
    },
    {
        hbt::PrimitiveChannelId::KPlusKPlus,
        hbt::SpeciesId::KPlus,
        hbt::SpeciesId::KPlus,
        "k_plus_k_plus"
    },
    {
        hbt::PrimitiveChannelId::KMinusKMinus,
        hbt::SpeciesId::KMinus,
        hbt::SpeciesId::KMinus,
        "k_minus_k_minus"
    },
    {
        hbt::PrimitiveChannelId::KZeroKZero,
        hbt::SpeciesId::KZero,
        hbt::SpeciesId::KZero,
        "k_zero_k_zero"
    },
    {
        hbt::PrimitiveChannelId::KZeroBarKZeroBar,
        hbt::SpeciesId::KZeroBar,
        hbt::SpeciesId::KZeroBar,
        "k_zero_bar_k_zero_bar"
    },
    {
        hbt::PrimitiveChannelId::ProtonProton,
        hbt::SpeciesId::Proton,
        hbt::SpeciesId::Proton,
        "p_p"
    },
    {
        hbt::PrimitiveChannelId::ProtonBarProtonBar,
        hbt::SpeciesId::ProtonBar,
        hbt::SpeciesId::ProtonBar,
        "p_bar_p_bar"
    },
    {
        hbt::PrimitiveChannelId::LambdaLambda,
        hbt::SpeciesId::Lambda,
        hbt::SpeciesId::Lambda,
        "lambda_lambda"
    },
    {
        hbt::PrimitiveChannelId::LambdaBarLambdaBar,
        hbt::SpeciesId::LambdaBar,
        hbt::SpeciesId::LambdaBar,
        "lambda_bar_lambda_bar"
    },
    {
        hbt::PrimitiveChannelId::KMinusProton,
        hbt::SpeciesId::KMinus,
        hbt::SpeciesId::Proton,
        "k_minus_p"
    },
    {
        hbt::PrimitiveChannelId::KPlusProton,
        hbt::SpeciesId::KPlus,
        hbt::SpeciesId::Proton,
        "k_plus_p"
    },
    {
        hbt::PrimitiveChannelId::KZeroBarNeutron,
        hbt::SpeciesId::KZeroBar,
        hbt::SpeciesId::Neutron,
        "k_zero_bar_n"
    },
    {
        hbt::PrimitiveChannelId::KPlusProtonBar,
        hbt::SpeciesId::KPlus,
        hbt::SpeciesId::ProtonBar,
        "k_plus_p_bar"
    },
    {
        hbt::PrimitiveChannelId::KMinusProtonBar,
        hbt::SpeciesId::KMinus,
        hbt::SpeciesId::ProtonBar,
        "k_minus_p_bar"
    },
    {
        hbt::PrimitiveChannelId::PiPlusProton,
        hbt::SpeciesId::PiPlus,
        hbt::SpeciesId::Proton,
        "pi_plus_p"
    },
    {
        hbt::PrimitiveChannelId::PiMinusProtonBar,
        hbt::SpeciesId::PiMinus,
        hbt::SpeciesId::ProtonBar,
        "pi_minus_p_bar"
    },
    {
        hbt::PrimitiveChannelId::PiMinusSigmaPlus,
        hbt::SpeciesId::PiMinus,
        hbt::SpeciesId::SigmaPlus,
        "pi_minus_sigma_plus"
    },
    {
        hbt::PrimitiveChannelId::PiPlusSigmaBarMinus,
        hbt::SpeciesId::PiPlus,
        hbt::SpeciesId::SigmaBarMinus,
        "pi_plus_sigma_bar_minus"
    },
    {
        hbt::PrimitiveChannelId::PiZeroSigmaZero,
        hbt::SpeciesId::PiZero,
        hbt::SpeciesId::SigmaZero,
        "pi_zero_sigma_zero"
    }
}};

/**
 * @brief Report one catalogue mismatch.
 * @param index Index of the channel whose field differed.
 * @param field Name of the mismatching field.
 * @return Always false.
 */
bool fail_channel(std::size_t index, const char* field) {
    std::cerr
        << "channel_catalog_test: channel " << index
        << " has incorrect " << field << ".\n";
    return false;
}

/**
 * @brief Verify all canonical definitions and reverse name lookup.
 * @return true when all 21 definitions and lookups match exactly.
 */
bool verify_supported_channels() {
    bool success = true;

    for (std::size_t index = 0U; index < kExpectedChannels.size(); ++index) {
        const ExpectedChannel& expected = kExpectedChannels[index];
        const hbt::PrimitiveChannel& actual =
            hbt::primitive_channel_definition(expected.id);

        if (actual.id != expected.id) {
            success = fail_channel(index, "channel identifier") && success;
        }
        if (actual.species_a != expected.species_a) {
            success = fail_channel(index, "species_a") && success;
        }
        if (actual.species_b != expected.species_b) {
            success = fail_channel(index, "species_b") && success;
        }
        if (actual.canonical_name != expected.name) {
            success = fail_channel(index, "canonical_name") && success;
        }

        const auto resolved =
            hbt::primitive_channel_id_from_name(expected.name);
        if (!resolved.has_value() || resolved.value() != expected.id) {
            success = fail_channel(index, "name lookup") && success;
        }
    }

    return success;
}

/**
 * @brief Verify all canonical channel names are unique.
 * @return true when no two catalogue entries expose the same name.
 */
bool verify_unique_names() {
    for (std::size_t left = 0U; left < kExpectedChannels.size(); ++left) {
        const std::string_view left_name =
            hbt::primitive_channel_definition(
                kExpectedChannels[left].id
            ).canonical_name;

        for (std::size_t right = left + 1U;
             right < kExpectedChannels.size();
             ++right) {
            const std::string_view right_name =
                hbt::primitive_channel_definition(
                    kExpectedChannels[right].id
                ).canonical_name;

            if (left_name == right_name) {
                std::cerr
                    << "channel_catalog_test: duplicate canonical name.\n";
                return false;
            }
        }
    }

    return true;
}

/**
 * @brief Verify invalid identifiers and unknown names are rejected.
 * @return true when both invalid queries follow their documented contracts.
 */
bool verify_invalid_queries() {
    try {
        static_cast<void>(hbt::primitive_channel_definition(
            static_cast<hbt::PrimitiveChannelId>(999)
        ));
    } catch (const std::invalid_argument&) {
        return !hbt::primitive_channel_id_from_name(
            "not_a_primitive_channel"
        ).has_value();
    } catch (...) {
        std::cerr
            << "channel_catalog_test: invalid ID used wrong exception type.\n";
        return false;
    }

    std::cerr << "channel_catalog_test: invalid ID was not rejected.\n";
    return false;
}

}  // namespace

/**
 * @brief Run the primitive-channel catalogue unit tests.
 * @return EXIT_SUCCESS when every test passes, otherwise EXIT_FAILURE.
 */
int main() {
    bool success = true;
    success = verify_supported_channels() && success;
    success = verify_unique_names() && success;
    success = verify_invalid_queries() && success;
    return success ? EXIT_SUCCESS : EXIT_FAILURE;
}
