/**
 * @file species_test.cpp
 * @brief Unit tests for HBT particle-species identification and metadata.
 *
 * This test verifies the public contract exposed by species.h:
 * particle identification from electric charge and PDG code, canonical
 * metadata, unsupported-particle handling, and invalid SpeciesId handling.
 */

#include "hbt/species/species.h"

#include <array>
#include <cstddef>
#include <cstdlib>
#include <iostream>
#include <stdexcept>
#include <string_view>

namespace {

    /**
     * @brief Expected public metadata for one supported particle species.
     */
    struct ExpectedSpecies {
        hbt::SpeciesId id;  ///< Expected canonical species identifier.
        int pdg;  ///< Expected PDG particle code.
        int electric_charge;  ///< Expected electric charge in units of |e|.
        std::string_view ascii_token;  ///< Expected canonical ASCII token.
        std::string_view physical_symbol;  ///< Expected physical symbol.
    };

    /**
     * @brief Charge/PDG combination expected to be rejected.
     */
    struct UnsupportedParticle {
        int electric_charge;  ///< Electric charge in units of |e|.
        int pdg;              ///< PDG particle code.
    };

    /**
     * @brief Complete expected catalogue of particle species supported by HBT.
     */
    constexpr std::array<
        ExpectedSpecies,
        static_cast<std::size_t>(hbt::SpeciesId::Count)
    > kExpectedSpecies{{
        {
            hbt::SpeciesId::PiPlus,
            211,
            +1,
            "pi_plus",
            "π⁺"
        },
        {
            hbt::SpeciesId::PiMinus,
            -211,
            -1,
            "pi_minus",
            "π⁻"
        },
        {
            hbt::SpeciesId::PiZero,
            111,
            0,
            "pi_zero",
            "π⁰"
        },

        {
            hbt::SpeciesId::KPlus,
            321,
            +1,
            "k_plus",
            "K⁺"
        },
        {
            hbt::SpeciesId::KMinus,
            -321,
            -1,
            "k_minus",
            "K⁻"
        },
        {
            hbt::SpeciesId::KZero,
            311,
            0,
            "k_zero",
            "K⁰"
        },
        {
            hbt::SpeciesId::KZeroBar,
            -311,
            0,
            "k_zero_bar",
            "K̄⁰"
        },

        {
            hbt::SpeciesId::Proton,
            2212,
            +1,
            "p",
            "p"
        },
        {
            hbt::SpeciesId::ProtonBar,
            -2212,
            -1,
            "p_bar",
            "p̄"
        },

        {
            hbt::SpeciesId::Neutron,
            2112,
            0,
            "n",
            "n"
        },
        {
            hbt::SpeciesId::NeutronBar,
            -2112,
            0,
            "n_bar",
            "n̄"
        },

        {
            hbt::SpeciesId::SigmaPlus,
            3222,
            +1,
            "sigma_plus",
            "Σ⁺"
        },
        {
            hbt::SpeciesId::SigmaBarMinus,
            -3222,
            -1,
            "sigma_bar_minus",
            "Σ̄⁻"
        },
        {
            hbt::SpeciesId::SigmaZero,
            3212,
            0,
            "sigma_zero",
            "Σ⁰"
        },

        {
            hbt::SpeciesId::Lambda,
            3122,
            0,
            "lambda",
            "Λ"
        },
        {
            hbt::SpeciesId::LambdaBar,
            -3122,
            0,
            "lambda_bar",
            "Λ̄"
        }
    }};

    /**
     * @brief Report a failed test condition.
     *
     * @param message Human-readable description of the failed condition.
     *
     * @return Always returns false so the function can be used directly from
     *         verification helpers.
     */
    bool report_failure(std::string_view message) {
        std::cerr << "species_test failure: " << message << '\n';
        return false;
    }

    /**
     * @brief Verify identification and metadata for all supported species.
     *
     * For every expected species, this test checks:
     * - contiguous SpeciesId value,
     * - charge/PDG identification,
     * - SpeciesId,
     * - canonical PDG code,
     * - electric charge,
     * - ASCII token,
     * - physical symbol.
     *
     * @return true if every supported species satisfies the public contract.
     */
    bool verify_supported_species() {
        for (std::size_t index = 0; index < kExpectedSpecies.size(); ++index) {
            const ExpectedSpecies& expected = kExpectedSpecies[index];

            if (static_cast<std::size_t>(expected.id) != index) {
                return report_failure(
                    "SpeciesId values are not contiguous from zero"
                );
            }

            const auto identified = hbt::identify_species(
                expected.electric_charge,
                expected.pdg
            );

            if (!identified.has_value()) {
                return report_failure(
                    "supported charge/PDG combination was not identified"
                );
            }

            if (*identified != expected.id) {
                return report_failure(
                    "charge/PDG combination resolved to the wrong SpeciesId"
                );
            }

            const hbt::SpeciesMetadata& metadata =
                hbt::species_metadata(expected.id);

            if (metadata.id != expected.id) {
                return report_failure(
                    "species metadata contains the wrong SpeciesId"
                );
            }

            if (metadata.pdg != expected.pdg) {
                return report_failure(
                    "species metadata contains the wrong PDG code"
                );
            }

            if (metadata.electric_charge != expected.electric_charge) {
                return report_failure(
                    "species metadata contains the wrong electric charge"
                );
            }

            if (metadata.ascii_token != expected.ascii_token) {
                return report_failure(
                    "species metadata contains the wrong ASCII token"
                );
            }

            if (metadata.physical_symbol != expected.physical_symbol) {
                return report_failure(
                    "species metadata contains the wrong physical symbol"
                );
            }
        }

        return true;
    }

    /**
     * @brief Verify that canonical metadata can identify its own species.
     *
     * This checks the round-trip relation:
     *
     * SpeciesId -> SpeciesMetadata -> charge/PDG -> SpeciesId.
     *
     * @return true if the round-trip succeeds for every supported species.
     */
    bool verify_metadata_round_trip() {
        for (const ExpectedSpecies& expected : kExpectedSpecies) {
            const hbt::SpeciesMetadata& metadata =
                hbt::species_metadata(expected.id);

            const auto identified =
                hbt::identify_species(
                    metadata.electric_charge,
                    metadata.pdg
                );

            if (!identified.has_value()) {
                return report_failure(
                    "metadata charge/PDG pair could not be identified"
                );
            }

            if (*identified != expected.id) {
                return report_failure(
                    "metadata round-trip returned the wrong SpeciesId"
                );
            }
        }

        return true;
    }

    /**
     * @brief Verify that unsupported charge/PDG combinations are rejected.
     *
     * The cases include unsupported PDG codes and deliberately inconsistent
     * charge/PDG combinations for otherwise supported particles.
     *
     * @return true if every unsupported combination returns std::nullopt.
     */
    bool verify_unsupported_species() {
        constexpr std::array<UnsupportedParticle, 8> unsupported{{
            {0, 0},
            {0, 999999},
            {+1, -211},
            {-1, 211},
            {0, 2212},
            {+1, -2212},
            {-1, 3222},
            {+1, -3222}
        }};

        for (const UnsupportedParticle& particle : unsupported) {
            if (hbt::identify_species(
                particle.electric_charge,
                particle.pdg
            ).has_value()) {
                return report_failure(
                    "unsupported charge/PDG combination was accepted"
                );
            }
        }

        return true;
    }

    /**
     * @brief Verify that SpeciesId::Count is not a physical species.
     *
     * @return true if species_metadata() rejects the Count sentinel with
     *         std::invalid_argument.
     */
    bool verify_count_is_not_physical_species() {
        try {
            static_cast<void>(
                hbt::species_metadata(hbt::SpeciesId::Count)
            );
        } catch (const std::invalid_argument&) {
            return true;
        } catch (...) {
            return report_failure(
                "SpeciesId::Count threw an unexpected exception type"
            );
        }

        return report_failure(
            "SpeciesId::Count was accepted as a physical species"
        );
    }

    /**
     * @brief Verify error handling for an invalid SpeciesId value.
     *
     * @return true if species_metadata() throws std::invalid_argument for a
     *         value outside the defined SpeciesId enumeration.
     */
    bool verify_invalid_species_id() {
        const auto invalid_species =
            static_cast<hbt::SpeciesId>(999);

        try {
            static_cast<void>(
                hbt::species_metadata(invalid_species)
            );
        } catch (const std::invalid_argument&) {
            return true;
        } catch (...) {
            return report_failure(
                "invalid SpeciesId threw an unexpected exception type"
            );
        }

        return report_failure(
            "invalid SpeciesId did not throw std::invalid_argument"
        );
    }

}  // namespace

/**
 * @brief Execute all unit tests for the particle-species module.
 *
 * @return EXIT_SUCCESS when every test passes, otherwise EXIT_FAILURE.
 */
int main() {
    if (!verify_supported_species()) {
        return EXIT_FAILURE;
    }

    if (!verify_metadata_round_trip()) {
        return EXIT_FAILURE;
    }

    if (!verify_unsupported_species()) {
        return EXIT_FAILURE;
    }

    if (!verify_count_is_not_physical_species()) {
        return EXIT_FAILURE;
    }

    if (!verify_invalid_species_id()) {
        return EXIT_FAILURE;
    }

    return EXIT_SUCCESS;
}
