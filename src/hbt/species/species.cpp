/**
 * @file species.cpp
 * @brief Canonical metadata and identification rules for HBT particle species.
 *
 * This file implements the canonical particle-species metadata used by the
 * HBT analysis.
 *
 * Particle identification is defined by electric charge and PDG code.
 * Species names, ASCII tokens, and physical symbols follow the physical
 * identity represented by the PDG code.
 */

#include "species.h"

#include <array>
#include <cstddef>
#include <stdexcept>

namespace hbt {

    namespace {

        /// Number of physical species represented by SpeciesId.
        constexpr std::size_t kSpeciesCount =
            static_cast<std::size_t>(SpeciesId::Count);

        /**
         * @brief Canonical metadata table for all HBT particle species.
         *
         * This table is the single source of truth for the relationship
         * between:
         * - SpeciesId,
         * - PDG code,
         * - electric charge,
         * - canonical ASCII token,
         * - physical symbol.
         *
         * All particle names and symbols used by the HBT analysis originate
         * from this table.
         */
        constexpr std::array<SpeciesMetadata, kSpeciesCount>
            kSpeciesMetadata{{
                {
                    SpeciesId::PiPlus,
                    211,
                    +1,
                    "pi_plus",
                    "π⁺"
                },
                {
                    SpeciesId::PiMinus,
                    -211,
                    -1,
                    "pi_minus",
                    "π⁻"
                },
                {
                    SpeciesId::PiZero,
                    111,
                    0,
                    "pi_zero",
                    "π⁰"
                },

                {
                    SpeciesId::KPlus,
                    321,
                    +1,
                    "k_plus",
                    "K⁺"
                },
                {
                    SpeciesId::KMinus,
                    -321,
                    -1,
                    "k_minus",
                    "K⁻"
                },
                {
                    SpeciesId::KZero,
                    311,
                    0,
                    "k_zero",
                    "K⁰"
                },
                {
                    SpeciesId::KZeroBar,
                    -311,
                    0,
                    "k_zero_bar",
                    "K̄⁰"
                },

                {
                    SpeciesId::Proton,
                    2212,
                    +1,
                    "p",
                    "p"
                },
                {
                    SpeciesId::ProtonBar,
                    -2212,
                    -1,
                    "p_bar",
                    "p̄"
                },

                {
                    SpeciesId::Neutron,
                    2112,
                    0,
                    "n",
                    "n"
                },
                {
                    SpeciesId::NeutronBar,
                    -2112,
                    0,
                    "n_bar",
                    "n̄"
                },

                {
                    SpeciesId::SigmaPlus,
                    3222,
                    +1,
                    "sigma_plus",
                    "Σ⁺"
                },
                {
                    SpeciesId::SigmaBarMinus,
                    -3222,
                    -1,
                    "sigma_bar_minus",
                    "Σ̄⁻"
                },
                {
                    SpeciesId::SigmaZero,
                    3212,
                    0,
                    "sigma_zero",
                    "Σ⁰"
                },

                {
                    SpeciesId::Lambda,
                    3122,
                    0,
                    "lambda",
                    "Λ"
                },
                {
                    SpeciesId::LambdaBar,
                    -3122,
                    0,
                    "lambda_bar",
                    "Λ̄"
                }
            }};

        /**
         * @brief Verify metadata order against contiguous SpeciesId values.
         *
         * @return true when every physical SpeciesId maps to the metadata entry
         *         at the same zero-based index.
         */
        constexpr bool species_metadata_matches_enum_order() noexcept {
            for (std::size_t index = 0; index < kSpeciesCount; ++index) {
                if (static_cast<std::size_t>(kSpeciesMetadata[index].id) !=
                    index) {
                    return false;
                }
            }

            return true;
        }

        static_assert(
            species_metadata_matches_enum_order(),
            "SpeciesId values and species metadata must remain contiguous"
        );

    }  // namespace

    std::optional<SpeciesId> identify_species(
        int electric_charge,
        int pdg
    ) noexcept {
        for (const SpeciesMetadata& metadata : kSpeciesMetadata) {
            if (metadata.electric_charge == electric_charge &&
                metadata.pdg == pdg) {
                return metadata.id;
            }
        }

        return std::nullopt;
    }

    const SpeciesMetadata& species_metadata(SpeciesId species) {
        for (const SpeciesMetadata& metadata : kSpeciesMetadata) {
            if (metadata.id == species) {
                return metadata;
            }
        }

        throw std::invalid_argument(
            "species_metadata(): invalid SpeciesId"
        );
    }

}  // namespace hbt
