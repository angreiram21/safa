/**
 * @file species.h
 * @brief Definition of the particle species used by the HBT analysis.
 *
 * This file defines the canonical particle identities used throughout the
 * modular HBT implementation.
 *
 * Particle identity is represented by SpeciesId. Human-readable names,
 * physical symbols, PDG codes, electric charges, and ASCII configuration
 * tokens are exposed through SpeciesMetadata.
 *
 * Species names follow the physical identity implied by the PDG code.
 * Alternative textual labels are not a source of truth for nomenclature.
 */

#ifndef HBT_DOMAIN_SPECIES_H
#define HBT_DOMAIN_SPECIES_H

#include <optional>
#include <string_view>

namespace hbt {

    /**
     * @brief Canonical particle species supported by the HBT analysis.
     *
     * Particle/antiparticle identity is represented explicitly.
     *
     * Names such as ProtonBar, NeutronBar, SigmaBarMinus, and LambdaBar
     * describe physical antiparticle identity and must not be replaced by
     * ambiguous charge-based aliases such as ProtonMinus or LambdaMinus.
     *
     * Physical species values are contiguous and start at zero. Count is a
     * sentinel equal to the number of physical species and is not itself a
     * particle species.
     *
     * The canonical identity of each species is determined by its PDG code.
     */
    enum class SpeciesId {
        PiPlus,        ///< Positively charged pion, π⁺.
        PiMinus,       ///< Negatively charged pion, π⁻.
        PiZero,        ///< Neutral pion, π⁰.

        KPlus,         ///< Positively charged kaon, K⁺.
        KMinus,        ///< Negatively charged kaon, K⁻.
        KZero,         ///< Neutral kaon, K⁰.
        KZeroBar,      ///< Neutral antikaon, K̄⁰.

        Proton,        ///< Proton, p.
        ProtonBar,     ///< Antiproton, p̄.

        Neutron,       ///< Neutron, n.
        NeutronBar,    ///< Antineutron, n̄.

        SigmaPlus,     ///< Sigma-plus baryon, Σ⁺.
        SigmaBarMinus, ///< Anti-Sigma-plus baryon, Σ̄⁻.
        SigmaZero,     ///< Neutral sigma baryon, Σ⁰.

        Lambda,        ///< Lambda baryon, Λ.
        LambdaBar,     ///< Antilambda baryon, Λ̄.

        Count          ///< Number of physical species; not a species.
    };

    /**
     * @brief Immutable canonical metadata associated with one particle species.
     *
     * All user-visible and machine-readable species names must originate from
     * this metadata rather than from independent lookup tables in parsers,
     * writers, or analysis modules.
     *
     * This structure contains only canonical information. Alternative textual
     * aliases must not be stored here.
     */
    struct SpeciesMetadata {
        SpeciesId id;                  ///< Canonical species identifier.
        int pdg;                       ///< Canonical PDG particle code.
        int electric_charge;           ///< Electric charge in units of |e|.
        std::string_view ascii_token;  ///< Canonical ASCII token.
        std::string_view physical_symbol;  ///< Physical symbol, e.g. "p̄".
    };

    /**
     * @brief Identify a supported particle from electric charge and PDG code.
     *
     * Particle identification uses both electric charge and PDG code.
     *
     * Once the charge/PDG pair has been recognized, the corresponding canonical
     * SpeciesId determines the physical identity used throughout the HBT
     * analysis.
     *
     * A future simplification to PDG-only identification would require an
     * explicit equivalence test and must not be introduced silently.
     *
     * @param electric_charge Electric charge in units of |e|.
     * @param pdg PDG particle code.
     *
     * @return The corresponding SpeciesId when the combination is supported,
     *         or std::nullopt when the charge/PDG combination is not part of
     *         the HBT species catalogue.
     */
    std::optional<SpeciesId> identify_species(
        int electric_charge,
        int pdg
    ) noexcept;

    /**
     * @brief Obtain the canonical metadata for a particle species.
     *
     * This function provides the single source of truth for the PDG code,
     * electric charge, ASCII token, and physical symbol associated with a
     * SpeciesId.
     *
     * Parsers, writers, channel catalogues, and analysis modules must obtain
     * species metadata through this function rather than maintaining
     * independent lookup tables.
     *
     * @param species Canonical particle species.
     *
     * @return Immutable metadata associated with the requested species.
     *
     * @throws std::invalid_argument If species is not a physical SpeciesId.
     */
    const SpeciesMetadata& species_metadata(SpeciesId species);

}  // namespace hbt

#endif  // HBT_DOMAIN_SPECIES_H
