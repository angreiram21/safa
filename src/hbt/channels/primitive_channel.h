/**
 * @file primitive_channel.h
 * @brief Definition of the primitive particle-pair channels used by HBT.
 *
 * This file defines the canonical identities of the primitive HBT channels,
 * their ordered particle-species roles, and their canonical ASCII names.
 *
 * A primitive channel represents one elementary particle-pair classification
 * used by the HBT analysis. Aggregate analysis products may later combine one
 * or more primitive channels, but that composition is outside the
 * responsibility of this file.
 *
 * For cross-species channels, the order of the two species is canonical and
 * must be preserved. The members species_a and species_b therefore represent
 * ordered roles, not an unordered set of particle species.
 */

#ifndef HBT_CHANNELS_PRIMITIVE_CHANNEL_H
#define HBT_CHANNELS_PRIMITIVE_CHANNEL_H

#include "hbt/species/species.h"

#include <string_view>

namespace hbt {

/**
 * @brief Canonical identifier of a primitive HBT particle-pair channel.
 *
 * Each enumerator represents exactly one primitive channel. Identical
 * channels contain two particles of the same species. Cross-species channels
 * encode their canonical species-A/species-B ordering directly in the
 * enumerator name.
 */
enum class PrimitiveChannelId {
    PiPlusPiPlus,          ///< π⁺π⁺.
    PiMinusPiMinus,        ///< π⁻π⁻.
    PiZeroPiZero,          ///< π⁰π⁰.

    KPlusKPlus,            ///< K⁺K⁺.
    KMinusKMinus,          ///< K⁻K⁻.
    KZeroKZero,            ///< K⁰K⁰.
    KZeroBarKZeroBar,      ///< K̄⁰K̄⁰.

    ProtonProton,          ///< pp.
    ProtonBarProtonBar,    ///< p̄p̄.

    LambdaLambda,          ///< ΛΛ.
    LambdaBarLambdaBar,    ///< Λ̄Λ̄.

    KMinusProton,          ///< K⁻p.
    KPlusProton,           ///< K⁺p.
    KZeroBarNeutron,       ///< K̄⁰n.
    KPlusProtonBar,        ///< K⁺p̄.
    KMinusProtonBar,       ///< K⁻p̄.

    PiPlusProton,          ///< π⁺p.
    PiMinusProtonBar,      ///< π⁻p̄.
    PiMinusSigmaPlus,      ///< π⁻Σ⁺.
    PiPlusSigmaBarMinus,   ///< π⁺Σ̄⁻.
    PiZeroSigmaZero        ///< π⁰Σ⁰.
};

/**
 * @brief Canonical definition of one primitive HBT channel.
 *
 * The channel identifier, canonical ordered species roles, and canonical ASCII
 * name form one immutable conceptual relationship supplied by the channel
 * catalogue.
 *
 * For identical channels, species_a and species_b are equal. For cross-species
 * channels, their order must not be changed merely because a later calculation
 * is symmetric under particle exchange.
 */
struct PrimitiveChannel {
    PrimitiveChannelId id;          ///< Canonical primitive-channel identifier.
    SpeciesId species_a;            ///< Canonical species occupying role A.
    SpeciesId species_b;            ///< Canonical species occupying role B.
    /// Stable ASCII configuration/output name.
    std::string_view canonical_name;
};

}  // namespace hbt

#endif  // HBT_CHANNELS_PRIMITIVE_CHANNEL_H
