/**
 * @file naming_catalog_test.cpp
 * @brief Cross-cutting contract test for canonical HBT nomenclature.
 *
 * This gate verifies that every canonical species and primitive channel
 * exposes stable ASCII naming derived from the domain catalogue, that physical
 * symbols remain available through species metadata, and that forbidden legacy
 * antiparticle aliases do not re-enter the modular domain.
 */

#include "hbt/channels/channel_catalog.h"
#include "hbt/species/species.h"

#include <array>
#include <cstddef>
#include <cstdlib>
#include <iostream>
#include <string>
#include <string_view>

namespace {

/** Expected canonical naming for one physical species. */
struct SpeciesNamingExpectation {
    hbt::SpeciesId id;  ///< Canonical species identifier.
    std::string_view ascii_token;  ///< Stable canonical ASCII token.
    std::string_view physical_symbol;  ///< Canonical physical symbol.
};

/** Complete canonical nomenclature for every physical SpeciesId. */
constexpr std::array<SpeciesNamingExpectation, 16> kSpeciesNaming{{
    {hbt::SpeciesId::PiPlus, "pi_plus", "π⁺"},
    {hbt::SpeciesId::PiMinus, "pi_minus", "π⁻"},
    {hbt::SpeciesId::PiZero, "pi_zero", "π⁰"},
    {hbt::SpeciesId::KPlus, "k_plus", "K⁺"},
    {hbt::SpeciesId::KMinus, "k_minus", "K⁻"},
    {hbt::SpeciesId::KZero, "k_zero", "K⁰"},
    {hbt::SpeciesId::KZeroBar, "k_zero_bar", "K̄⁰"},
    {hbt::SpeciesId::Proton, "p", "p"},
    {hbt::SpeciesId::ProtonBar, "p_bar", "p̄"},
    {hbt::SpeciesId::Neutron, "n", "n"},
    {hbt::SpeciesId::NeutronBar, "n_bar", "n̄"},
    {hbt::SpeciesId::SigmaPlus, "sigma_plus", "Σ⁺"},
    {hbt::SpeciesId::SigmaBarMinus, "sigma_bar_minus", "Σ̄⁻"},
    {hbt::SpeciesId::SigmaZero, "sigma_zero", "Σ⁰"},
    {hbt::SpeciesId::Lambda, "lambda", "Λ"},
    {hbt::SpeciesId::LambdaBar, "lambda_bar", "Λ̄"}
}};

/** Complete set of primitive channels exposed by PrimitiveChannelId. */
constexpr std::array<hbt::PrimitiveChannelId, 21> kChannelIds{{
    hbt::PrimitiveChannelId::PiPlusPiPlus,
    hbt::PrimitiveChannelId::PiMinusPiMinus,
    hbt::PrimitiveChannelId::PiZeroPiZero,
    hbt::PrimitiveChannelId::KPlusKPlus,
    hbt::PrimitiveChannelId::KMinusKMinus,
    hbt::PrimitiveChannelId::KZeroKZero,
    hbt::PrimitiveChannelId::KZeroBarKZeroBar,
    hbt::PrimitiveChannelId::ProtonProton,
    hbt::PrimitiveChannelId::ProtonBarProtonBar,
    hbt::PrimitiveChannelId::LambdaLambda,
    hbt::PrimitiveChannelId::LambdaBarLambdaBar,
    hbt::PrimitiveChannelId::KMinusProton,
    hbt::PrimitiveChannelId::KPlusProton,
    hbt::PrimitiveChannelId::KZeroBarNeutron,
    hbt::PrimitiveChannelId::KPlusProtonBar,
    hbt::PrimitiveChannelId::KMinusProtonBar,
    hbt::PrimitiveChannelId::PiPlusProton,
    hbt::PrimitiveChannelId::PiMinusProtonBar,
    hbt::PrimitiveChannelId::PiMinusSigmaPlus,
    hbt::PrimitiveChannelId::PiPlusSigmaBarMinus,
    hbt::PrimitiveChannelId::PiZeroSigmaZero
}};

/** Legacy antiparticle aliases forbidden in the modular domain. */
constexpr std::array<std::string_view, 3> kForbiddenAliases{{
    "p_minus",
    "n_minus",
    "lambda_minus"
}};

/**
 * @brief Determine whether a token is a stable lowercase ASCII identifier.
 * @param token Token to inspect.
 * @return true when token is non-empty and uses only [a-z0-9_].
 */
bool is_canonical_ascii_token(std::string_view token) {
    if (token.empty()) {
        return false;
    }

    for (const char character : token) {
        const bool lowercase = character >= 'a' && character <= 'z';
        const bool digit = character >= '0' && character <= '9';
        if (!lowercase && !digit && character != '_') {
            return false;
        }
    }

    return true;
}

/**
 * @brief Detect a forbidden legacy antiparticle alias in a token.
 * @param token Canonical-domain token to inspect.
 * @return true when one forbidden legacy alias occurs in @p token.
 */
bool contains_forbidden_alias(std::string_view token) {
    for (const std::string_view alias : kForbiddenAliases) {
        if (token.find(alias) != std::string_view::npos) {
            return true;
        }
    }

    return false;
}

/**
 * @brief Report a failed nomenclature contract.
 * @param message Short diagnostic describing the failed invariant.
 * @return false for convenient accumulation by callers.
 */
bool fail(std::string_view message) {
    std::cerr << "naming_catalog_test: " << message << ".\n";
    return false;
}

/**
 * @brief Verify canonical metadata for every physical SpeciesId.
 * @return true when all species naming contracts hold.
 */
bool verify_species_naming() {
    if (kSpeciesNaming.size() !=
        static_cast<std::size_t>(hbt::SpeciesId::Count)) {
        return fail("species gate does not cover every SpeciesId");
    }

    for (std::size_t index = 0U; index < kSpeciesNaming.size(); ++index) {
        const SpeciesNamingExpectation& expected = kSpeciesNaming[index];
        const hbt::SpeciesMetadata& metadata =
            hbt::species_metadata(expected.id);

        if (metadata.id != expected.id) {
            return fail("species metadata returned a different SpeciesId");
        }
        if (metadata.ascii_token != expected.ascii_token) {
            return fail("species exposes the wrong canonical ASCII token");
        }
        if (metadata.physical_symbol != expected.physical_symbol) {
            return fail("species exposes the wrong canonical physical symbol");
        }
        if (!is_canonical_ascii_token(metadata.ascii_token)) {
            return fail("species exposes a non-canonical ASCII token");
        }
        if (contains_forbidden_alias(metadata.ascii_token)) {
            return fail("legacy antiparticle alias entered species metadata");
        }

        const auto identified = hbt::identify_species(
            metadata.electric_charge,
            metadata.pdg
        );
        if (!identified.has_value() || *identified != expected.id) {
            return fail("species metadata does not round-trip by charge/PDG");
        }

        for (std::size_t other = index + 1U;
             other < kSpeciesNaming.size();
             ++other) {
            const hbt::SpeciesMetadata& other_metadata =
                hbt::species_metadata(kSpeciesNaming[other].id);
            if (metadata.ascii_token == other_metadata.ascii_token) {
                return fail("species ASCII tokens are not unique");
            }
        }
    }

    return true;
}

/**
 * @brief Build the canonical channel token from endpoint species metadata.
 * @param channel Primitive channel definition to encode.
 * @return Canonical ASCII token implied by the endpoint species.
 */
std::string derived_channel_token(const hbt::PrimitiveChannel& channel) {
    const hbt::SpeciesMetadata& species_a =
        hbt::species_metadata(channel.species_a);
    const hbt::SpeciesMetadata& species_b =
        hbt::species_metadata(channel.species_b);

    std::string token{species_a.ascii_token};
    token += '_';
    token += species_b.ascii_token;
    return token;
}

/**
 * @brief Build the physical channel symbol from endpoint species metadata.
 * @param channel Primitive channel definition to render.
 * @return Physical symbol implied by the canonical endpoint ordering.
 */
std::string derived_channel_symbol(const hbt::PrimitiveChannel& channel) {
    const hbt::SpeciesMetadata& species_a =
        hbt::species_metadata(channel.species_a);
    const hbt::SpeciesMetadata& species_b =
        hbt::species_metadata(channel.species_b);

    std::string symbol{species_a.physical_symbol};
    symbol += species_b.physical_symbol;
    return symbol;
}

/**
 * @brief Verify canonical naming for every PrimitiveChannelId.
 * @return true when channel names derive from canonical species metadata.
 */
bool verify_channel_naming() {
    for (std::size_t index = 0U; index < kChannelIds.size(); ++index) {
        const hbt::PrimitiveChannelId id = kChannelIds[index];
        const hbt::PrimitiveChannel& channel =
            hbt::primitive_channel_definition(id);

        if (channel.id != id) {
            return fail("channel catalogue returned a different channel ID");
        }
        if (!is_canonical_ascii_token(channel.canonical_name)) {
            return fail("channel exposes a non-canonical ASCII token");
        }
        if (contains_forbidden_alias(channel.canonical_name)) {
            return fail("legacy antiparticle alias entered channel naming");
        }
        if (channel.canonical_name != derived_channel_token(channel)) {
            return fail("channel token is not derived from species metadata");
        }
        if (derived_channel_symbol(channel).empty()) {
            return fail("channel has no physical symbol through its species");
        }

        const auto resolved =
            hbt::primitive_channel_id_from_name(channel.canonical_name);
        if (!resolved.has_value() || *resolved != id) {
            return fail("channel name does not round-trip through catalogue");
        }

        for (std::size_t other = index + 1U;
             other < kChannelIds.size();
             ++other) {
            const hbt::PrimitiveChannel& other_channel =
                hbt::primitive_channel_definition(kChannelIds[other]);
            if (channel.canonical_name == other_channel.canonical_name) {
                return fail("primitive-channel ASCII tokens are not unique");
            }
        }
    }

    return true;
}

/**
 * @brief Verify that representative legacy aliases are not catalogue names.
 * @return true when all prohibited aliases are rejected.
 */
bool verify_legacy_alias_rejection() {
    constexpr std::array<std::string_view, 4> aliases{{
        "p_minus_p_minus",
        "n_minus_n_minus",
        "lambda_minus_lambda_minus",
        "pi_minus_p_minus"
    }};

    for (const std::string_view alias : aliases) {
        if (hbt::primitive_channel_id_from_name(alias).has_value()) {
            return fail("legacy primitive-channel alias was accepted");
        }
    }

    return true;
}

}  // namespace

/**
 * @brief Run the cross-cutting canonical-naming contract test.
 * @return EXIT_SUCCESS when all naming invariants hold, otherwise EXIT_FAILURE.
 */
int main() {
    bool success = true;
    success = verify_species_naming() && success;
    success = verify_channel_naming() && success;
    success = verify_legacy_alias_rejection() && success;
    return success ? EXIT_SUCCESS : EXIT_FAILURE;
}
