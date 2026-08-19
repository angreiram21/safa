/**
 * @file channel_catalog.cpp
 * @brief Canonical definitions of the primitive HBT channels.
 */

#include "hbt/channels/channel_catalog.h"

#include <stdexcept>

namespace hbt {
namespace {

/**
 * @brief Canonical definitions of all primitive HBT channels.
 *
 * This table is the single production source of truth for channel identity,
 * canonical A/B species roles, and stable ASCII channel names.
 */
constexpr PrimitiveChannel kPrimitiveChannels[] = {
    {
        PrimitiveChannelId::PiPlusPiPlus,
        SpeciesId::PiPlus,
        SpeciesId::PiPlus,
        "pi_plus_pi_plus"
    },
    {
        PrimitiveChannelId::PiMinusPiMinus,
        SpeciesId::PiMinus,
        SpeciesId::PiMinus,
        "pi_minus_pi_minus"
    },
    {
        PrimitiveChannelId::PiZeroPiZero,
        SpeciesId::PiZero,
        SpeciesId::PiZero,
        "pi_zero_pi_zero"
    },
    {
        PrimitiveChannelId::KPlusKPlus,
        SpeciesId::KPlus,
        SpeciesId::KPlus,
        "k_plus_k_plus"
    },
    {
        PrimitiveChannelId::KMinusKMinus,
        SpeciesId::KMinus,
        SpeciesId::KMinus,
        "k_minus_k_minus"
    },
    {
        PrimitiveChannelId::KZeroKZero,
        SpeciesId::KZero,
        SpeciesId::KZero,
        "k_zero_k_zero"
    },
    {
        PrimitiveChannelId::KZeroBarKZeroBar,
        SpeciesId::KZeroBar,
        SpeciesId::KZeroBar,
        "k_zero_bar_k_zero_bar"
    },
    {
        PrimitiveChannelId::ProtonProton,
        SpeciesId::Proton,
        SpeciesId::Proton,
        "p_p"
    },
    {
        PrimitiveChannelId::ProtonBarProtonBar,
        SpeciesId::ProtonBar,
        SpeciesId::ProtonBar,
        "p_bar_p_bar"
    },
    {
        PrimitiveChannelId::LambdaLambda,
        SpeciesId::Lambda,
        SpeciesId::Lambda,
        "lambda_lambda"
    },
    {
        PrimitiveChannelId::LambdaBarLambdaBar,
        SpeciesId::LambdaBar,
        SpeciesId::LambdaBar,
        "lambda_bar_lambda_bar"
    },
    {
        PrimitiveChannelId::KMinusProton,
        SpeciesId::KMinus,
        SpeciesId::Proton,
        "k_minus_p"
    },
    {
        PrimitiveChannelId::KPlusProton,
        SpeciesId::KPlus,
        SpeciesId::Proton,
        "k_plus_p"
    },
    {
        PrimitiveChannelId::KZeroBarNeutron,
        SpeciesId::KZeroBar,
        SpeciesId::Neutron,
        "k_zero_bar_n"
    },
    {
        PrimitiveChannelId::KPlusProtonBar,
        SpeciesId::KPlus,
        SpeciesId::ProtonBar,
        "k_plus_p_bar"
    },
    {
        PrimitiveChannelId::KMinusProtonBar,
        SpeciesId::KMinus,
        SpeciesId::ProtonBar,
        "k_minus_p_bar"
    },
    {
        PrimitiveChannelId::PiPlusProton,
        SpeciesId::PiPlus,
        SpeciesId::Proton,
        "pi_plus_p"
    },
    {
        PrimitiveChannelId::PiMinusProtonBar,
        SpeciesId::PiMinus,
        SpeciesId::ProtonBar,
        "pi_minus_p_bar"
    },
    {
        PrimitiveChannelId::PiMinusSigmaPlus,
        SpeciesId::PiMinus,
        SpeciesId::SigmaPlus,
        "pi_minus_sigma_plus"
    },
    {
        PrimitiveChannelId::PiPlusSigmaBarMinus,
        SpeciesId::PiPlus,
        SpeciesId::SigmaBarMinus,
        "pi_plus_sigma_bar_minus"
    },
    {
        PrimitiveChannelId::PiZeroSigmaZero,
        SpeciesId::PiZero,
        SpeciesId::SigmaZero,
        "pi_zero_sigma_zero"
    }
};

}  // namespace

const PrimitiveChannel& primitive_channel_definition(
    PrimitiveChannelId channel
) {
    for (const PrimitiveChannel& definition : kPrimitiveChannels) {
        if (definition.id == channel) {
            return definition;
        }
    }

    throw std::invalid_argument(
        "primitive_channel_definition(): invalid PrimitiveChannelId"
    );
}

std::optional<PrimitiveChannelId> primitive_channel_id_from_name(
    std::string_view canonical_name
) noexcept {
    for (const PrimitiveChannel& definition : kPrimitiveChannels) {
        if (definition.canonical_name == canonical_name) {
            return definition.id;
        }
    }

    return std::nullopt;
}

}  // namespace hbt
