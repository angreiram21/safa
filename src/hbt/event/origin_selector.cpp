/**
 * @file origin_selector.cpp
 * @brief Implementation of HBT origin classification and mode eligibility.
 */

#include "hbt/event/origin_selector.h"

#include <stdexcept>

namespace hbt {

OriginFlags classify_origin(
    int pdg_mother1,
    int pdg_mother2
) noexcept {
    const bool first_is_zero = pdg_mother1 == 0;
    const bool second_is_zero = pdg_mother2 == 0;

    const bool primordial = first_is_zero && second_is_zero;
    const bool rescattering =
        !first_is_zero && !second_is_zero;

    return {
        primordial,
        primordial || rescattering,
        true
    };
}

bool is_origin_eligible(
    const OriginFlags& flags,
    OriginMode mode
) {
    switch (mode) {
        case OriginMode::Primordial:
            return flags.primordial;

        case OriginMode::PrimordialRescattering:
            return flags.primordial_rescattering;

        case OriginMode::PrimordialRescatteringDecay:
            return flags.primordial_rescattering_decay;

        case OriginMode::All:
            return flags.primordial_rescattering_decay;
    }

    throw std::invalid_argument(
        "is_origin_eligible(): invalid OriginMode"
    );
}

}  // namespace hbt
