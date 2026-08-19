/**
 * @file origin_selector.h
 * @brief Nested HBT origin classification and mode eligibility.
 *
 * This file declares the event-preparation operations that translate the raw
 * Afterburner mother-PDG fields into nested HBT origin flags and test those
 * flags against the configured OriginMode.
 *
 * Kinematic acceptance occurs before these operations. Sampler lookup,
 * emission-point resolution, final Particle construction, and pair processing
 * occur after them.
 */

#ifndef SMASH_AFTERBURNER_ANALYSIS_HBT_EVENT_ORIGIN_SELECTOR_H
#define SMASH_AFTERBURNER_ANALYSIS_HBT_EVENT_ORIGIN_SELECTOR_H

#include "hbt/config/hbt_config.h"
#include "hbt/event/origin_flags.h"

namespace hbt {

/**
 * @brief Classify one particle from its two raw mother-PDG fields.
 *
 * Classification preserves the legacy mother-field rule exactly:
 *
 * - both mothers zero:
 *   primordial, primordial+rescattering, and widest selection;
 * - both mothers non-zero:
 *   primordial+rescattering and widest selection;
 * - exactly one mother non-zero:
 *   widest selection only.
 *
 * The numerical PDG identities are not interpreted here. Only the zero/non-zero
 * pattern of the two raw mother fields participates in this classification.
 *
 * @param pdg_mother1 Raw first-mother PDG field from Afterburner.
 * @param pdg_mother2 Raw second-mother PDG field from Afterburner.
 * @return Nested origin flags for the particle.
 */
[[nodiscard]] OriginFlags classify_origin(
    int pdg_mother1,
    int pdg_mother2
) noexcept;

/**
 * @brief Test whether origin flags are eligible for one configured mode.
 *
 * Primordial tests the most restrictive flag. PrimordialRescattering tests the
 * corresponding inclusive flag. PrimordialRescatteringDecay tests the widest
 * flag.
 *
 * OriginMode::All requests all three nested selections simultaneously. For
 * event-preparation eligibility it therefore accepts a particle whenever it
 * belongs to the widest selection; the individual flags remain available for
 * later routing to each requested origin slice.
 *
 * @param flags Nested origin memberships of the particle.
 * @param mode Configured HBT origin mode.
 * @return `true` when the particle is eligible for the requested mode.
 * @throws std::invalid_argument If mode is not a valid OriginMode value.
 */
[[nodiscard]] bool is_origin_eligible(
    const OriginFlags& flags,
    OriginMode mode
);

}  // namespace hbt

#endif  // SMASH_AFTERBURNER_ANALYSIS_HBT_EVENT_ORIGIN_SELECTOR_H
