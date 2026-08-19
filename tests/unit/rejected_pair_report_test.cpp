/**
 * @file rejected_pair_report_test.cpp
 * @brief Unit tests for in-memory numerical pair-rejection reporting.
 */

#include "hbt/reporting/rejected_pair_report.h"

#include <cstdlib>
#include <iostream>
#include <stdexcept>

namespace {

/**
 * @brief Report one failed rejected-pair report condition.
 * @param message Description of the violated contract.
 * @return Always false.
 */
bool fail(const char* message) {
    std::cerr << "rejected_pair_report_test: " << message << ".\n";
    return false;
}

/**
 * @brief Construct one accepted particle with recognizable diagnostic data.
 * @param species Canonical HBT species.
 * @param marker Numeric marker copied into the particle data.
 * @return Complete accepted particle for reporting tests.
 */
hbt::Particle make_particle(hbt::SpeciesId species, double marker) {
    return {
        species,
        {marker, marker + 0.1, marker + 0.2, marker + 0.3},
        {marker + 1.0, marker + 1.1, marker + 1.2, marker + 1.3},
        marker + 2.0,
        {true, true, true},
        species == hbt::SpeciesId::PiPlus ? 211 : 2212,
        1
    };
}

/**
 * @brief Construct one complete rejected-pair record.
 * @param reason Numerical reason assigned to the rejected pair.
 * @param ordinal One-based pair ordinal within the primitive channel.
 * @return Complete record suitable for report tests.
 */
hbt::RejectedPairRecord make_record(
    hbt::PairRejectionReason reason,
    std::uint64_t ordinal
) {
    const hbt::Particle particle_a =
        make_particle(hbt::SpeciesId::PiPlus, 1.0);
    const hbt::Particle particle_b =
        make_particle(hbt::SpeciesId::Proton, 2.0);

    return {
        3U,
        17,
        hbt::PrimitiveChannelId::PiPlusProton,
        ordinal,
        hbt::make_rejected_pair_particle_snapshot(particle_a),
        hbt::make_rejected_pair_particle_snapshot(particle_b),
        {{9.0, 1.0, 2.0, 3.0}, 0.5, 1.0, 4.5, 6.5},
        reason
    };
}

/**
 * @brief Verify particle snapshots preserve pair-kinematics input diagnostics.
 * @return true when every selected particle field is copied exactly.
 */
bool verify_particle_snapshot() {
    const hbt::Particle particle =
        make_particle(hbt::SpeciesId::PiPlus, 4.0);
    const hbt::RejectedPairParticleSnapshot snapshot =
        hbt::make_rejected_pair_particle_snapshot(particle);

    if (snapshot.species != particle.species ||
        snapshot.momentum.x0 != particle.momentum.x0 ||
        snapshot.momentum.x1 != particle.momentum.x1 ||
        snapshot.momentum.x2 != particle.momentum.x2 ||
        snapshot.momentum.x3 != particle.momentum.x3 ||
        snapshot.invariant_mass_gev != particle.invariant_mass_gev ||
        snapshot.raw_pdg != particle.raw_pdg ||
        snapshot.raw_charge != particle.raw_charge) {
        return fail("particle snapshot did not preserve diagnostic data");
    }
    return true;
}

/**
 * @brief Verify complete records and exact rejection-reason counts are stored.
 * @return true when insertion order, fields, size and counts are preserved.
 */
bool verify_storage_and_counts() {
    hbt::RejectedPairReport report;
    report.add(make_record(hbt::PairRejectionReason::NonFiniteKt, 4U));
    report.add(make_record(hbt::PairRejectionReason::NonFiniteMt, 9U));
    report.add(make_record(hbt::PairRejectionReason::NonFiniteKt, 12U));

    if (report.empty() || report.size() != 3U) {
        return fail("stored rejected-pair size is incorrect");
    }
    if (report.count(hbt::PairRejectionReason::NonFiniteKt) != 2U ||
        report.count(hbt::PairRejectionReason::NonFiniteMt) != 1U) {
        return fail("pair-rejection reason counts are incorrect");
    }

    const auto& records = report.records();
    if (records[0].pair_ordinal_in_channel != 4U ||
        records[1].pair_ordinal_in_channel != 9U ||
        records[2].pair_ordinal_in_channel != 12U) {
        return fail("rejected-pair insertion order was not preserved");
    }
    if (records[0].outer_event_number != 3U ||
        records[0].subevent_id != 17 ||
        records[0].channel != hbt::PrimitiveChannelId::PiPlusProton ||
        records[0].kinematics.kt_gev != 4.5 ||
        records[0].kinematics.mt_gev != 6.5) {
        return fail("complete rejected-pair record was not preserved");
    }
    return true;
}

/**
 * @brief Verify an empty report exposes zero records and zero reason counts.
 * @return true when the empty-state contract is preserved.
 */
bool verify_empty_report() {
    const hbt::RejectedPairReport report;

    if (!report.empty() || report.size() != 0U || !report.records().empty()) {
        return fail("new rejected-pair report is not empty");
    }
    if (report.count(hbt::PairRejectionReason::NonFiniteKt) != 0U ||
        report.count(hbt::PairRejectionReason::NonFiniteMt) != 0U) {
        return fail("empty rejected-pair report has nonzero reason counts");
    }
    return true;
}

/**
 * @brief Verify Count and out-of-range reasons are rejected structurally.
 * @return true when both invalid reason forms throw std::invalid_argument.
 */
bool verify_invalid_reasons_are_rejected() {
    hbt::RejectedPairReport report;

    try {
        report.add(make_record(hbt::PairRejectionReason::Count, 1U));
        return fail("PairRejectionReason::Count was accepted");
    } catch (const std::invalid_argument&) {
    } catch (...) {
        return fail("Count reason produced unexpected exception type");
    }

    const auto invalid = static_cast<hbt::PairRejectionReason>(999);
    try {
        static_cast<void>(report.count(invalid));
        return fail("invalid pair-rejection reason was accepted");
    } catch (const std::invalid_argument&) {
    } catch (...) {
        return fail("invalid reason produced unexpected exception type");
    }

    return true;
}

}  // namespace

/**
 * @brief Run all rejected-pair report unit tests.
 * @return EXIT_SUCCESS when every test passes, otherwise EXIT_FAILURE.
 */
int main() {
    bool success = true;
    success = verify_particle_snapshot() && success;
    success = verify_storage_and_counts() && success;
    success = verify_empty_report() && success;
    success = verify_invalid_reasons_are_rejected() && success;
    return success ? EXIT_SUCCESS : EXIT_FAILURE;
}
