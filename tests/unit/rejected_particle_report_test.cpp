/**
 * @file rejected_particle_report_test.cpp
 * @brief Unit tests for in-memory rejected-particle reporting.
 */

#include "hbt/reporting/rejected_particle_report.h"

#include <cstdlib>
#include <iostream>
#include <stdexcept>

namespace {

/**
 * @brief Construct one recognizable rejected-particle record.
 * @param id Raw particle identifier.
 * @param reason Numerical rejection reason.
 * @param diagnostic Optional diagnostic quantity.
 * @return Complete rejection record.
 */
hbt::RejectedParticleRecord make_record(
    int id,
    hbt::ParticleRejectionReason reason,
    std::optional<double> diagnostic
) {
    return {
        7U,
        31,
        id,
        -211,
        -1,
        hbt::SpeciesId::PiMinus,
        {5.0, 1.0, 2.0, 3.0},
        {10.0, 4.0, 5.0, 6.0},
        0.138,
        2,
        7.5,
        reason,
        diagnostic,
        std::nullopt
    };
}

/**
 * @brief Verify complete record preservation and aggregate reason counts.
 * @return true when insertion order, fields, and counts are exact.
 */
bool verify_storage_and_counts() {
    hbt::RejectedParticleReport report;

    report.add(make_record(
        628,
        hbt::ParticleRejectionReason::NonPositiveInvariantMassSquared,
        -0.003
    ));
    report.add(make_record(
        629,
        hbt::ParticleRejectionReason::NonFiniteMomentum,
        std::nullopt
    ));
    report.add(make_record(
        630,
        hbt::ParticleRejectionReason::NonPositiveEnergy,
        -1.0
    ));
    report.add(make_record(
        631,
        hbt::ParticleRejectionReason::NonPositiveInvariantMassSquared,
        0.0
    ));

    if (report.size() != 4U || report.empty()) {
        std::cerr
            << "rejected_particle_report_test: wrong total record count.\n";
        return false;
    }

    if (report.count(
            hbt::ParticleRejectionReason::
                NonPositiveInvariantMassSquared) != 2U ||
        report.count(
            hbt::ParticleRejectionReason::NonFiniteMomentum) != 1U ||
        report.count(
            hbt::ParticleRejectionReason::NonPositiveEnergy) != 1U) {
        std::cerr
            << "rejected_particle_report_test: wrong reason counts.\n";
        return false;
    }

    const auto& records = report.records();
    const auto& first = records.front();

    if (first.outer_event_number != 7U || first.subevent_id != 31 ||
        first.particle_id != 628 || first.raw_pdg != -211 ||
        first.raw_charge != -1 || first.species != hbt::SpeciesId::PiMinus ||
        first.momentum.x0 != 5.0 || first.momentum.x3 != 3.0 ||
        first.raw_position.x0 != 10.0 || first.raw_position.x3 != 6.0 ||
        first.raw_mass_gev != 0.138 || first.ncoll != 2 ||
        first.time_last_coll != 7.5 ||
        !first.diagnostic_value.has_value() ||
        first.diagnostic_value.value() != -0.003) {
        std::cerr
            << "rejected_particle_report_test: stored fields changed.\n";
        return false;
    }

    return true;
}

/**
 * @brief Verify the non-physical Count sentinel cannot be stored or queried.
 * @return true when both invalid operations throw std::invalid_argument.
 */
bool verify_count_sentinel_rejected() {
    hbt::RejectedParticleReport report;
    bool add_threw = false;
    bool count_threw = false;

    try {
        report.add(make_record(
            1,
            hbt::ParticleRejectionReason::Count,
            std::nullopt
        ));
    } catch (const std::invalid_argument&) {
        add_threw = true;
    }

    try {
        static_cast<void>(report.count(
            hbt::ParticleRejectionReason::Count
        ));
    } catch (const std::invalid_argument&) {
        count_threw = true;
    }

    return add_threw && count_threw && report.empty();
}

}  // namespace

/**
 * @brief Run all RejectedParticleReport unit tests.
 * @return EXIT_SUCCESS when every test passes, otherwise EXIT_FAILURE.
 */
int main() {
    bool success = true;
    success = verify_storage_and_counts() && success;
    success = verify_count_sentinel_rejected() && success;
    return success ? EXIT_SUCCESS : EXIT_FAILURE;
}
