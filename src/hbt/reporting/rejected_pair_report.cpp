/**
 * @file rejected_pair_report.cpp
 * @brief Implementation of rejected-pair in-memory storage.
 */

#include "hbt/reporting/rejected_pair_report.h"

#include <stdexcept>
#include <utility>

namespace hbt {
namespace {

/**
 * @brief Convert a physical pair-rejection reason to its array index.
 * @param reason Pair-rejection reason to convert.
 * @return Zero-based physical-reason index.
 * @throws std::invalid_argument If reason is Count or invalid.
 */
std::size_t reason_index(PairRejectionReason reason) {
    const std::size_t index = static_cast<std::size_t>(reason);
    const std::size_t count =
        static_cast<std::size_t>(PairRejectionReason::Count);

    if (index >= count) {
        throw std::invalid_argument(
            "RejectedPairReport: invalid rejection reason"
        );
    }

    return index;
}

}  // namespace

RejectedPairParticleSnapshot make_rejected_pair_particle_snapshot(
    const Particle& particle
) noexcept {
    return {
        particle.species,
        particle.momentum,
        particle.invariant_mass_gev,
        particle.raw_pdg,
        particle.raw_charge
    };
}

void RejectedPairReport::add(RejectedPairRecord record) {
    const std::size_t index = reason_index(record.reason);
    records_.push_back(std::move(record));
    ++counts_[index];
}

const std::vector<RejectedPairRecord>& RejectedPairReport::records()
    const noexcept {
    return records_;
}

std::size_t RejectedPairReport::size() const noexcept {
    return records_.size();
}

bool RejectedPairReport::empty() const noexcept {
    return records_.empty();
}

std::size_t RejectedPairReport::count(PairRejectionReason reason) const {
    return counts_[reason_index(reason)];
}

}  // namespace hbt
