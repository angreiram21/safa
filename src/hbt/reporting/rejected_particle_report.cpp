/**
 * @file rejected_particle_report.cpp
 * @brief Implementation of rejected-particle in-memory storage.
 */

#include "hbt/reporting/rejected_particle_report.h"

#include <stdexcept>
#include <utility>

namespace hbt {
namespace {

/**
 * @brief Convert a physical rejection reason to its array index.
 * @param reason Rejection reason to convert.
 * @return Zero-based physical-reason index.
 * @throws std::invalid_argument If reason is Count or invalid.
 */
std::size_t reason_index(ParticleRejectionReason reason) {
    const std::size_t index = static_cast<std::size_t>(reason);
    const std::size_t count =
        static_cast<std::size_t>(ParticleRejectionReason::Count);

    if (index >= count) {
        throw std::invalid_argument(
            "RejectedParticleReport: invalid rejection reason"
        );
    }

    return index;
}

}  // namespace

void RejectedParticleReport::add(RejectedParticleRecord record) {
    const std::size_t index = reason_index(record.reason);
    records_.push_back(std::move(record));
    ++counts_[index];
}

const std::vector<RejectedParticleRecord>&
RejectedParticleReport::records() const noexcept {
    return records_;
}

std::size_t RejectedParticleReport::size() const noexcept {
    return records_.size();
}

bool RejectedParticleReport::empty() const noexcept {
    return records_.empty();
}

std::size_t RejectedParticleReport::count(
    ParticleRejectionReason reason) const {
    return counts_[reason_index(reason)];
}

}  // namespace hbt
