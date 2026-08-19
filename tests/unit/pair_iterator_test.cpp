/**
 * @file pair_iterator_test.cpp
 * @brief Unit tests for streaming primitive-channel pair iteration.
 */

#include "hbt/pair/pair_iterator.h"

#include <algorithm>
#include <cstdint>
#include <cstdlib>
#include <iostream>
#include <stdexcept>
#include <string>
#include <utility>
#include <vector>

namespace {

/**
 * @brief Compact particle-marker pair used to inspect iterator output.
 */
struct MarkerPair {
    double marker_a;  ///< Marker from the particle occupying canonical role A.
    double marker_b;  ///< Marker from the particle occupying canonical role B.
};

/**
 * @brief Construct one accepted particle with a recognizable marker.
 * @param species Canonical particle species.
 * @param marker Value stored in position.x0 for pair-identification checks.
 * @return Complete particle suitable for EventBuffers tests.
 */
hbt::Particle make_particle(hbt::SpeciesId species, double marker) {
    return {
        species,
        {marker, marker + 0.1, marker + 0.2, marker + 0.3},
        {marker + 1.0, marker + 1.1, marker + 1.2, marker + 1.3},
        marker + 2.0,
        {true, true, true},
        211,
        1
    };
}

/**
 * @brief Report one failed unit-test condition.
 * @param message Description of the violated iterator contract.
 * @return Always false, allowing direct return from a failed test.
 */
bool fail(const char* message) {
    std::cerr << "pair_iterator_test: " << message << ".\n";
    return false;
}

/**
 * @brief Verify an empty required-channel list performs no work.
 * @return true when the summary is empty and the consumer is never invoked.
 */
bool verify_empty_channel_list() {
    const hbt::EventBuffers buffers;
    std::uint64_t calls = 0U;

    const hbt::PairCountSummary summary = hbt::for_each_pair(
        buffers,
        {},
        [&calls](
            std::size_t,
            hbt::PrimitiveChannelId,
            const hbt::Particle&,
            const hbt::Particle&
        ) {
            ++calls;
        }
    );

    if (!summary.channels.empty() || calls != 0U) {
        return fail("empty channel list produced pair work");
    }
    return true;
}

/**
 * @brief Verify identical-channel counts for zero, one, two and four particles.
 * @return true when each N produces N(N-1)/2 delivered pairs.
 */
bool verify_identical_pair_counts() {
    /**
     * @brief Expected pair count for one identical-channel multiplicity.
     */
    struct CountCase {
        /// Number of particles supplied to the iterator.
        std::size_t particle_count;
        /// Expected number of unique unordered pairs.
        std::uint64_t expected_pairs;
    };

    const std::vector<CountCase> cases{
        {0U, 0U},
        {1U, 0U},
        {2U, 1U},
        {4U, 6U}
    };

    for (const CountCase& test_case : cases) {
        hbt::EventBuffers buffers;
        for (std::size_t index = 0; index < test_case.particle_count; ++index) {
            buffers.add(
                make_particle(
                    hbt::SpeciesId::PiPlus,
                    static_cast<double>(index + 1U)
                )
            );
        }

        std::uint64_t calls = 0U;
        const hbt::PairCountSummary summary = hbt::for_each_pair(
            buffers,
            {hbt::PrimitiveChannelId::PiPlusPiPlus},
            [&calls](
                std::size_t,
                hbt::PrimitiveChannelId,
                const hbt::Particle&,
                const hbt::Particle&
            ) {
                ++calls;
            }
        );

        if (summary.channels.size() != 1U) {
            return fail("identical channel summary has incorrect size");
        }
        if (summary.channels[0].pair_count != test_case.expected_pairs ||
            calls != test_case.expected_pairs) {
            return fail("identical channel pair count is incorrect");
        }
    }
    return true;
}

/**
 * @brief Verify identical particles are paired once with the lower index first.
 * @return true when the six expected i<j combinations are delivered in order.
 */
bool verify_identical_pair_membership() {
    hbt::EventBuffers buffers;
    for (int marker = 1; marker <= 4; ++marker) {
        buffers.add(
            make_particle(hbt::SpeciesId::PiPlus, static_cast<double>(marker))
        );
    }

    std::vector<MarkerPair> actual;
    const hbt::PairCountSummary summary = hbt::for_each_pair(
        buffers,
        {hbt::PrimitiveChannelId::PiPlusPiPlus},
        [&actual](
            std::size_t channel_index,
            hbt::PrimitiveChannelId channel,
            const hbt::Particle& particle_a,
            const hbt::Particle& particle_b
        ) {
            if (channel_index != 0U ||
                channel != hbt::PrimitiveChannelId::PiPlusPiPlus) {
                throw std::runtime_error("unexpected channel/index");
            }
            actual.push_back({particle_a.position.x0, particle_b.position.x0});
        }
    );

    const std::vector<MarkerPair> expected{
        {1.0, 2.0},
        {1.0, 3.0},
        {1.0, 4.0},
        {2.0, 3.0},
        {2.0, 4.0},
        {3.0, 4.0}
    };

    if (summary.channels[0].pair_count != expected.size()) {
        return fail("identical membership count is incorrect");
    }
    if (actual.size() != expected.size()) {
        return fail("identical membership delivered incorrect pair total");
    }

    for (std::size_t index = 0; index < expected.size(); ++index) {
        if (actual[index].marker_a != expected[index].marker_a ||
            actual[index].marker_b != expected[index].marker_b) {
            return fail("identical channel did not preserve i<j membership");
        }
    }
    return true;
}

/**
 * @brief Compare marker pairs lexicographically for deterministic set checks.
 * @param left First marker pair.
 * @param right Second marker pair.
 * @return true when left sorts before right.
 */
bool marker_pair_less(const MarkerPair& left, const MarkerPair& right) {
    if (left.marker_a != right.marker_a) {
        return left.marker_a < right.marker_a;
    }
    return left.marker_b < right.marker_b;
}

/**
 * @brief Verify a cross-species channel produces the complete A x B product.
 * @return true when all six expected pairs appear exactly once with canonical
 *         A/B species roles.
 */
bool verify_cross_species_product_and_roles() {
    hbt::EventBuffers buffers;
    buffers.add(make_particle(hbt::SpeciesId::KMinus, 1.0));
    buffers.add(make_particle(hbt::SpeciesId::KMinus, 2.0));
    buffers.add(make_particle(hbt::SpeciesId::KMinus, 3.0));
    buffers.add(make_particle(hbt::SpeciesId::Proton, 10.0));
    buffers.add(make_particle(hbt::SpeciesId::Proton, 20.0));

    std::vector<MarkerPair> actual;
    bool roles_are_canonical = true;
    const hbt::PairCountSummary summary = hbt::for_each_pair(
        buffers,
        {hbt::PrimitiveChannelId::KMinusProton},
        [&actual, &roles_are_canonical](
            std::size_t channel_index,
            hbt::PrimitiveChannelId channel,
            const hbt::Particle& particle_a,
            const hbt::Particle& particle_b
        ) {
            if (channel_index != 0U ||
                channel != hbt::PrimitiveChannelId::KMinusProton ||
                particle_a.species != hbt::SpeciesId::KMinus ||
                particle_b.species != hbt::SpeciesId::Proton) {
                roles_are_canonical = false;
            }
            actual.push_back({particle_a.position.x0, particle_b.position.x0});
        }
    );

    std::vector<MarkerPair> expected{
        {1.0, 10.0},
        {1.0, 20.0},
        {2.0, 10.0},
        {2.0, 20.0},
        {3.0, 10.0},
        {3.0, 20.0}
    };

    std::sort(actual.begin(), actual.end(), marker_pair_less);
    std::sort(expected.begin(), expected.end(), marker_pair_less);

    if (!roles_are_canonical) {
        return fail("cross-species canonical A/B roles were not preserved");
    }
    if (summary.channels[0].pair_count != expected.size() ||
        actual.size() != expected.size()) {
        return fail("cross-species Cartesian-product count is incorrect");
    }

    for (std::size_t index = 0; index < expected.size(); ++index) {
        if (actual[index].marker_a != expected[index].marker_a ||
            actual[index].marker_b != expected[index].marker_b) {
            return fail("cross-species Cartesian product is incomplete");
        }
    }
    return true;
}

/**
 * @brief Verify cross-species channels with an empty side produce zero pairs.
 * @return true when both empty-side cases retain a zero-count summary entry.
 */
bool verify_cross_species_empty_side() {
    hbt::EventBuffers only_a;
    only_a.add(make_particle(hbt::SpeciesId::KMinus, 1.0));

    hbt::EventBuffers only_b;
    only_b.add(make_particle(hbt::SpeciesId::Proton, 2.0));

    const auto no_op = [](
        std::size_t,
        hbt::PrimitiveChannelId,
        const hbt::Particle&,
        const hbt::Particle&
    ) {};

    const hbt::PairCountSummary missing_b = hbt::for_each_pair(
        only_a,
        {hbt::PrimitiveChannelId::KMinusProton},
        no_op
    );
    const hbt::PairCountSummary missing_a = hbt::for_each_pair(
        only_b,
        {hbt::PrimitiveChannelId::KMinusProton},
        no_op
    );

    if (missing_b.channels[0].pair_count != 0U ||
        missing_a.channels[0].pair_count != 0U) {
        return fail("cross-species empty buffer did not produce zero pairs");
    }
    return true;
}

/**
 * @brief Verify summary entries preserve required-channel order.
 * @return true when all channel IDs and pair counts remain in input order.
 */
bool verify_multiple_channel_summary_order() {
    hbt::EventBuffers buffers;
    buffers.add(make_particle(hbt::SpeciesId::PiPlus, 1.0));
    buffers.add(make_particle(hbt::SpeciesId::PiPlus, 2.0));
    buffers.add(make_particle(hbt::SpeciesId::KMinus, 3.0));
    buffers.add(make_particle(hbt::SpeciesId::Proton, 4.0));
    buffers.add(make_particle(hbt::SpeciesId::Proton, 5.0));

    const std::vector<hbt::PrimitiveChannelId> channels{
        hbt::PrimitiveChannelId::PiPlusPiPlus,
        hbt::PrimitiveChannelId::KMinusProton,
        hbt::PrimitiveChannelId::ProtonProton
    };

    std::vector<std::size_t> delivered_indices;
    const hbt::PairCountSummary summary = hbt::for_each_pair(
        buffers,
        channels,
        [&delivered_indices](
            std::size_t channel_index,
            hbt::PrimitiveChannelId,
            const hbt::Particle&,
            const hbt::Particle&
        ) {
            delivered_indices.push_back(channel_index);
        }
    );

    if (summary.channels.size() != channels.size()) {
        return fail("multiple-channel summary has incorrect size");
    }

    const std::vector<std::uint64_t> expected_counts{1U, 2U, 1U};
    const std::vector<std::size_t> expected_indices{0U, 1U, 1U, 2U};
    if (delivered_indices != expected_indices) {
        return fail("consumer channel indices are incorrect");
    }

    for (std::size_t index = 0; index < channels.size(); ++index) {
        if (summary.channels[index].channel != channels[index] ||
            summary.channels[index].pair_count != expected_counts[index]) {
            return fail("multiple-channel summary order is incorrect");
        }
    }
    return true;
}

/**
 * @brief Verify a repeated required channel is rejected before any pair work.
 * @return true when std::invalid_argument is thrown and the consumer is unused.
 */
bool verify_duplicate_channel_rejected() {
    hbt::EventBuffers buffers;
    buffers.add(make_particle(hbt::SpeciesId::PiPlus, 1.0));
    buffers.add(make_particle(hbt::SpeciesId::PiPlus, 2.0));
    std::uint64_t calls = 0U;

    try {
        static_cast<void>(hbt::for_each_pair(
            buffers,
            {
                hbt::PrimitiveChannelId::PiPlusPiPlus,
                hbt::PrimitiveChannelId::PiPlusPiPlus
            },
            [&calls](
                std::size_t,
                hbt::PrimitiveChannelId,
                const hbt::Particle&,
                const hbt::Particle&
            ) {
                ++calls;
            }
        ));
    } catch (const std::invalid_argument&) {
        return calls == 0U;
    } catch (...) {
        return fail("duplicate channel produced unexpected exception type");
    }

    return fail("duplicate channel was not rejected");
}

/**
 * @brief Verify an invalid channel ID is rejected before any pair work.
 * @return true when std::invalid_argument is thrown and the consumer is unused.
 */
bool verify_invalid_channel_rejected() {
    hbt::EventBuffers buffers;
    buffers.add(make_particle(hbt::SpeciesId::PiPlus, 1.0));
    buffers.add(make_particle(hbt::SpeciesId::PiPlus, 2.0));
    const auto invalid = static_cast<hbt::PrimitiveChannelId>(999);
    std::uint64_t calls = 0U;

    try {
        static_cast<void>(hbt::for_each_pair(
            buffers,
            {hbt::PrimitiveChannelId::PiPlusPiPlus, invalid},
            [&calls](
                std::size_t,
                hbt::PrimitiveChannelId,
                const hbt::Particle&,
                const hbt::Particle&
            ) {
                ++calls;
            }
        ));
    } catch (const std::invalid_argument&) {
        return calls == 0U;
    } catch (...) {
        return fail("invalid channel produced unexpected exception type");
    }

    return fail("invalid channel was not rejected");
}

/**
 * @brief Verify pair-consumer exceptions are not caught by the iterator.
 * @return true when the exact std::runtime_error reaches the caller.
 */
bool verify_consumer_error_propagates() {
    hbt::EventBuffers buffers;
    buffers.add(make_particle(hbt::SpeciesId::PiPlus, 1.0));
    buffers.add(make_particle(hbt::SpeciesId::PiPlus, 2.0));

    try {
        static_cast<void>(hbt::for_each_pair(
            buffers,
            {hbt::PrimitiveChannelId::PiPlusPiPlus},
            [](
                std::size_t,
                hbt::PrimitiveChannelId,
                const hbt::Particle&,
                const hbt::Particle&
            ) {
                throw std::runtime_error("consumer failure");
            }
        ));
    } catch (const std::runtime_error& error) {
        return std::string(error.what()) == "consumer failure";
    } catch (...) {
        return fail("consumer error changed exception type");
    }

    return fail("consumer error was suppressed");
}

}  // namespace

/**
 * @brief Run the complete pair-iterator unit-test collection.
 * @return EXIT_SUCCESS when every test passes, otherwise EXIT_FAILURE.
 */
int main() {
    bool success = true;
    success = verify_empty_channel_list() && success;
    success = verify_identical_pair_counts() && success;
    success = verify_identical_pair_membership() && success;
    success = verify_cross_species_product_and_roles() && success;
    success = verify_cross_species_empty_side() && success;
    success = verify_multiple_channel_summary_order() && success;
    success = verify_duplicate_channel_rejected() && success;
    success = verify_invalid_channel_rejected() && success;
    success = verify_consumer_error_propagates() && success;
    return success ? EXIT_SUCCESS : EXIT_FAILURE;
}
