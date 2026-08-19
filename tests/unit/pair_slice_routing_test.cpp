/**
 * @file pair_slice_routing_test.cpp
 * @brief Unit tests for validated pair routing to configured kT/mT slices.
 */

#include "hbt/pair/pair_slice_routing.h"

#include <cstdlib>
#include <iostream>
#include <limits>
#include <optional>
#include <stdexcept>
#include <vector>

namespace {

/**
 * @brief Build one slicing axis for concise test setup.
 * @param enabled Whether the axis participates in routing.
 * @param edges Configured bin edges in GeV.
 * @return Axis configuration containing the supplied values.
 */
hbt::PairSlicingAxisConfig axis(
    bool enabled,
    std::vector<double> edges
) {
    return {enabled, std::move(edges)};
}

/**
 * @brief Build pair kinematics with selected kT and mT values.
 * @param kt_gev kT value to route.
 * @param mt_gev mT value to route.
 * @return Complete PairKinematics for slice-routing tests.
 */
hbt::PairKinematics make_kinematics(double kt_gev, double mt_gev) {
    return {{2.0, 0.0, 0.0, 0.0}, 0.0, 0.0, kt_gev, mt_gev};
}

/**
 * @brief Compare a route against expected optional kT/mT indices.
 * @param route Actual route returned by production code.
 * @param expected_kt Expected kT index.
 * @param expected_mt Expected mT index.
 * @param expected_flat Expected already resolved flat slice index.
 * @return `true` when all route indices match exactly.
 */
bool route_matches(
    const std::optional<hbt::PairSliceRoute>& route,
    std::optional<std::size_t> expected_kt,
    std::optional<std::size_t> expected_mt,
    std::size_t expected_flat
) {
    return route.has_value() &&
           route->kt_slice_index == expected_kt &&
           route->mt_slice_index == expected_mt &&
           route->flat_slice_index == expected_flat;
}

/**
 * @brief Verify disabled slicing creates no inclusive dummy route.
 * @return `true` when both-disabled slicing always returns std::nullopt.
 */
bool verify_both_disabled_have_no_route() {
    const hbt::PairSlicingConfig slicing{
        axis(false, {0.2, 0.4, 0.6}),
        axis(false, {0.5, 0.8, 1.1})
    };
    const double nan = std::numeric_limits<double>::quiet_NaN();

    if (hbt::route_pair_to_slices(
            make_kinematics(nan, nan),
            slicing
        ).has_value()) {
        std::cerr
            << "pair_slice_routing_test: disabled slicing made a route.\n";
        return false;
    }

    return true;
}

/**
 * @brief Verify kT-only half-open boundaries and out-of-range behavior.
 * @return `true` when every tested kT value follows the configured bins.
 */
bool verify_kt_only_half_open_routing() {
    const hbt::PairSlicingConfig slicing{
        axis(true, {0.2, 0.4, 0.6}),
        axis(false, {0.5, 0.8, 1.1})
    };
    const double nan = std::numeric_limits<double>::quiet_NaN();

    const struct {
        double kt;
        std::optional<std::size_t> expected;
    } cases[] = {
        {0.19, std::nullopt},
        {0.20, 0U},
        {0.399, 0U},
        {0.40, 1U},
        {0.599, 1U},
        {0.60, std::nullopt},
        {0.90, std::nullopt}
    };

    for (const auto& test_case : cases) {
        const auto route = hbt::route_pair_to_slices(
            make_kinematics(test_case.kt, nan),
            slicing
        );

        if (!test_case.expected.has_value()) {
            if (route.has_value()) {
                std::cerr
                    << "pair_slice_routing_test: out-of-range kT routed.\n";
                return false;
            }
            continue;
        }

        if (!route_matches(
                route, test_case.expected, std::nullopt,
                test_case.expected.value()
            )) {
            std::cerr
                << "pair_slice_routing_test: wrong kT-only slice.\n";
            return false;
        }
    }

    return true;
}

/**
 * @brief Verify mT-only routing ignores the disabled kT axis completely.
 * @return `true` when mT routes correctly with non-finite unused kT.
 */
bool verify_mt_only_routing() {
    const hbt::PairSlicingConfig slicing{
        axis(false, {0.2, 0.4, 0.6}),
        axis(true, {0.5, 0.8, 1.1})
    };
    const double nan = std::numeric_limits<double>::quiet_NaN();

    if (!route_matches(
            hbt::route_pair_to_slices(make_kinematics(nan, 0.50), slicing),
            std::nullopt,
            0U,
            0U)) {
        std::cerr
            << "pair_slice_routing_test: wrong first mT slice.\n";
        return false;
    }
    if (!route_matches(
            hbt::route_pair_to_slices(make_kinematics(nan, 0.80), slicing),
            std::nullopt,
            1U,
            1U)) {
        std::cerr
            << "pair_slice_routing_test: wrong second mT slice.\n";
        return false;
    }
    if (hbt::route_pair_to_slices(
            make_kinematics(nan, 1.10),
            slicing
        ).has_value()) {
        std::cerr
            << "pair_slice_routing_test: upper mT edge was included.\n";
        return false;
    }

    return true;
}

/**
 * @brief Verify Cartesian routing requires membership in both enabled axes.
 * @return `true` when one unique kT x mT cell is returned only in range.
 */
bool verify_cartesian_routing() {
    const hbt::PairSlicingConfig slicing{
        axis(true, {0.2, 0.4, 0.6}),
        axis(true, {0.5, 0.8, 1.1})
    };

    if (!route_matches(
            hbt::route_pair_to_slices(make_kinematics(0.45, 0.90), slicing),
            1U,
            1U,
            3U)) {
        std::cerr
            << "pair_slice_routing_test: wrong Cartesian cell.\n";
        return false;
    }
    if (hbt::route_pair_to_slices(
            make_kinematics(0.45, 1.10),
            slicing
        ).has_value()) {
        std::cerr
            << "pair_slice_routing_test: mT-out pair got a cell.\n";
        return false;
    }
    if (hbt::route_pair_to_slices(
            make_kinematics(0.60, 0.90),
            slicing
        ).has_value()) {
        std::cerr
            << "pair_slice_routing_test: kT-out pair got a cell.\n";
        return false;
    }

    return true;
}

/**
 * @brief Verify direct misuse with non-finite enabled kinematics is rejected.
 * @return `true` when required non-finite kinematics throw invalid_argument.
 */
bool verify_nonfinite_enabled_axis_is_rejected() {
    const double nan = std::numeric_limits<double>::quiet_NaN();
    const hbt::PairSlicingConfig slicing{
        axis(true, {0.2, 0.4}),
        axis(false, {})
    };

    try {
        static_cast<void>(hbt::route_pair_to_slices(
            make_kinematics(nan, 0.7),
            slicing
        ));
    } catch (const std::invalid_argument&) {
        return true;
    }

    std::cerr
        << "pair_slice_routing_test: non-finite enabled kT was accepted.\n";
    return false;
}

/**
 * @brief Verify a too-short enabled axis is rejected defensively.
 * @return `true` when fewer than two enabled edges throw invalid_argument.
 */
bool verify_undersized_enabled_axis_is_rejected() {
    const hbt::PairSlicingConfig slicing{
        axis(true, {0.2}),
        axis(false, {})
    };

    try {
        static_cast<void>(hbt::route_pair_to_slices(
            make_kinematics(0.3, 0.7),
            slicing
        ));
    } catch (const std::invalid_argument&) {
        return true;
    }

    std::cerr
        << "pair_slice_routing_test: undersized enabled axis accepted.\n";
    return false;
}

/**
 * @brief Verify binary lookup preserves half-open semantics over many bins.
 * @return `true` when interior edges route to the bin beginning at that edge.
 */
bool verify_many_bin_binary_lookup_boundaries() {
    const hbt::PairSlicingConfig slicing{
        axis(true, {0.0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7}),
        axis(false, {})
    };

    for (std::size_t index = 0U; index < 7U; ++index) {
        const double value = slicing.kt.bin_edges_gev[index];
        if (!route_matches(
                hbt::route_pair_to_slices(make_kinematics(value, 0.0), slicing),
                index,
                std::nullopt,
                index)) {
            std::cerr
                << "pair_slice_routing_test: binary boundary mismatch.\n";
            return false;
        }
    }

    if (hbt::route_pair_to_slices(
            make_kinematics(0.7, 0.0),
            slicing
        ).has_value()) {
        std::cerr
            << "pair_slice_routing_test: final edge was included.\n";
        return false;
    }

    return true;
}

}  // namespace

/**
 * @brief Run the complete pair-slice-routing unit-test collection.
 * @return `EXIT_SUCCESS` when every test passes, otherwise `EXIT_FAILURE`.
 */
int main() {
    if (!verify_both_disabled_have_no_route() ||
        !verify_kt_only_half_open_routing() ||
        !verify_mt_only_routing() ||
        !verify_cartesian_routing() ||
        !verify_nonfinite_enabled_axis_is_rejected() ||
        !verify_undersized_enabled_axis_is_rejected() ||
        !verify_many_bin_binary_lookup_boundaries()) {
        return EXIT_FAILURE;
    }

    return EXIT_SUCCESS;
}
