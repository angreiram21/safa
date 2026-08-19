/**
 * @file pair_origin_routing_test.cpp
 * @brief Unit tests for pair-origin membership and requested-origin routing.
 */

#include "hbt/pair/pair_origin_routing.h"

#include <array>
#include <cstdlib>
#include <iostream>
#include <stdexcept>
#include <string_view>

namespace {

/**
 * @brief Exclusive origin label used only to define the normative test table.
 */
enum class TestOrigin {
    Primordial,   ///< Primordial test particle.
    Rescattering, ///< Rescattering test particle.
    Decay         ///< Decay test particle.
};

/**
 * @brief Convert one exclusive test label to production nested OriginFlags.
 * @param origin Exclusive origin label used by the test matrix.
 * @return Corresponding nested production flags.
 */
hbt::OriginFlags flags_for(TestOrigin origin) {
    switch (origin) {
    case TestOrigin::Primordial:
        return {true, true, true};
    case TestOrigin::Rescattering:
        return {false, true, true};
    case TestOrigin::Decay:
        return {false, false, true};
    }

    throw std::invalid_argument("invalid test origin");
}

/**
 * @brief One row of the complete P/R/D pair-origin membership table.
 */
struct MembershipCase {
    TestOrigin origin_a;  ///< Exclusive origin of particle A.
    TestOrigin origin_b;  ///< Exclusive origin of particle B.
    hbt::PairOriginMemberships expected;  ///< Required pair memberships.
    std::string_view label;  ///< Diagnostic row label.
};

/**
 * @brief Compare pair-origin memberships field by field.
 * @param actual Production memberships.
 * @param expected Required memberships.
 * @return `true` when all memberships are equal.
 */
bool equal_memberships(
    const hbt::PairOriginMemberships& actual,
    const hbt::PairOriginMemberships& expected
) {
    return actual.primordial == expected.primordial &&
           actual.primordial_rescattering ==
               expected.primordial_rescattering &&
           actual.primordial_rescattering_decay ==
               expected.primordial_rescattering_decay;
}

/**
 * @brief Compare requested pair-origin routes field by field.
 * @param actual Production routes.
 * @param expected Required routes.
 * @return `true` when all routes are equal.
 */
bool equal_routes(
    const hbt::PairOriginRoutes& actual,
    const hbt::PairOriginRoutes& expected
) {
    return actual.primordial == expected.primordial &&
           actual.primordial_rescattering ==
               expected.primordial_rescattering &&
           actual.primordial_rescattering_decay ==
               expected.primordial_rescattering_decay;
}

/**
 * @brief Return the normative nine-row P/R/D pair-origin table.
 * @return All ordered particle-origin combinations and required memberships.
 */
std::array<MembershipCase, 9> membership_cases() {
    using Origin = TestOrigin;

    return {{
        {Origin::Primordial,
         Origin::Primordial,
         {true, true, true},
         "P-P"},
        {Origin::Primordial,
         Origin::Rescattering,
         {false, true, true},
         "P-R"},
        {Origin::Rescattering,
         Origin::Primordial,
         {false, true, true},
         "R-P"},
        {Origin::Rescattering,
         Origin::Rescattering,
         {false, true, true},
         "R-R"},
        {Origin::Primordial,
         Origin::Decay,
         {false, false, true},
         "P-D"},
        {Origin::Decay,
         Origin::Primordial,
         {false, false, true},
         "D-P"},
        {Origin::Rescattering,
         Origin::Decay,
         {false, false, true},
         "R-D"},
        {Origin::Decay,
         Origin::Rescattering,
         {false, false, true},
         "D-R"},
        {Origin::Decay,
         Origin::Decay,
         {false, false, true},
         "D-D"}
    }};
}

/**
 * @brief Verify every row of the normative physical-membership table.
 * @return `true` when all nine ordered combinations match exactly.
 */
bool verify_complete_membership_table() {
    for (const MembershipCase& test_case : membership_cases()) {
        const hbt::PairOriginMemberships actual =
            hbt::calculate_pair_origin_memberships(
                flags_for(test_case.origin_a),
                flags_for(test_case.origin_b)
            );

        if (!equal_memberships(actual, test_case.expected)) {
            std::cerr
                << "pair_origin_routing_test: wrong memberships for "
                << test_case.label
                << ".\n";
            return false;
        }
    }

    return true;
}

/**
 * @brief Verify individual modes activate only their requested origin route.
 * @return `true` when all nine table rows obey every individual mode.
 */
bool verify_individual_mode_routing() {
    for (const MembershipCase& test_case : membership_cases()) {
        const hbt::PairOriginMemberships memberships =
            hbt::calculate_pair_origin_memberships(
                flags_for(test_case.origin_a),
                flags_for(test_case.origin_b)
            );

        const std::array<hbt::PairOriginRoutes, 3> expected{{
            {test_case.expected.primordial, false, false},
            {false, test_case.expected.primordial_rescattering, false},
            {false, false,
             test_case.expected.primordial_rescattering_decay}
        }};

        const std::array<hbt::OriginMode, 3> modes{{
            hbt::OriginMode::Primordial,
            hbt::OriginMode::PrimordialRescattering,
            hbt::OriginMode::PrimordialRescatteringDecay
        }};

        for (std::size_t index = 0; index < modes.size(); ++index) {
            const hbt::PairOriginRoutes actual =
                hbt::route_pair_origin_memberships(
                    memberships,
                    modes[index]
                );

            if (!equal_routes(actual, expected[index])) {
                std::cerr
                    << "pair_origin_routing_test: wrong individual routing for "
                    << test_case.label
                    << ".\n";
                return false;
            }
        }
    }

    return true;
}

/**
 * @brief Verify OriginMode::All follows the complete normative table exactly.
 * @return `true` when every pair is routed to all and only compatible slices.
 */
bool verify_all_mode_routes_every_compatible_slice() {
    for (const MembershipCase& test_case : membership_cases()) {
        const hbt::PairOriginMemberships memberships =
            hbt::calculate_pair_origin_memberships(
                flags_for(test_case.origin_a),
                flags_for(test_case.origin_b)
            );
        const hbt::PairOriginRoutes actual =
            hbt::route_pair_origin_memberships(
                memberships,
                hbt::OriginMode::All
            );
        const hbt::PairOriginRoutes expected{
            test_case.expected.primordial,
            test_case.expected.primordial_rescattering,
            test_case.expected.primordial_rescattering_decay
        };

        if (!equal_routes(actual, expected)) {
            std::cerr
                << "pair_origin_routing_test: All routing violated table for "
                << test_case.label
                << ".\n";
            return false;
        }
    }

    return true;
}

/**
 * @brief Verify defensive rejection of an invalid OriginMode value.
 * @return `true` when invalid routing mode throws std::invalid_argument.
 */
bool verify_invalid_mode_is_rejected() {
    try {
        static_cast<void>(hbt::route_pair_origin_memberships(
            {true, true, true},
            static_cast<hbt::OriginMode>(999)
        ));
    } catch (const std::invalid_argument&) {
        return true;
    }

    std::cerr
        << "pair_origin_routing_test: invalid OriginMode was accepted.\n";
    return false;
}

}  // namespace

/**
 * @brief Run the complete pair-origin-routing unit-test collection.
 * @return `EXIT_SUCCESS` when every test passes, otherwise `EXIT_FAILURE`.
 */
int main() {
    bool success = true;

    success = verify_complete_membership_table() && success;
    success = verify_individual_mode_routing() && success;
    success = verify_all_mode_routes_every_compatible_slice() && success;
    success = verify_invalid_mode_is_rejected() && success;

    return success ? EXIT_SUCCESS : EXIT_FAILURE;
}
