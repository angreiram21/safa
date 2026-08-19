/**
 * @file origin_selector_test.cpp
 * @brief Unit tests for nested HBT origin classification and mode eligibility.
 */

#include "hbt/event/origin_selector.h"

#include <array>
#include <cstdlib>
#include <iostream>
#include <stdexcept>
#include <string_view>

namespace {

/**
 * @brief One expected origin classification for a mother-field pair.
 */
struct ClassificationCase {
    int pdg_mother1;       ///< Raw first-mother PDG field.
    int pdg_mother2;       ///< Raw second-mother PDG field.
    hbt::OriginFlags flags;  ///< Expected nested origin flags.
    std::string_view label;  ///< Diagnostic label for the case.
};

/**
 * @brief Compare two OriginFlags values field by field.
 * @param actual Flags returned by production code.
 * @param expected Expected flags for the test case.
 * @return `true` when all three memberships are equal.
 */
bool equal_flags(
    const hbt::OriginFlags& actual,
    const hbt::OriginFlags& expected
) {
    return actual.primordial == expected.primordial &&
           actual.primordial_rescattering ==
               expected.primordial_rescattering &&
           actual.primordial_rescattering_decay ==
               expected.primordial_rescattering_decay;
}

/**
 * @brief Verify the complete legacy mother-field classification matrix.
 * @return `true` when every representative mother pattern is classified
 *         correctly.
 */
bool verify_origin_classification() {
    const std::array<ClassificationCase, 6> cases{{
        {0, 0, {true, true, true}, "both mothers zero"},
        {211, -211, {false, true, true}, "both mothers non-zero"},
        {-2212, 2112, {false, true, true}, "negative non-zero mother"},
        {113, 0, {false, false, true}, "first mother only"},
        {0, 113, {false, false, true}, "second mother only"},
        {-113, 0, {false, false, true}, "negative single mother"}
    }};

    for (const ClassificationCase& test_case : cases) {
        const hbt::OriginFlags actual = hbt::classify_origin(
            test_case.pdg_mother1,
            test_case.pdg_mother2
        );

        if (!equal_flags(actual, test_case.flags)) {
            std::cerr
                << "origin_selector_test: wrong classification for "
                << test_case.label
                << ".\n";
            return false;
        }
    }

    return true;
}

/**
 * @brief One expected eligibility result for a mode and fixed origin flags.
 */
struct EligibilityCase {
    hbt::OriginMode mode;   ///< Origin mode being tested.
    bool expected;          ///< Expected eligibility result.
    std::string_view label;  ///< Diagnostic label for the case.
};

/**
 * @brief Verify one complete OriginMode eligibility row.
 * @param flags Origin flags supplied to the production query.
 * @param expected Expected results for all four OriginMode values.
 * @return `true` when every mode produces the expected result.
 */
bool verify_eligibility_row(
    const hbt::OriginFlags& flags,
    const std::array<EligibilityCase, 4>& expected
) {
    for (const EligibilityCase& test_case : expected) {
        if (hbt::is_origin_eligible(flags, test_case.mode) !=
            test_case.expected) {
            std::cerr
                << "origin_selector_test: wrong eligibility for "
                << test_case.label
                << ".\n";
            return false;
        }
    }

    return true;
}

/**
 * @brief Verify eligibility of a primordial particle for every origin mode.
 * @return `true` when all four modes accept the primordial particle.
 */
bool verify_primordial_eligibility() {
    const hbt::OriginFlags flags = hbt::classify_origin(0, 0);

    const std::array<EligibilityCase, 4> expected{{
        {hbt::OriginMode::Primordial, true, "primordial / primordial"},
        {hbt::OriginMode::PrimordialRescattering,
         true,
         "primordial / primordial_rescattering"},
        {hbt::OriginMode::PrimordialRescatteringDecay,
         true,
         "primordial / widest"},
        {hbt::OriginMode::All, true, "primordial / all"}
    }};

    return verify_eligibility_row(flags, expected);
}

/**
 * @brief Verify eligibility of a rescattering particle for every origin mode.
 * @return `true` when only inclusive modes accept the particle.
 */
bool verify_rescattering_eligibility() {
    const hbt::OriginFlags flags = hbt::classify_origin(211, -211);

    const std::array<EligibilityCase, 4> expected{{
        {hbt::OriginMode::Primordial, false, "rescattering / primordial"},
        {hbt::OriginMode::PrimordialRescattering,
         true,
         "rescattering / primordial_rescattering"},
        {hbt::OriginMode::PrimordialRescatteringDecay,
         true,
         "rescattering / widest"},
        {hbt::OriginMode::All, true, "rescattering / all"}
    }};

    return verify_eligibility_row(flags, expected);
}

/**
 * @brief Verify eligibility of a widest-only particle for every origin mode.
 * @return `true` when only the widest mode and All accept the particle.
 */
bool verify_widest_only_eligibility() {
    const hbt::OriginFlags flags = hbt::classify_origin(113, 0);

    const std::array<EligibilityCase, 4> expected{{
        {hbt::OriginMode::Primordial, false, "widest-only / primordial"},
        {hbt::OriginMode::PrimordialRescattering,
         false,
         "widest-only / primordial_rescattering"},
        {hbt::OriginMode::PrimordialRescatteringDecay,
         true,
         "widest-only / widest"},
        {hbt::OriginMode::All, true, "widest-only / all"}
    }};

    return verify_eligibility_row(flags, expected);
}

/**
 * @brief Verify that All is an eligibility union rather than a fourth flag.
 * @return `true` when All reads only the widest nested membership.
 */
bool verify_all_uses_widest_membership() {
    const hbt::OriginFlags no_membership{false, false, false};
    const hbt::OriginFlags widest_membership{false, false, true};

    if (hbt::is_origin_eligible(
            no_membership,
            hbt::OriginMode::All)) {
        std::cerr
            << "origin_selector_test: All accepted no origin membership.\n";
        return false;
    }

    if (!hbt::is_origin_eligible(
            widest_membership,
            hbt::OriginMode::All)) {
        std::cerr
            << "origin_selector_test: All rejected widest membership.\n";
        return false;
    }

    return true;
}

/**
 * @brief Verify defensive rejection of an invalid OriginMode value.
 * @return `true` when the invalid enum value throws std::invalid_argument.
 */
bool verify_invalid_mode_is_rejected() {
    try {
        static_cast<void>(hbt::is_origin_eligible(
            {true, true, true},
            static_cast<hbt::OriginMode>(999)
        ));
    } catch (const std::invalid_argument&) {
        return true;
    }

    std::cerr
        << "origin_selector_test: invalid OriginMode was accepted.\n";
    return false;
}

}  // namespace

/**
 * @brief Run the complete origin-selector unit-test collection.
 * @return `EXIT_SUCCESS` when every test passes, otherwise `EXIT_FAILURE`.
 */
int main() {
    bool success = true;

    success = verify_origin_classification() && success;
    success = verify_primordial_eligibility() && success;
    success = verify_rescattering_eligibility() && success;
    success = verify_widest_only_eligibility() && success;
    success = verify_all_uses_widest_membership() && success;
    success = verify_invalid_mode_is_rejected() && success;

    return success ? EXIT_SUCCESS : EXIT_FAILURE;
}
