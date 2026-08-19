/**
 * @file required_channels_test.cpp
 * @brief Unit tests for primitive-channel dependency derivation.
 *
 * This test verifies that required_primitive_channels() derives the unique
 * primitive HBT channels required by an HBTSelection while preserving their
 * order of first occurrence.
 */

#include "hbt/selection/required_channels.h"

#include <cstddef>
#include <cstdlib>
#include <initializer_list>
#include <iostream>
#include <vector>

namespace {

    /**
     * @brief Construct an analysis product for a unit-test case.
     *
     * @param channels Primitive channels contributing to the product.
     *
     * @return Analysis product containing the supplied primitive channels.
     */
    hbt::AnalysisProduct make_product(
        std::initializer_list<hbt::PrimitiveChannelId> channels
    ) {
        return hbt::AnalysisProduct{
            std::vector<hbt::PrimitiveChannelId>(channels)
        };
    }

    /**
     * @brief Verify the primitive channels derived for one test selection.
     *
     * The result must contain exactly the expected channels in exactly the
     * expected order.
     *
     * @param test_name Human-readable name of the test case.
     * @param selection HBT selection to inspect.
     * @param expected Expected unique primitive channels.
     *
     * @return true when the derived result matches the expectation, otherwise
     *         false.
     */
    bool verify_result(
        const char* test_name,
        const hbt::HBTSelection& selection,
        const std::vector<hbt::PrimitiveChannelId>& expected
    ) {
        const std::vector<hbt::PrimitiveChannelId> actual =
            hbt::required_primitive_channels(selection);

        if (actual.size() != expected.size()) {
            std::cerr
                << "required_channels_test: "
                << test_name
                << " returned "
                << actual.size()
                << " channels; expected "
                << expected.size()
                << ".\n";
            return false;
        }

        bool success = true;

        for (std::size_t index = 0; index < expected.size(); ++index) {
            if (actual[index] != expected[index]) {
                std::cerr
                    << "required_channels_test: "
                    << test_name
                    << " has an incorrect channel at index "
                    << index
                    << ".\n";
                success = false;
            }
        }

        return success;
    }

    /**
     * @brief Verify that an empty selection requires no primitive channels.
     *
     * @return true when the derived result is empty, otherwise false.
     */
    bool verify_empty_selection() {
        const hbt::HBTSelection selection{};

        const std::vector<hbt::PrimitiveChannelId> expected{};

        return verify_result(
            "empty selection",
            selection,
            expected
        );
    }

    /**
     * @brief Verify that products containing no channels require no channels.
     *
     * @return true when the derived result is empty, otherwise false.
     */
    bool verify_products_without_channels() {
        const hbt::HBTSelection selection{{
            make_product({}),
            make_product({}),
            make_product({})
        }};

        const std::vector<hbt::PrimitiveChannelId> expected{};

        return verify_result(
            "products without channels",
            selection,
            expected
        );
    }

    /**
     * @brief Verify preservation of channel order for one product.
     *
     * @return true when all channels are returned in their stored order,
     *         otherwise false.
     */
    bool verify_single_product() {
        const hbt::HBTSelection selection{{
            make_product({
                hbt::PrimitiveChannelId::PiPlusPiPlus,
                hbt::PrimitiveChannelId::KMinusProton,
                hbt::PrimitiveChannelId::LambdaLambda
            })
        }};

        const std::vector<hbt::PrimitiveChannelId> expected{
            hbt::PrimitiveChannelId::PiPlusPiPlus,
            hbt::PrimitiveChannelId::KMinusProton,
            hbt::PrimitiveChannelId::LambdaLambda
        };

        return verify_result(
            "single product",
            selection,
            expected
        );
    }

    /**
     * @brief Verify concatenation of multiple products without overlap.
     *
     * @return true when all channels are returned in product and channel order,
     *         otherwise false.
     */
    bool verify_multiple_products_without_overlap() {
        const hbt::HBTSelection selection{{
            make_product({
                hbt::PrimitiveChannelId::PiPlusPiPlus,
                hbt::PrimitiveChannelId::KMinusProton
            }),
            make_product({
                hbt::PrimitiveChannelId::LambdaLambda,
                hbt::PrimitiveChannelId::PiZeroSigmaZero
            })
        }};

        const std::vector<hbt::PrimitiveChannelId> expected{
            hbt::PrimitiveChannelId::PiPlusPiPlus,
            hbt::PrimitiveChannelId::KMinusProton,
            hbt::PrimitiveChannelId::LambdaLambda,
            hbt::PrimitiveChannelId::PiZeroSigmaZero
        };

        return verify_result(
            "multiple products without overlap",
            selection,
            expected
        );
    }

    /**
     * @brief Verify removal of repeated channels within one product.
     *
     * The first occurrence of every primitive channel must be retained and
     * later occurrences must not alter the resulting order.
     *
     * @return true when duplicates are removed correctly, otherwise false.
     */
    bool verify_duplicates_within_product() {
        const hbt::HBTSelection selection{{
            make_product({
                hbt::PrimitiveChannelId::PiPlusPiPlus,
                hbt::PrimitiveChannelId::PiMinusPiMinus,
                hbt::PrimitiveChannelId::PiPlusPiPlus,
                hbt::PrimitiveChannelId::KPlusKPlus,
                hbt::PrimitiveChannelId::PiMinusPiMinus,
                hbt::PrimitiveChannelId::PiPlusPiPlus
            })
        }};

        const std::vector<hbt::PrimitiveChannelId> expected{
            hbt::PrimitiveChannelId::PiPlusPiPlus,
            hbt::PrimitiveChannelId::PiMinusPiMinus,
            hbt::PrimitiveChannelId::KPlusKPlus
        };

        return verify_result(
            "duplicates within product",
            selection,
            expected
        );
    }

    /**
     * @brief Verify deduplication across products and first-occurrence order.
     *
     * Channels appearing in several products must occur only once in the
     * result. Their position must be determined by their first occurrence
     * while traversing products and their stored channels.
     *
     * @return true when overlap and ordering are handled correctly, otherwise
     *         false.
     */
    bool verify_overlap_between_products() {
        const hbt::HBTSelection selection{{
            make_product({
                hbt::PrimitiveChannelId::KMinusProton,
                hbt::PrimitiveChannelId::PiPlusProton
            }),
            make_product({
                hbt::PrimitiveChannelId::PiPlusProton,
                hbt::PrimitiveChannelId::LambdaLambda,
                hbt::PrimitiveChannelId::KMinusProton,
                hbt::PrimitiveChannelId::PiPlusSigmaBarMinus
            }),
            make_product({
                hbt::PrimitiveChannelId::LambdaLambda,
                hbt::PrimitiveChannelId::PiMinusPiMinus
            })
        }};

        const std::vector<hbt::PrimitiveChannelId> expected{
            hbt::PrimitiveChannelId::KMinusProton,
            hbt::PrimitiveChannelId::PiPlusProton,
            hbt::PrimitiveChannelId::LambdaLambda,
            hbt::PrimitiveChannelId::PiPlusSigmaBarMinus,
            hbt::PrimitiveChannelId::PiMinusPiMinus
        };

        return verify_result(
            "overlap between products",
            selection,
            expected
        );
    }

}  // namespace

/**
 * @brief Run the required primitive-channel unit tests.
 *
 * @return EXIT_SUCCESS when every test passes, otherwise EXIT_FAILURE.
 */
int main() {
    bool success = true;

    if (!verify_empty_selection()) {
        success = false;
    }

    if (!verify_products_without_channels()) {
        success = false;
    }

    if (!verify_single_product()) {
        success = false;
    }

    if (!verify_multiple_products_without_overlap()) {
        success = false;
    }

    if (!verify_duplicates_within_product()) {
        success = false;
    }

    if (!verify_overlap_between_products()) {
        success = false;
    }

    return success ? EXIT_SUCCESS : EXIT_FAILURE;
}
