/**
 * @file hbt_startup_builder_test.cpp
 * @brief Unit tests for construction of resolved HBT startup state.
 */

#include "hbt/startup/hbt_startup_builder.h"

#include "hbt/channels/primitive_channel.h"
#include "hbt/config/hbt_config.h"
#include "hbt/selection/analysis_product.h"
#include "hbt/selection/hbt_selection.h"
#include "hbt/species/species.h"

#include <cstdlib>
#include <iostream>
#include <utility>
#include <vector>

namespace {

/**
 * @brief Construct a complete HBTConfig for startup-builder tests.
 *
 * The startup builder currently derives its state from selection. The
 * remaining scientific configuration members are populated with valid,
 * representative values so that each test provides a complete HBTConfig.
 *
 * @param selection HBT selection to place in the configuration.
 *
 * @return Complete HBT configuration containing @p selection.
 */
hbt::HBTConfig make_hbt_config(
    hbt::HBTSelection selection
) {
    return hbt::HBTConfig{
        std::move(selection),
        hbt::ParticleAcceptanceConfig{
            hbt::LongitudinalVariable::Pseudorapidity,
            hbt::ParticleAcceptanceCuts{
                0.8,
                0.14,
                4.0
            },
            hbt::ParticleAcceptanceCuts{
                0.8,
                0.4,
                1.4
            },
            hbt::ParticleAcceptanceCuts{
                0.8,
                0.5,
                4.05
            },
            hbt::ParticleAcceptanceCuts{
                0.8,
                1.0,
                10000.0
            },
            hbt::ParticleAcceptanceCuts{
                0.8,
                0.3,
                10000.0
            }
        },
        hbt::PairSlicingConfig{
            hbt::PairSlicingAxisConfig{false, {}},
            hbt::PairSlicingAxisConfig{false, {}}
        },
        hbt::HBTHistogramConfig{
            {10U, 0.0, 10.0, 1.0},
            {10U, 0.0, 10.0, 1.0},
            {20U, -10.0, 10.0, 1.0}
        },
        hbt::OriginMode::All,
        hbt::FitEstimatorMode::All
    };
}

/**
 * @brief Verify startup rejects one invalid final-product selection.
 * @param selection Invalid selection to place in a complete HBT config.
 * @param label Human-readable failure label.
 * @return true when startup throws std::invalid_argument.
 */
bool verify_invalid_selection(
    hbt::HBTSelection selection,
    const char* label
) {
    try {
        static_cast<void>(hbt::build_hbt_startup_state(
            make_hbt_config(std::move(selection))
        ));
    } catch (const std::invalid_argument&) {
        return true;
    } catch (const std::exception& exception) {
        std::cerr
            << "hbt_startup_builder_test: "
            << label
            << " threw unexpected exception: "
            << exception.what()
            << ".\n";
        return false;
    }

    std::cerr
        << "hbt_startup_builder_test: "
        << label
        << " was not rejected.\n";
    return false;
}

/**
 * @brief Verify an empty final-product selection is rejected.
 * @return true when startup rejects the empty selection explicitly.
 */
bool verify_empty_selection_rejected() {
    return verify_invalid_selection(
        hbt::HBTSelection{},
        "empty selection"
    );
}

/**
 * @brief Verify one product cannot repeat the same primitive channel.
 * @return true when A+A is rejected instead of double-weighting pairs.
 */
bool verify_duplicate_channel_rejected() {
    return verify_invalid_selection(
        hbt::HBTSelection{{
            hbt::AnalysisProduct{{
                hbt::PrimitiveChannelId::PiPlusPiPlus,
                hbt::PrimitiveChannelId::PiPlusPiPlus
            }}
        }},
        "duplicate primitive channel"
    );
}

/**
 * @brief Verify one final product cannot be requested twice verbatim.
 * @return true when duplicate product state is rejected in startup.
 */
bool verify_duplicate_product_rejected() {
    const hbt::AnalysisProduct product{{
        hbt::PrimitiveChannelId::PiPlusPiPlus,
        hbt::PrimitiveChannelId::PiMinusPiMinus
    }};
    return verify_invalid_selection(
        hbt::HBTSelection{{product, product}},
        "duplicate final product"
    );
}

/**
 * @brief Verify final-product identity is independent of channel order.
 * @return true when A+B and B+A are rejected as semantic duplicates.
 */
bool verify_reordered_duplicate_product_rejected() {
    return verify_invalid_selection(
        hbt::HBTSelection{{
            hbt::AnalysisProduct{{
                hbt::PrimitiveChannelId::PiPlusPiPlus,
                hbt::PrimitiveChannelId::PiMinusPiMinus
            }},
            hbt::AnalysisProduct{{
                hbt::PrimitiveChannelId::PiMinusPiMinus,
                hbt::PrimitiveChannelId::PiPlusPiPlus
            }}
        }},
        "reordered duplicate final product"
    );
}

/**
 * @brief Verify startup-state construction for one primitive channel.
 *
 * @return true when the configured selection is preserved and its primitive
 *         channel and particle-species requirements are derived correctly.
 */
bool verify_single_channel_selection() {
    const hbt::HBTConfig config =
        make_hbt_config(
            hbt::HBTSelection{
                {
                    hbt::AnalysisProduct{
                        {
                            hbt::PrimitiveChannelId::PiPlusPiPlus
                        }
                    }
                }
            }
        );

    const hbt::HBTStartupState state =
        hbt::build_hbt_startup_state(config);

    if (state.selection.products.size() != 1U) {
        std::cerr
            << "hbt_startup_builder_test: single-channel configuration "
            << "produced the wrong number of products.\n";
        return false;
    }

    const hbt::AnalysisProduct& product =
        state.selection.products.front();

    if (product.primitive_channels.size() != 1U) {
        std::cerr
            << "hbt_startup_builder_test: preserved product contains the "
            << "wrong number of primitive channels.\n";
        return false;
    }

    if (
        product.primitive_channels.front() !=
        hbt::PrimitiveChannelId::PiPlusPiPlus
    ) {
        std::cerr
            << "hbt_startup_builder_test: configured primitive channel was "
            << "not preserved.\n";
        return false;
    }

    const std::vector<hbt::PrimitiveChannelId> expected_channels{
        hbt::PrimitiveChannelId::PiPlusPiPlus
    };

    if (state.required_primitive_channels != expected_channels) {
        std::cerr
            << "hbt_startup_builder_test: single-channel configuration "
            << "produced incorrect primitive-channel requirements.\n";
        return false;
    }

    const std::vector<hbt::SpeciesId> expected_species{
        hbt::SpeciesId::PiPlus
    };

    if (state.required_species != expected_species) {
        std::cerr
            << "hbt_startup_builder_test: single-channel configuration "
            << "produced incorrect species requirements.\n";
        return false;
    }

    return true;
}

/**
 * @brief Verify startup-state construction for overlapping analysis products.
 *
 * The configured products deliberately share KPlusProton. The startup state
 * must preserve both original products while deriving a unique, stable list
 * of primitive channels and particle species.
 *
 * @return true when selection preservation, duplicate removal, and stable
 *         requirement ordering are correct.
 */
bool verify_overlapping_products() {
    const hbt::HBTConfig config =
        make_hbt_config(
            hbt::HBTSelection{
                {
                    hbt::AnalysisProduct{
                        {
                            hbt::PrimitiveChannelId::KPlusProton,
                            hbt::PrimitiveChannelId::PiPlusProton
                        }
                    },
                    hbt::AnalysisProduct{
                        {
                            hbt::PrimitiveChannelId::KPlusProton,
                            hbt::PrimitiveChannelId::KMinusProtonBar
                        }
                    }
                }
            }
        );

    const hbt::HBTStartupState state =
        hbt::build_hbt_startup_state(config);

    if (state.selection.products.size() != 2U) {
        std::cerr
            << "hbt_startup_builder_test: overlapping configuration "
            << "produced the wrong number of preserved products.\n";
        return false;
    }

    const std::vector<hbt::PrimitiveChannelId> expected_first_product{
        hbt::PrimitiveChannelId::KPlusProton,
        hbt::PrimitiveChannelId::PiPlusProton
    };

    if (
        state.selection.products[0].primitive_channels !=
        expected_first_product
    ) {
        std::cerr
            << "hbt_startup_builder_test: first configured product was not "
            << "preserved.\n";
        return false;
    }

    const std::vector<hbt::PrimitiveChannelId> expected_second_product{
        hbt::PrimitiveChannelId::KPlusProton,
        hbt::PrimitiveChannelId::KMinusProtonBar
    };

    if (
        state.selection.products[1].primitive_channels !=
        expected_second_product
    ) {
        std::cerr
            << "hbt_startup_builder_test: second configured product was not "
            << "preserved.\n";
        return false;
    }

    const std::vector<hbt::PrimitiveChannelId> expected_channels{
        hbt::PrimitiveChannelId::KPlusProton,
        hbt::PrimitiveChannelId::PiPlusProton,
        hbt::PrimitiveChannelId::KMinusProtonBar
    };

    if (state.required_primitive_channels != expected_channels) {
        std::cerr
            << "hbt_startup_builder_test: overlapping products produced "
            << "incorrect primitive-channel requirements.\n";
        return false;
    }

    const std::vector<hbt::SpeciesId> expected_species{
        hbt::SpeciesId::KPlus,
        hbt::SpeciesId::Proton,
        hbt::SpeciesId::PiPlus,
        hbt::SpeciesId::KMinus,
        hbt::SpeciesId::ProtonBar
    };

    if (state.required_species != expected_species) {
        std::cerr
            << "hbt_startup_builder_test: overlapping products produced "
            << "incorrect species requirements.\n";
        return false;
    }

    return true;
}

}  // namespace

/**
 * @brief Run the HBT startup-state builder unit tests.
 *
 * @return EXIT_SUCCESS when every test passes, otherwise EXIT_FAILURE.
 */
int main() {
    bool success = true;

    if (!verify_empty_selection_rejected()) {
        success = false;
    }
    if (!verify_duplicate_channel_rejected()) {
        success = false;
    }
    if (!verify_duplicate_product_rejected()) {
        success = false;
    }
    if (!verify_reordered_duplicate_product_rejected()) {
        success = false;
    }

    if (!verify_single_channel_selection()) {
        success = false;
    }

    if (!verify_overlapping_products()) {
        success = false;
    }

    return success ? EXIT_SUCCESS : EXIT_FAILURE;
}
