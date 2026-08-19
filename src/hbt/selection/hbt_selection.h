/**
 * @file hbt_selection.h
 * @brief Definition of the requested HBT analysis selection.
 *
 * An HBTSelection contains the final analysis products requested for one HBT
 * analysis configuration.
 *
 * This type represents only the resulting selection. Configuration parsing,
 * validation, primitive-channel dependency resolution, particle filtering,
 * pair construction, observable calculation, histogram allocation, fitting,
 * and output generation are outside its responsibility.
 */

#ifndef HBT_SELECTION_HBT_SELECTION_H
#define HBT_SELECTION_HBT_SELECTION_H

#include "hbt/selection/analysis_product.h"

#include <vector>

namespace hbt {

    /**
     * @brief Collection of requested final HBT analysis products.
     *
     * Each AnalysisProduct describes one final result and the primitive
     * channels contributing to it.
     *
     * Primitive channels may occur in more than one product. Such overlap
     * describes shared product composition and must not imply repeated
     * calculation of the corresponding particle pairs.
     *
     * Validation of the selection, including rejection of invalid or empty
     * configurations, belongs to the selection-construction layer.
     */
    struct HBTSelection {
        /**
         * @brief Final HBT analysis products requested by the selection.
         *
         * Each entry represents one independently requested final product.
         */
        std::vector<AnalysisProduct> products;
    };

}  // namespace hbt

#endif  // HBT_SELECTION_HBT_SELECTION_H
