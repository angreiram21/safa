/**
 * @file pair_origin_route_count_summary.h
 * @brief Pair counts routed to requested nested origin selections.
 */

#ifndef SMASH_AFTERBURNER_ANALYSIS_HBT_PAIR_PAIR_ORIGIN_ROUTE_COUNT_SUMMARY_H
#define SMASH_AFTERBURNER_ANALYSIS_HBT_PAIR_PAIR_ORIGIN_ROUTE_COUNT_SUMMARY_H

#include "hbt/config/hbt_config.h"
#include "hbt/pair/pair_count_summary.h"

namespace hbt {

/**
 * @brief Pair counts routed to requested origin selections.
 *
 * The three count summaries are nested output routes, not disjoint physical
 * categories. Under OriginMode::All, one valid physical pair may increment
 * more than one route. origin_mode records the routing policy that produced
 * the counts so local and run-total summaries cannot be mixed silently.
 */
struct PairOriginRouteCountSummary {
    OriginMode origin_mode;  ///< Configured routing policy.
    PairCountSummary routed_P;  ///< Routed primordial pairs by channel.
    PairCountSummary routed_PR;  ///< Routed primordial+rescattering pairs.
    PairCountSummary routed_PRD;  ///< Routed widest-origin pairs by channel.
};

}  // namespace hbt

#endif  // SMASH_AFTERBURNER_ANALYSIS_HBT_PAIR_PAIR_ORIGIN_ROUTE_COUNT_SUMMARY_H
