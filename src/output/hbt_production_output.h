/**
 * @file hbt_production_output.h
 * @brief Filesystem boundary for post-sample HBT production output.
 */

#ifndef OUTPUT_HBT_PRODUCTION_OUTPUT_H
#define OUTPUT_HBT_PRODUCTION_OUTPUT_H

#include "app/analysis_runner.h"

#include <filesystem>

namespace output {

/**
 * @brief Serialize canonical post-sample HBT production files under one root.
 * @param result Completed analysis result containing raw and derived HBT state.
 * @param output_root Explicit caller-owned production output root.
 * @throws std::logic_error If enabled HBT result state is incomplete or its
 *         raw/derived layout is inconsistent with startup configuration.
 * @throws std::invalid_argument If @p output_root is empty.
 * @throws std::runtime_error If the root already contains entries or a file
 *         cannot be opened or written.
 * @throws std::filesystem::filesystem_error If filesystem operations fail.
 *
 * The caller owns @p output_root and all state in @p result. This synchronous
 * function retains no references after return. The root must be absent or
 * empty so stale files from an earlier run can never masquerade as current
 * output. Scientific calculations are not performed here: selected regions,
 * normalized bins, fitted densities, parameters, and statistics are consumed
 * exactly as already present in the derived state. One root-level
 * `product_catalog.csv` records configured product expressions together with
 * their resolved canonical channels and species for run traceability.
 */
void write_hbt_production_output(
    const app::AnalysisRunSummary& result,
    const std::filesystem::path& output_root
);

}  // namespace output

#endif  // OUTPUT_HBT_PRODUCTION_OUTPUT_H
