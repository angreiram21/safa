/**
 * @file analysis_output_writer.h
 * @brief Serialization boundary for analysis outputs.
 */

#ifndef SMASH_AFTERBURNER_ANALYSIS_OUTPUT_ANALYSIS_OUTPUT_WRITER_H
#define SMASH_AFTERBURNER_ANALYSIS_OUTPUT_ANALYSIS_OUTPUT_WRITER_H

#include "app/analysis_runner.h"
#include "hbt/reporting/rejected_pair_report.h"
#include "hbt/reporting/rejected_particle_report.h"

#include <iosfwd>

namespace output {

/**
 * @brief Serialize a rejected-particle report to an output stream.
 *
 * The writer does not modify the report or perform scientific calculations.
 *
 * @param report Complete in-memory numerical rejection report.
 * @param output Destination stream supplied by the application boundary.
 */
void write_rejected_particle_report(
    const hbt::RejectedParticleReport& report,
    std::ostream& output
);

/**
 * @brief Serialize a rejected-pair report to an output stream.
 *
 * The writer does not modify the report or perform scientific calculations.
 *
 * @param report Complete in-memory numerical pair-rejection report.
 * @param output Destination stream supplied by the application boundary.
 */
void write_rejected_pair_report(
    const hbt::RejectedPairReport& report,
    std::ostream& output
);

/**
 * @brief Serialize all currently available analysis-run outputs.
 *
 * This function is the analysis-output boundary. Scientific modules
 * provide data structures and never choose output formats or write files.
 *
 * @param result Completed analysis-run result to serialize.
 * @param output Destination stream supplied by the application boundary.
 * @throws std::logic_error If a present pair summary or raw histogram state
 *         is inconsistent with its resolved startup configuration.
 */
void write_analysis_output(
    const app::AnalysisRunSummary& result,
    std::ostream& output
);

/**
 * @brief Serialize filesystem production output for one completed run.
 *
 * The configured output root is taken from result.startup.run_config. This
 * boundary does not perform scientific calculations or select an alternative
 * path. Module-specific filesystem writers remain responsible for their own
 * canonical hierarchy beneath that root.
 *
 * @param result Completed analysis-run result to serialize. The function
 *        borrows this state synchronously and retains no references after
 *        return.
 * @pre result.startup.run_config.output_path is the explicitly resolved
 *      production-output root supplied by the run configuration.
 * @throws std::logic_error If enabled module state is incomplete or internally
 *         inconsistent.
 * @throws std::invalid_argument If the configured output root is empty.
 * @throws std::runtime_error If production output cannot be written according
 *         to the filesystem-output contract.
 * @throws std::filesystem::filesystem_error If filesystem operations fail.
 *
 * Ownership of the result state and configured path remains with the caller.
 * This function performs no physical recalculation and does not alter the
 * canonical identity or hierarchy of module output.
 */
void write_production_output(
    const app::AnalysisRunSummary& result
);

}  // namespace output

#endif
