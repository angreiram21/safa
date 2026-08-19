/**
 * @file analysis_runner.h
 * @brief Application startup, preparation, pair processing, and raw-histogram
 *        orchestration.
 */

#ifndef APP_ANALYSIS_RUNNER_H
#define APP_ANALYSIS_RUNNER_H

#include "app/analysis_progress.h"
#include "config/run_config.h"
#include "hbt/config/hbt_config.h"
#include "hbt/fits/fit_results.h"
#include "hbt/histograms/raw_histograms.h"
#include "hbt/pair/pair_processing_summary.h"
#include "hbt/reporting/rejected_particle_report.h"
#include "hbt/species/species.h"
#include "hbt/startup/hbt_startup_state.h"

#include <cstddef>
#include <cstdint>
#include <filesystem>
#include <optional>
#include <string>
#include <vector>

namespace app {

/**
 * @brief Resolved application state available before event processing.
 *
 * run_config contains the resolved global run controls.
 *
 * hbt_config contains the resolved HBT scientific configuration when HBT is
 * enabled and is std::nullopt when HBT is disabled.
 *
 * hbt_startup_state contains the HBT selection and its derived primitive
 * channel and species requirements when HBT is enabled and is std::nullopt
 * when HBT is disabled.
 *
 * The two HBT optionals are engaged together. Event readers, event paths,
 * subevent buffers, particle selection, emission-point resolution, pair
 * construction, and accumulation are intentionally absent from this startup
 * state.
 */
struct AnalysisStartupState {
    config::RunConfig run_config;  ///< Resolved global run controls.
    std::optional<hbt::HBTConfig> hbt_config;  ///< Resolved HBT config.
    /// Resolved HBT startup requirements when HBT is enabled.
    std::optional<hbt::HBTStartupState> hbt_startup_state;
};

/**
 * @brief Accepted-particle counts for one required HBT species.
 *
 * Origin counts are inclusive and follow OriginFlags. A primordial particle
 * therefore contributes to all three origin counters.
 */
struct HBTSpeciesPreparationCounts {
    hbt::SpeciesId species;                    ///< Canonical required species.
    std::size_t accepted;                      ///< Accepted particles.
    std::size_t primordial;                    ///< Primordial membership count.
    std::size_t primordial_rescattering;       ///< Middle membership count.
    std::size_t primordial_rescattering_decay; ///< Widest membership count.
};

/**
 * @brief Counts of successful emission-position resolution branches.
 */
struct HBTEmissionPointCounts {
    std::size_t sampler;      ///< Positions obtained from Sampler.
    std::size_t propagation;  ///< Positions obtained by propagation.
    std::size_t afterburner;  ///< Raw Afterburner positions retained.
};


/**
 * @brief Terminal execution state of one configured outer event.
 *
 * Failure states are reserved for explicitly classified recoverable failures.
 * Fatal conditions from earlier phases remain exceptions and never appear as
 * skipped states. This value type has no mutable shared state and requires no
 * synchronization when owned by one event result.
 */
enum class EventStatus {
    Processed,                 ///< Event contained processable input.
    SkippedDueToEventEmpty,    ///< Valid event contained only empty subevents.
    SkippedDueToEventFailure   ///< Explicit recoverable event-local failure.
};

/**
 * @brief Terminal execution state of one configured independent subevent.
 *
 * Empty is a valid zero-contribution condition. Failure is reserved for an
 * explicitly classified recoverable subevent-local failure; unrecoverable
 * input and invariant failures remain exceptions. This value type has no
 * mutable shared state and requires no synchronization when event-confined.
 */
enum class SubeventStatus {
    Processed,                    ///< Valid non-empty subevent was processed.
    SkippedDueToSubeventEmpty,    ///< Valid subevent contained zero particles.
    SkippedDueToSubeventFailure   ///< Explicit recoverable subevent failure.
};

/**
 * @brief Aggregate outer-event execution-status counts.
 *
 * The type is not internally synchronized. Parallel workers update only
 * event-local instances; run-total instances are reduced after worker join.
 */
struct HBTEventStatusCounts {
    std::uint64_t processed{0U};  ///< Successfully processed events.
    std::uint64_t skipped_due_to_event_empty{0U};  ///< Empty valid events.
    std::uint64_t skipped_due_to_event_failure{0U};  ///< Recoverable failures.
};

/**
 * @brief Aggregate subevent execution-status counts.
 *
 * The type is not internally synchronized. Parallel workers update only
 * event-local instances; run-total instances are reduced after worker join.
 */
struct HBTSubeventStatusCounts {
    std::uint64_t processed{0U};  ///< Successfully processed subevents.
    std::uint64_t skipped_due_to_subevent_empty{0U};  ///< Empty subevents.
    /// Recoverable subevent failures.
    std::uint64_t skipped_due_to_subevent_failure{0U};
};

/**
 * @brief Ordered execution status for one completed outer event.
 *
 * A worker owns the value until its event result is committed. The main
 * orchestration thread consumes it only after all workers have joined.
 */
struct HBTEventExecutionSummary {
    std::size_t outer_event_number;  ///< One-based outer-event number.
    EventStatus status{EventStatus::Processed};  ///< Terminal event status.
    /// Non-empty diagnostic for recoverable failure, otherwise empty.
    std::string diagnostic{};
};

/**
 * @brief Preparation summary for one independent subevent.
 *
 * The producing event worker owns this value. Run-total ordering is rebuilt
 * after join in canonical outer-event/subevent order.
 */
struct HBTSubeventPreparationSummary {
    std::size_t outer_event_number;  ///< One-based outer-event number.
    int subevent_id;                 ///< Subevent identifier from Afterburner.
    std::size_t accepted_particles;  ///< Accepted particles in this subevent.
    /// Terminal subevent execution status.
    SubeventStatus status{SubeventStatus::Processed};
    /// Non-empty diagnostic for recoverable failure, otherwise empty.
    std::string diagnostic{};
};

/**
 * @brief Aggregate HBT event-preparation diagnostics for one completed run.
 *
 * The summary contains aggregate counters plus complete recoverable numerical
 * rejection records. EventBuffers remain local to each subevent and are not
 * retained here.
 *
 * The type is not thread-safe for concurrent mutation. During Phase 8 each
 * worker mutates only its event-local summary. The global summary is assembled
 * by the orchestration thread after all workers have joined, preserving
 * canonical event-major order.
 */
struct HBTEventPreparationSummary {
    std::size_t outer_events_processed;  ///< Completed outer events.
    std::size_t subevents_processed;     ///< Completed independent subevents.
    std::size_t raw_particles;           ///< Raw Afterburner rows inspected.
    std::size_t unsupported_species;     ///< Unsupported charge/PDG rows.
    std::size_t unrequired_species;      ///< Supported but unrequired rows.
    /// Required-species rows rejected by configured physical acceptance.
    std::size_t particle_acceptance_rejections;
    /// Complete records for recoverable numerical particle rejections.
    hbt::RejectedParticleReport numerical_rejections;
    std::size_t origin_rejections;       ///< Rows rejected by OriginMode.
    std::size_t accepted_particles;      ///< Final particles stored.
    HBTEmissionPointCounts emission_points;  ///< Resolution-source counts.
    /// Accepted counts for each required species in startup order.
    std::vector<HBTSpeciesPreparationCounts> species;
    /// Per-subevent accepted counts in canonical event-major order.
    std::vector<HBTSubeventPreparationSummary> subevents;
    /// Aggregate event-status counts.
    HBTEventStatusCounts event_status_counts{};
    /// Aggregate subevent-status counts.
    HBTSubeventStatusCounts subevent_status_counts{};
    /// Per-event execution status in canonical outer-event order after merge.
    std::vector<HBTEventExecutionSummary> events{};
};

/**
 * @brief Result of one completed application run at the current rewrite stage.
 *
 * startup preserves the fully resolved startup state. When HBT is enabled,
 * event preparation, pair processing, complete-sample raw histograms, and
 * post-sample derived results are stored in dedicated members.
 */
struct AnalysisRunSummary {
    AnalysisStartupState startup;  ///< Fully resolved startup state.
    /// HBT preparation diagnostics, absent when HBT is disabled.
    std::optional<HBTEventPreparationSummary> hbt_event_preparation;
    /// HBT pair-processing results, absent when HBT is disabled.
    std::optional<hbt::HBTPairProcessingSummary> hbt_pair_processing;
    /// Complete-sample raw histograms, absent when HBT is disabled.
    std::optional<hbt::RawHistogramState> hbt_raw_histograms{};
    /// Post-sample histogram analysis, absent when HBT is disabled.
    std::optional<hbt::HistogramAnalysisState> hbt_histogram_analysis{};
};

/**
 * @brief Coordinate application startup and current event preparation.
 *
 * The runner is the application boundary that joins global configuration,
 * HBT startup state, event-specific input paths, raw input readers, and the
 * tested HBT event-preparation operations.
 *
 * It delegates pair processing after each completed subevent, but does not
 * implement pair formulas, numerical validation, slicing, histogram binning,
 * result writing. It owns the complete-run raw histogram lifetime, delegates
 * raw histogram accumulation and invokes post-sample analysis once after all
 * events finish.
 *
 * One run() call owns its complete worker-pool lifetime. Concurrent run() calls
 * on the same AnalysisRunner instance are outside the supported thread-safety
 * contract. The optional progress observer is borrowed synchronously and is
 * invoked only from the calling orchestration thread.
 */
class AnalysisRunner {
public:
    /**
     * @brief Construct a runner for one global run-configuration file.
     * @param run_config_path Path to the global run-configuration YAML file.
     */
    explicit AnalysisRunner(
        std::filesystem::path run_config_path
    );

    /**
     * @brief Prepare all resolved state required before event processing.
     *
     * Global configuration is loaded first. HBT scientific configuration is
     * then loaded only when enabled, after which its startup requirements are
     * derived. The resolved events_path and output_path are retained but not
     * accessed here.
     *
     * @return Resolved application startup state.
     * @throws std::runtime_error If global or enabled HBT configuration cannot
     *         be loaded or violates its YAML structure contract.
     * @throws std::invalid_argument If enabled HBT scientific configuration is
     *         invalid or an internal enabled-HBT invariant is violated.
     */
    [[nodiscard]] AnalysisStartupState prepare_startup() const;

    /**
     * @brief Execute startup, event preparation, and current pair processing.
     *
     * When HBT is disabled, the method returns after startup without accessing
     * events_path. When HBT is enabled, outer-event directories are addressed
     * as one-based numeric children of events_path. Each outer event uses a
     * fresh SamplerReader and one streaming AfterburnerReader from the same
     * outer-event directory.
     *
     * A fresh EventBuffers instance is created for every Afterburner subevent.
     * After the subevent input is structurally complete, required primitive
     * pairs are processed from that buffer before it is destroyed. The HBT pair
     * module calculates and validates kT/mT once per formed pair and, when the
     * configured path requires frames, validates the calculated frame
     * observables before committing valid counts. With one worker the runner
     * owns the serial complete-sample raw state. With multiple workers each
     * worker owns one private state with the identical layout; after join these
     * states are reduced bin-for-bin into one complete-sample raw state before
     * post-sample analysis.
     *
     * @param progress Optional non-owning orchestration progress observer. The
     *        pointer is never retained and workers never invoke it directly.
     * @return Startup, event-preparation, pair-processing, raw-histogram, and
     *         post-sample derived results.
     * The configured threads control applies only to outer-event processing.
     * A value of zero resolves hardware concurrency, one preserves the serial
     * path, and larger values enable a dynamically scheduled worker pool.
     * Post-sample analysis and output remain serial.
     *
     * Recoverable event-local input failures and mandatory Sampler misses are
     * returned through the explicit event/subevent failure statuses instead
     * of escaping this boundary. Unrecoverable run-global failures still
     * propagate as exceptions.
     *
     * @throws std::overflow_error If pair or raw-histogram accounting would
     *         overflow.
     * @throws std::logic_error If internal pair or raw-histogram invariants
     *         fail.
     */
    [[nodiscard]] AnalysisRunSummary run(
        AnalysisProgressObserver* progress = nullptr
    ) const;

private:
    /// Global run-configuration file used by prepare_startup() and run().
    std::filesystem::path run_config_path_;
};

}  // namespace app

#endif  // APP_ANALYSIS_RUNNER_H
