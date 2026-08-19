/**
 * @file analysis_progress.h
 * @brief Application-level progress reporting interface.
 */

#ifndef APP_ANALYSIS_PROGRESS_H
#define APP_ANALYSIS_PROGRESS_H

#include <cstddef>

namespace app {

/**
 * @brief Observer for non-scientific analysis progress events.
 *
 * Progress notifications describe completed orchestration work only. They do
 * not participate in particle selection, pair processing, histogramming, or
 * fitting decisions.
 */
class AnalysisProgressObserver {
public:
    /** @brief Destroy the observer through its interface. */
    virtual ~AnalysisProgressObserver() = default;

    /**
     * @brief Start progress reporting for an enabled-HBT run.
     * @param total_events Configured number of outer events.
     * @param total_subevents Total configured subevents across all events.
     */
    virtual void begin(
        std::size_t total_events,
        std::size_t total_subevents
    ) = 0;

    /**
     * @brief Report one fully completed subevent.
     * @param completed_subevents Completed subevents across the whole run.
     * @param outer_event_number One-based canonical event-major progress
     *        position. In serial execution this is the completed event; in
     *        parallel execution it does not expose scheduler completion order.
     * @param subevent_number One-based canonical subevent position in the
     *        reported outer-event progress position.
     */
    virtual void subevent_completed(
        std::size_t completed_subevents,
        std::size_t outer_event_number,
        std::size_t subevent_number
    ) = 0;

    /** @brief Report entry into post-sample statistical analysis. */
    virtual void begin_postprocessing() = 0;

    /** @brief Report successful completion of in-memory analysis. */
    virtual void analysis_complete() = 0;

    /** @brief Report entry into production-output serialization. */
    virtual void begin_output() = 0;

    /** @brief Report successful completion of the complete run. */
    virtual void finish() = 0;

    /**
     * @brief End any active display before an exception is printed.
     *
     * This operation must not throw.
     */
    virtual void fail() noexcept = 0;
};

}  // namespace app

#endif  // APP_ANALYSIS_PROGRESS_H
