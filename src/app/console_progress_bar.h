/**
 * @file console_progress_bar.h
 * @brief Adaptive terminal progress display for application orchestration.
 */

#ifndef APP_CONSOLE_PROGRESS_BAR_H
#define APP_CONSOLE_PROGRESS_BAR_H

#include "app/analysis_progress.h"

#include <cstddef>
#include <iosfwd>
#include <optional>
#include <string>

namespace app {

/**
 * @brief Render global HBT progress to a terminal or plain-text log.
 *
 * The production constructor writes to stderr, detects whether stderr is an
 * interactive terminal, and queries the terminal width for every redraw so a
 * resize is reflected without restarting the analysis. Interactive output is
 * rendered in ANSI bright green. Redirected output is plain text and is
 * throttled to coarse percentage milestones.
 */
class ConsoleProgressBar final : public AnalysisProgressObserver {
public:
    /**
     * @brief Construct the production stderr progress display.
     */
    ConsoleProgressBar();

    /**
     * @brief Construct a deterministic display for tests or custom sinks.
     * @param stream Destination stream.
     * @param interactive Whether terminal-style redraws and color are enabled.
     * @param terminal_width Fixed visible width used for formatting.
     * @throws std::invalid_argument If terminal_width is zero.
     */
    ConsoleProgressBar(
        std::ostream& stream,
        bool interactive,
        std::size_t terminal_width
    );

    /** @brief End an active interactive block without throwing. */
    ~ConsoleProgressBar() override;

    /** @copydoc AnalysisProgressObserver::begin */
    void begin(
        std::size_t total_events,
        std::size_t total_subevents
    ) override;

    /** @copydoc AnalysisProgressObserver::subevent_completed */
    void subevent_completed(
        std::size_t completed_subevents,
        std::size_t outer_event_number,
        std::size_t subevent_number
    ) override;

    /** @copydoc AnalysisProgressObserver::begin_postprocessing */
    void begin_postprocessing() override;

    /** @copydoc AnalysisProgressObserver::analysis_complete */
    void analysis_complete() override;

    /** @copydoc AnalysisProgressObserver::begin_output */
    void begin_output() override;

    /** @copydoc AnalysisProgressObserver::finish */
    void finish() override;

    /** @copydoc AnalysisProgressObserver::fail */
    void fail() noexcept override;

private:
    /** @brief Return the current or fixed visible terminal width. */
    [[nodiscard]] std::size_t terminal_width() const noexcept;

    /** @brief Erase the currently displayed interactive block. */
    void erase_interactive_block();

    /** @brief Render the current progress state for an interactive terminal. */
    void render_interactive(const std::string& detail);

    /** @brief Emit a throttled plain-text progress record. */
    void render_noninteractive_progress();

    /** @brief Emit one plain-text stage record. */
    void render_noninteractive_stage(const std::string& stage);

    /** @brief End the active interactive block with a newline. */
    void end_interactive_block() noexcept;

    std::ostream* stream_;  ///< Output stream; never null after construction.
    bool interactive_;      ///< Whether redraw control sequences are enabled.
    bool color_enabled_;    ///< Whether ANSI bright-green output is enabled.
    /// Fixed width for deterministic sinks; absent for live terminal queries.
    std::optional<std::size_t> fixed_terminal_width_;
    std::size_t total_events_;       ///< Configured outer-event count.
    std::size_t total_subevents_;    ///< Configured global subevent count.
    std::size_t completed_subevents_;  ///< Completed global subevent count.
    std::size_t outer_event_number_; ///< Current one-based outer event.
    std::size_t subevent_number_;    ///< Current one-based local subevent.
    std::size_t rendered_lines_;     ///< Lines in the active terminal block.
    unsigned int logged_decile_;     ///< Last plain-text 10% milestone.
    bool begun_;                     ///< Whether enabled-HBT progress began.
    bool block_active_;              ///< Whether a terminal block is active.
};

}  // namespace app

#endif  // APP_CONSOLE_PROGRESS_BAR_H
