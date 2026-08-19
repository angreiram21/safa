/**
 * @file console_progress_bar.cpp
 * @brief Adaptive terminal progress display implementation.
 */

#include "app/console_progress_bar.h"

#include <algorithm>
#include <cmath>
#include <cstdlib>
#include <iomanip>
#include <iostream>
#include <sstream>
#include <stdexcept>
#include <string>
#include <vector>

#include <sys/ioctl.h>
#include <unistd.h>

namespace app {
namespace {

constexpr std::size_t kFallbackTerminalWidth = 80U;
constexpr std::size_t kMaximumBarWidth = 60U;
constexpr std::size_t kMinimumBarLineWidth = 14U;
constexpr const char* kBrightGreen = "\033[92m";
constexpr const char* kResetColor = "\033[0m";
constexpr const char* kClearLine = "\r\033[2K";
constexpr const char* kCursorUp = "\033[1A";

/**
 * @brief Test whether stderr supports interactive terminal rendering.
 * @return true when stderr is attached to a terminal.
 */
bool stderr_is_interactive() noexcept {
    return ::isatty(STDERR_FILENO) == 1;
}

/**
 * @brief Test whether ANSI color should be used for an interactive terminal.
 * @return true unless TERM explicitly declares a dumb terminal.
 */
bool stderr_supports_color() noexcept {
    const char* term = std::getenv("TERM");
    return stderr_is_interactive() &&
           (term == nullptr || std::string(term) != "dumb");
}

/**
 * @brief Query the current stderr terminal width.
 * @return Positive terminal width, or the fallback width when unavailable.
 */
std::size_t current_stderr_width() noexcept {
    struct winsize size {};

    if (::ioctl(STDERR_FILENO, TIOCGWINSZ, &size) == 0 &&
        size.ws_col > 0U) {
        return static_cast<std::size_t>(size.ws_col);
    }

    return kFallbackTerminalWidth;
}

/**
 * @brief Format one percentage with one decimal place.
 * @param completed Completed work units.
 * @param total Total work units.
 * @return Percentage text including the percent sign.
 */
std::string percentage_text(
    std::size_t completed,
    std::size_t total
) {
    const double ratio = total == 0U
        ? 0.0
        : static_cast<double>(completed) /
              static_cast<double>(total);

    std::ostringstream stream;
    stream
        << std::fixed
        << std::setprecision(1)
        << std::setw(5)
        << (100.0 * ratio)
        << '%';
    return stream.str();
}

/**
 * @brief Build the adaptive first progress line.
 * @param completed Completed subevents.
 * @param total Total subevents.
 * @param width Visible terminal width.
 * @return A line whose visible length does not exceed width.
 */
std::string progress_line(
    std::size_t completed,
    std::size_t total,
    std::size_t width
) {
    const std::string percent = percentage_text(completed, total);

    if (width < kMinimumBarLineWidth) {
        const std::string compact = "HBT " + percent;
        return compact.substr(0U, width);
    }

    const std::size_t available_bar_width = width - 13U;
    const std::size_t bar_width = std::min(
        kMaximumBarWidth,
        available_bar_width
    );
    const double ratio = total == 0U
        ? 0.0
        : static_cast<double>(completed) /
              static_cast<double>(total);
    const std::size_t filled = std::min(
        bar_width,
        static_cast<std::size_t>(
            std::floor(ratio * static_cast<double>(bar_width))
        )
    );

    std::string bar(bar_width, '-');

    if (completed >= total && total > 0U) {
        std::fill(bar.begin(), bar.end(), '=');
    } else if (filled < bar_width) {
        std::fill(bar.begin(), bar.begin() + filled, '=');
        bar[filled] = '>';
    }

    return "HBT [" + bar + "] " + percent;
}

/**
 * @brief Choose a subevent-count detail line that fits the terminal.
 * @param completed Completed subevents.
 * @param total Total subevents.
 * @param width Visible terminal width.
 * @return Fitting detail line, or an empty string when no form fits.
 */
std::string subevent_detail(
    std::size_t completed,
    std::size_t total,
    std::size_t width
) {
    const std::string spaced =
        "    " + std::to_string(completed) + " / " +
        std::to_string(total) + " subevents";

    if (spaced.size() <= width) {
        return spaced;
    }

    const std::string compact =
        "    " + std::to_string(completed) + "/" +
        std::to_string(total) + " subevents";

    if (compact.size() <= width) {
        return compact;
    }

    const std::string counts =
        "    " + std::to_string(completed) + "/" +
        std::to_string(total);

    return counts.size() <= width ? counts : std::string{};
}

/**
 * @brief Choose an event detail line that fits the terminal.
 * @param current Current one-based outer-event number.
 * @param total Total outer-event count.
 * @param width Visible terminal width.
 * @return Fitting event line, or an empty string when no form fits.
 */
std::string event_detail(
    std::size_t current,
    std::size_t total,
    std::size_t width
) {
    const std::string spaced =
        "    event " + std::to_string(current) + " / " +
        std::to_string(total);

    if (spaced.size() <= width) {
        return spaced;
    }

    const std::string compact =
        "    event " + std::to_string(current) + "/" +
        std::to_string(total);

    return compact.size() <= width ? compact : std::string{};
}

/**
 * @brief Choose a stage detail line that fits the terminal.
 * @param stage Stage description.
 * @param width Visible terminal width.
 * @return Fitting stage line, or an empty string when it cannot fit.
 */
std::string stage_detail(
    const std::string& stage,
    std::size_t width
) {
    const std::string line = "    " + stage;
    return line.size() <= width ? line : std::string{};
}

}  // namespace

ConsoleProgressBar::ConsoleProgressBar()
    : stream_(&std::cerr),
      interactive_(stderr_is_interactive()),
      color_enabled_(stderr_supports_color()),
      fixed_terminal_width_(std::nullopt),
      total_events_(0U),
      total_subevents_(0U),
      completed_subevents_(0U),
      outer_event_number_(0U),
      subevent_number_(0U),
      rendered_lines_(0U),
      logged_decile_(0U),
      begun_(false),
      block_active_(false) {
}

ConsoleProgressBar::ConsoleProgressBar(
    std::ostream& stream,
    bool interactive,
    std::size_t terminal_width
) : stream_(&stream),
    interactive_(interactive),
    color_enabled_(interactive),
    fixed_terminal_width_(terminal_width),
    total_events_(0U),
    total_subevents_(0U),
    completed_subevents_(0U),
    outer_event_number_(0U),
    subevent_number_(0U),
    rendered_lines_(0U),
    logged_decile_(0U),
    begun_(false),
    block_active_(false) {
    if (terminal_width == 0U) {
        throw std::invalid_argument(
            "ConsoleProgressBar: terminal width must be positive"
        );
    }
}

ConsoleProgressBar::~ConsoleProgressBar() {
    end_interactive_block();
}

void ConsoleProgressBar::begin(
    std::size_t total_events,
    std::size_t total_subevents
) {
    total_events_ = total_events;
    total_subevents_ = total_subevents;
    completed_subevents_ = 0U;
    outer_event_number_ = total_events > 0U ? 1U : 0U;
    subevent_number_ = 0U;
    logged_decile_ = 0U;
    begun_ = true;

    if (interactive_) {
        render_interactive(
            event_detail(
                outer_event_number_,
                total_events_,
                terminal_width()
            )
        );
    } else {
        *stream_ << "HBT progress: 0% (0/"
                 << total_subevents_
                 << " subevents)\n";
    }
}

void ConsoleProgressBar::subevent_completed(
    std::size_t completed_subevents,
    std::size_t outer_event_number,
    std::size_t subevent_number
) {
    if (!begun_) {
        return;
    }

    completed_subevents_ = completed_subevents;
    outer_event_number_ = outer_event_number;
    subevent_number_ = subevent_number;

    if (interactive_) {
        render_interactive(
            event_detail(
                outer_event_number_,
                total_events_,
                terminal_width()
            )
        );
    } else {
        render_noninteractive_progress();
    }
}

void ConsoleProgressBar::begin_postprocessing() {
    if (!begun_) {
        return;
    }

    completed_subevents_ = total_subevents_;

    if (interactive_) {
        render_interactive(
            stage_detail(
                "post-sample statistical analysis",
                terminal_width()
            )
        );
    } else {
        render_noninteractive_progress();
        render_noninteractive_stage(
            "post-sample statistical analysis"
        );
    }
}

void ConsoleProgressBar::analysis_complete() {
    if (!begun_) {
        return;
    }

    if (interactive_) {
        render_interactive(
            stage_detail("analysis complete", terminal_width())
        );
        end_interactive_block();
    } else {
        render_noninteractive_stage("analysis complete");
    }
}

void ConsoleProgressBar::begin_output() {
    if (!begun_) {
        return;
    }

    if (interactive_) {
        render_interactive(
            stage_detail(
                "writing production output",
                terminal_width()
            )
        );
    } else {
        render_noninteractive_stage("writing production output");
    }
}

void ConsoleProgressBar::finish() {
    if (!begun_) {
        return;
    }

    if (interactive_) {
        render_interactive(stage_detail("complete", terminal_width()));
        end_interactive_block();
    } else {
        render_noninteractive_stage("complete");
    }

    begun_ = false;
}

void ConsoleProgressBar::fail() noexcept {
    end_interactive_block();
    begun_ = false;
}

std::size_t ConsoleProgressBar::terminal_width() const noexcept {
    if (fixed_terminal_width_.has_value()) {
        return fixed_terminal_width_.value();
    }

    return current_stderr_width();
}

void ConsoleProgressBar::erase_interactive_block() {
    if (!block_active_) {
        return;
    }

    for (std::size_t line = 0U; line < rendered_lines_; ++line) {
        *stream_ << kClearLine;

        if (line + 1U < rendered_lines_) {
            *stream_ << kCursorUp;
        }
    }
}

void ConsoleProgressBar::render_interactive(const std::string& detail) {
    erase_interactive_block();

    const std::size_t width = terminal_width();
    std::vector<std::string> lines;
    lines.push_back(
        progress_line(
            completed_subevents_,
            total_subevents_,
            width
        )
    );

    const std::string subevents = subevent_detail(
        completed_subevents_,
        total_subevents_,
        width
    );

    if (!subevents.empty()) {
        lines.push_back(subevents);
    }

    if (!detail.empty()) {
        lines.push_back(detail);
    }

    if (color_enabled_) {
        *stream_ << kBrightGreen;
    }

    for (std::size_t index = 0U; index < lines.size(); ++index) {
        if (index > 0U) {
            *stream_ << '\n';
        }
        *stream_ << lines[index];
    }

    if (color_enabled_) {
        *stream_ << kResetColor;
    }

    stream_->flush();
    rendered_lines_ = lines.size();
    block_active_ = true;
}

void ConsoleProgressBar::render_noninteractive_progress() {
    if (total_subevents_ == 0U) {
        return;
    }

    const long double scaled =
        10.0L * static_cast<long double>(completed_subevents_) /
        static_cast<long double>(total_subevents_);
    const unsigned int decile =
        static_cast<unsigned int>(scaled);

    if (decile <= logged_decile_) {
        return;
    }

    logged_decile_ = decile;
    *stream_ << "HBT progress: "
             << (10U * decile)
             << "% ("
             << completed_subevents_
             << '/'
             << total_subevents_
             << " subevents)\n";
}

void ConsoleProgressBar::render_noninteractive_stage(
    const std::string& stage
) {
    *stream_ << "HBT stage: " << stage << '\n';
}

void ConsoleProgressBar::end_interactive_block() noexcept {
    if (!block_active_) {
        return;
    }

    try {
        *stream_ << '\n';
        stream_->flush();
    } catch (...) {
    }

    block_active_ = false;
    rendered_lines_ = 0U;
}

}  // namespace app
