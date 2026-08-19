/**
 * @file afterburner_reader.cpp
 * @brief Implementation of SMASH Afterburner OSCAR input streaming.
 */

#include "input/afterburner_reader.h"

#include <charconv>
#include <sstream>
#include <stdexcept>
#include <string>
#include <string_view>
#include <system_error>

namespace input {
namespace {

constexpr std::string_view expected_column_header =
    "#!OSCAR2013Extended particle_lists "
    "t x y z mass p0 px py pz pdg ID charge ncoll form_time xsecfac "
    "proc_id_origin proc_type_origin time_last_coll pdg_mother1 pdg_mother2 "
    "baryon_number strangeness";

constexpr std::string_view expected_unit_header =
    "# Units: fm fm fm fm GeV GeV GeV GeV GeV none none e none fm none "
    "none none fm none none none none";

/**
 * @brief Reads one mandatory line from the Afterburner stream.
 * @param input Input stream positioned before the required line.
 * @param description Human-readable description used in failure diagnostics.
 * @return The complete line without the trailing newline character.
 * @throws std::runtime_error If the required line cannot be read.
 */
std::string read_required_line(
    std::ifstream& input,
    std::string_view description) {
    std::string line;

    if (!std::getline(input, line)) {
        throw std::runtime_error(
            "Missing or unreadable Afterburner " + std::string(description) +
            ".");
    }

    return line;
}

/**
 * @brief Normalizes whitespace between tokens in one input line.
 * @param line Raw input line.
 * @return Tokens joined by one ASCII space without surrounding whitespace.
 */
std::string normalize_whitespace(const std::string& line) {
    std::istringstream stream(line);
    std::ostringstream normalized;
    std::string token;
    bool first = true;

    while (stream >> token) {
        if (!first) {
            normalized << ' ';
        }

        normalized << token;
        first = false;
    }

    return normalized.str();
}

/**
 * @brief Tests whether a line declares a syntactically valid SMASH version.
 * @param line Candidate version line from the OSCAR file header.
 * @return `true` only for `# SMASH-<version>` with no trailing tokens.
 */
bool is_smash_version_line(const std::string& line) {
    std::istringstream stream(line);
    std::string marker;
    std::string version;
    std::string extra;

    if (!(stream >> marker >> version) || (stream >> extra)) {
        return false;
    }

    return marker == "#" && version.size() > 6 &&
           version.rfind("SMASH-", 0) == 0;
}

/**
 * @brief Validates and consumes the three-line OSCAR file header.
 * @param input Newly opened Afterburner input stream.
 * @throws std::runtime_error If the file header is missing or incompatible.
 */
void validate_file_header(std::ifstream& input) {
    const std::string columns =
        read_required_line(input, "column header");

    if (normalize_whitespace(columns) != expected_column_header) {
        throw std::runtime_error(
            "Unsupported Afterburner OSCAR column header.");
    }

    const std::string units = read_required_line(input, "unit header");

    if (normalize_whitespace(units) != expected_unit_header) {
        throw std::runtime_error(
            "Unsupported Afterburner OSCAR unit header.");
    }

    const std::string version =
        read_required_line(input, "SMASH version header");

    if (!is_smash_version_line(version)) {
        throw std::runtime_error(
            "Malformed Afterburner SMASH version header.");
    }
}

/**
 * @brief Parses a complete non-negative integer token.
 * @param token Input token to parse.
 * @param field_name Field name used in failure diagnostics.
 * @return Parsed integer value.
 * @throws std::runtime_error If the token is invalid or negative.
 */
int parse_nonnegative_int(
    std::string_view token,
    std::string_view field_name) {
    int value = 0;
    const char* const begin = token.data();
    const char* const end = begin + token.size();
    const auto result = std::from_chars(begin, end, value);

    if (result.ec != std::errc{} || result.ptr != end || value < 0) {
        throw std::runtime_error(
            "Invalid Afterburner " + std::string(field_name) + ".");
    }

    return value;
}

/**
 * @brief Parses a complete particle-count token as `std::size_t`.
 * @param token Input token to parse.
 * @return Parsed non-negative particle count.
 * @throws std::runtime_error If the token is invalid or out of range.
 */
std::size_t parse_particle_count(std::string_view token) {
    std::size_t value = 0;
    const char* const begin = token.data();
    const char* const end = begin + token.size();
    const auto result = std::from_chars(begin, end, value);

    if (result.ec != std::errc{} || result.ptr != end) {
        throw std::runtime_error(
            "Invalid Afterburner subevent particle count.");
    }

    return value;
}

/**
 * @brief Parses one `# event ... ensemble ... out ...` marker.
 * @param line Complete candidate subevent-opening line.
 * @return Parsed raw subevent metadata.
 * @throws std::runtime_error If the opening marker is malformed.
 */
AfterburnerSubeventHeader parse_subevent_header(const std::string& line) {
    std::istringstream stream(line);
    std::string marker;
    std::string event_keyword;
    std::string subevent_token;
    std::string ensemble_keyword;
    std::string ensemble_token;
    std::string out_keyword;
    std::string count_token;
    std::string extra;

    if (!(stream >> marker >> event_keyword >> subevent_token >>
          ensemble_keyword >> ensemble_token >> out_keyword >> count_token) ||
        (stream >> extra)) {
        throw std::runtime_error(
            "Malformed Afterburner subevent opening marker.");
    }

    if (marker != "#" || event_keyword != "event" ||
        ensemble_keyword != "ensemble" || out_keyword != "out") {
        throw std::runtime_error(
            "Malformed Afterburner subevent opening marker.");
    }

    return {
        parse_nonnegative_int(subevent_token, "subevent identifier"),
        parse_nonnegative_int(ensemble_token, "ensemble identifier"),
        parse_particle_count(count_token)
    };
}

/**
 * @brief Parses one complete 22-column Afterburner particle row.
 * @param line Complete candidate particle row.
 * @return Parsed raw particle record.
 * @throws std::runtime_error If the row is malformed or has extra fields.
 */
AfterburnerParticleRecord parse_particle_record(const std::string& line) {
    std::istringstream stream(line);
    AfterburnerParticleRecord record{};
    std::string extra;

    if (!(stream >> record.position.x0 >> record.position.x1 >>
          record.position.x2 >> record.position.x3 >> record.mass >>
          record.momentum.x0 >> record.momentum.x1 >> record.momentum.x2 >>
          record.momentum.x3 >> record.pdg >> record.id >> record.charge >>
          record.ncoll >> record.form_time >> record.xsecfac >>
          record.proc_id_origin >> record.proc_type_origin >>
          record.time_last_coll >> record.pdg_mother1 >> record.pdg_mother2 >>
          record.baryon_number >> record.strangeness) ||
        (stream >> extra)) {
        throw std::runtime_error(
            "Malformed Afterburner particle row.");
    }

    return record;
}

/**
 * @brief Validates one closing marker against the current subevent metadata.
 * @param line Complete candidate subevent-closing line.
 * @param expected Metadata of the currently open subevent.
 * @throws std::runtime_error If the marker is malformed or mismatched.
 */
void validate_subevent_end(
    const std::string& line,
    const AfterburnerSubeventHeader& expected) {
    std::istringstream stream(line);
    std::string marker;
    std::string event_keyword;
    std::string subevent_token;
    std::string ensemble_keyword;
    std::string ensemble_token;
    std::string end_keyword;
    int end_status = 0;
    std::string impact_keyword;
    double impact = 0.0;
    std::string scattering_keyword;
    std::string scattering_value;
    std::string extra;

    if (!(stream >> marker >> event_keyword >> subevent_token >>
          ensemble_keyword >> ensemble_token >> end_keyword >> end_status >>
          impact_keyword >> impact >> scattering_keyword >> scattering_value) ||
        (stream >> extra)) {
        throw std::runtime_error(
            "Malformed Afterburner subevent closing marker.");
    }

    if (marker != "#" || event_keyword != "event" ||
        ensemble_keyword != "ensemble" || end_keyword != "end" ||
        impact_keyword != "impact" ||
        scattering_keyword != "scattering_projectile_target") {
        throw std::runtime_error(
            "Malformed Afterburner subevent closing marker.");
    }

    const int subevent_id =
        parse_nonnegative_int(subevent_token, "closing subevent identifier");
    const int ensemble_id =
        parse_nonnegative_int(ensemble_token, "closing ensemble identifier");

    if (subevent_id != expected.subevent_id ||
        ensemble_id != expected.ensemble_id) {
        throw std::runtime_error(
            "Afterburner subevent closing marker does not match its opening "
            "marker.");
    }
}

}  // namespace

AfterburnerReader::AfterburnerReader(const std::filesystem::path& path)
    : input_(path) {
    if (!input_.is_open()) {
        throw std::runtime_error(
            "Unable to open Afterburner input file: " + path.string());
    }

    validate_file_header(input_);
}

std::optional<AfterburnerSubeventHeader>
AfterburnerReader::begin_next_subevent() {
    if (subevent_open_) {
        throw std::logic_error(
            "Cannot begin a new Afterburner subevent before closing the "
            "current subevent.");
    }

    std::string line;

    if (!std::getline(input_, line)) {
        if (input_.eof()) {
            return std::nullopt;
        }

        throw std::runtime_error(
            "Unable to read the next Afterburner subevent opening marker.");
    }

    current_header_ = parse_subevent_header(line);
    particles_remaining_ = current_header_.particle_count;
    subevent_open_ = true;

    return current_header_;
}

AfterburnerParticleRecord AfterburnerReader::read_particle() {
    if (!subevent_open_) {
        throw std::logic_error(
            "Cannot read an Afterburner particle without an open subevent.");
    }

    if (particles_remaining_ == 0) {
        throw std::logic_error(
            "Cannot read more Afterburner particles than declared by the "
            "current subevent.");
    }

    const std::string line = read_required_line(input_, "particle row");
    const AfterburnerParticleRecord record = parse_particle_record(line);
    --particles_remaining_;

    return record;
}

void AfterburnerReader::finish_subevent() {
    if (!subevent_open_) {
        throw std::logic_error(
            "Cannot finish an Afterburner subevent when none is open.");
    }

    if (particles_remaining_ != 0) {
        throw std::logic_error(
            "Cannot finish an Afterburner subevent before all declared "
            "particle rows have been read.");
    }

    const std::string line =
        read_required_line(input_, "subevent closing marker");
    validate_subevent_end(line, current_header_);

    current_header_ = {};
    particles_remaining_ = 0;
    subevent_open_ = false;
}

}  // namespace input
