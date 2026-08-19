/**
 * @file sampler_reader.cpp
 * @brief Implementation of indexed SMASH Sampler OSCAR input loading.
 */

#include "input/sampler_reader.h"

#include <charconv>
#include <fstream>
#include <functional>
#include <sstream>
#include <stdexcept>
#include <string>
#include <string_view>
#include <system_error>

namespace input {
namespace {

constexpr std::string_view expected_column_header =
    "#!OSCAR2013 particle_lists t x y z mass p0 px py pz pdg ID charge";

constexpr std::string_view expected_unit_header =
    "# Units: fm fm fm fm GeV GeV GeV GeV GeV none none e";

/**
 * @brief Metadata parsed from one Sampler subevent opening marker.
 */
struct SamplerSubeventHeader {
    int subevent_id;             ///< Integer following the `event` token.
    int ensemble_id;             ///< Integer following the `ensemble` token.
    std::size_t particle_count;  ///< Particle rows declared after `out`.
};

/**
 * @brief Sampler row fields retained for construction of the position index.
 *
 * The remaining Sampler columns are parsed and validated but intentionally not
 * returned because the HBT input contract uses Sampler only as a position
 * source.
 */
struct ParsedSamplerParticle {
    common::FourVector position;  ///< Raw `(t, x, y, z)` position in fm.
    int pdg;                      ///< Raw signed PDG particle code.
    int id;                       ///< Raw SMASH particle identifier.
};

/**
 * @brief Reads one mandatory line from the Sampler stream.
 * @param input Input stream positioned before the required line.
 * @param description Human-readable description used in failure diagnostics.
 * @return The complete line without its trailing newline character.
 * @throws std::runtime_error If the required line cannot be read.
 */
std::string read_required_line(
    std::ifstream& input,
    std::string_view description
) {
    std::string line;

    if (!std::getline(input, line)) {
        throw std::runtime_error(
            "Missing or unreadable Sampler " + std::string(description) +
            "."
        );
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
 * @brief Validates and consumes the three-line Sampler OSCAR file header.
 * @param input Newly opened Sampler input stream.
 * @throws std::runtime_error If the file header is missing or incompatible.
 */
void validate_file_header(std::ifstream& input) {
    const std::string columns = read_required_line(input, "column header");

    if (normalize_whitespace(columns) != expected_column_header) {
        throw std::runtime_error(
            "Unsupported Sampler OSCAR column header."
        );
    }

    const std::string units = read_required_line(input, "unit header");

    if (normalize_whitespace(units) != expected_unit_header) {
        throw std::runtime_error(
            "Unsupported Sampler OSCAR unit header."
        );
    }

    const std::string version =
        read_required_line(input, "SMASH version header");

    if (!is_smash_version_line(version)) {
        throw std::runtime_error(
            "Malformed Sampler SMASH version header."
        );
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
    std::string_view field_name
) {
    int value = 0;
    const char* const begin = token.data();
    const char* const end = begin + token.size();
    const auto result = std::from_chars(begin, end, value);

    if (result.ec != std::errc{} || result.ptr != end || value < 0) {
        throw std::runtime_error(
            "Invalid Sampler " + std::string(field_name) + "."
        );
    }

    return value;
}

/**
 * @brief Parses a complete Sampler particle-count token as `std::size_t`.
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
            "Invalid Sampler subevent particle count."
        );
    }

    return value;
}

/**
 * @brief Parses one `# event ... ensemble ... out ...` Sampler marker.
 * @param line Complete candidate subevent-opening line.
 * @return Parsed subevent metadata.
 * @throws std::runtime_error If the opening marker is malformed.
 */
SamplerSubeventHeader parse_subevent_header(const std::string& line) {
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
            "Malformed Sampler subevent opening marker."
        );
    }

    if (marker != "#" || event_keyword != "event" ||
        ensemble_keyword != "ensemble" || out_keyword != "out") {
        throw std::runtime_error(
            "Malformed Sampler subevent opening marker."
        );
    }

    return {
        parse_nonnegative_int(subevent_token, "subevent identifier"),
        parse_nonnegative_int(ensemble_token, "ensemble identifier"),
        parse_particle_count(count_token)
    };
}

/**
 * @brief Parses one complete 12-column Sampler particle row.
 * @param line Complete candidate particle row.
 * @return Position and legacy key fields required by the Sampler index.
 * @throws std::runtime_error If the row is malformed or has extra fields.
 */
ParsedSamplerParticle parse_particle_record(const std::string& line) {
    std::istringstream stream(line);
    ParsedSamplerParticle particle{};
    double mass = 0.0;
    common::FourVector momentum{};
    int charge = 0;
    std::string extra;

    if (!(stream >> particle.position.x0 >> particle.position.x1 >>
          particle.position.x2 >> particle.position.x3 >> mass >>
          momentum.x0 >> momentum.x1 >> momentum.x2 >> momentum.x3 >>
          particle.pdg >> particle.id >> charge) ||
        (stream >> extra)) {
        throw std::runtime_error(
            "Malformed Sampler particle row."
        );
    }

    return particle;
}

/**
 * @brief Validates one Sampler closing marker against its opening metadata.
 * @param line Complete candidate subevent-closing line.
 * @param expected Metadata of the subevent whose rows were just consumed.
 * @throws std::runtime_error If the marker is malformed or mismatched.
 */
void validate_subevent_end(
    const std::string& line,
    const SamplerSubeventHeader& expected
) {
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
            "Malformed Sampler subevent closing marker."
        );
    }

    if (marker != "#" || event_keyword != "event" ||
        ensemble_keyword != "ensemble" || end_keyword != "end" ||
        impact_keyword != "impact" ||
        scattering_keyword != "scattering_projectile_target") {
        throw std::runtime_error(
            "Malformed Sampler subevent closing marker."
        );
    }

    const int subevent_id =
        parse_nonnegative_int(subevent_token, "closing subevent identifier");
    const int ensemble_id =
        parse_nonnegative_int(ensemble_token, "closing ensemble identifier");

    if (subevent_id != expected.subevent_id ||
        ensemble_id != expected.ensemble_id) {
        throw std::runtime_error(
            "Sampler subevent closing marker does not match its opening "
            "marker."
        );
    }
}

}  // namespace

bool SamplerReader::Key::operator==(const Key& other) const noexcept {
    return subevent_id == other.subevent_id && id == other.id &&
           pdg == other.pdg;
}

std::size_t SamplerReader::KeyHash::operator()(
    const Key& key
) const noexcept {
    const std::size_t subevent_hash = std::hash<int>{}(key.subevent_id);
    const std::size_t id_hash = std::hash<int>{}(key.id);
    const std::size_t pdg_hash = std::hash<int>{}(key.pdg);

    return subevent_hash ^ (id_hash << 1U) ^ (pdg_hash << 2U);
}

SamplerReader::SamplerReader(const std::filesystem::path& path) {
    std::ifstream input(path);

    if (!input.is_open()) {
        throw std::runtime_error(
            "Unable to open Sampler input file: " + path.string()
        );
    }

    validate_file_header(input);

    std::string opening_line;

    while (std::getline(input, opening_line)) {
        const SamplerSubeventHeader header =
            parse_subevent_header(opening_line);
        ++subevent_count_;

        for (std::size_t index = 0; index < header.particle_count; ++index) {
            const std::string particle_line =
                read_required_line(input, "particle row");
            const ParsedSamplerParticle particle =
                parse_particle_record(particle_line);
            const Key key{header.subevent_id, particle.id, particle.pdg};

            const auto inserted = positions_.emplace(key, particle.position);

            if (!inserted.second) {
                throw std::runtime_error(
                    "Duplicate Sampler (subevent_id, ID, PDG) key."
                );
            }

            ++particle_count_;
        }

        const std::string closing_line =
            read_required_line(input, "subevent closing marker");
        validate_subevent_end(closing_line, header);
    }

    if (!input.eof()) {
        throw std::runtime_error(
            "Unable to read the next Sampler subevent opening marker."
        );
    }
}

std::optional<common::FourVector> SamplerReader::find_position(
    int subevent_id,
    int id,
    int pdg
) const {
    const auto iterator = positions_.find(Key{subevent_id, id, pdg});

    if (iterator == positions_.end()) {
        return std::nullopt;
    }

    return iterator->second;
}

std::size_t SamplerReader::subevent_count() const noexcept {
    return subevent_count_;
}

std::size_t SamplerReader::particle_count() const noexcept {
    return particle_count_;
}

}  // namespace input
