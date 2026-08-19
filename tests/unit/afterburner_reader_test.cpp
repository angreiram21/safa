/**
 * @file afterburner_reader_test.cpp
 * @brief Unit tests for streaming SMASH Afterburner OSCAR input.
 */

#include "input/afterburner_reader.h"

#include <chrono>
#include <cstdlib>
#include <exception>
#include <filesystem>
#include <fstream>
#include <iostream>
#include <stdexcept>
#include <string>
#include <string_view>
#include <system_error>

namespace {

constexpr std::string_view valid_column_header =
    "#!OSCAR2013Extended particle_lists "
    "t x y z mass p0 px py pz pdg ID charge ncoll form_time xsecfac "
    "proc_id_origin proc_type_origin time_last_coll pdg_mother1 pdg_mother2 "
    "baryon_number strangeness\n";

constexpr std::string_view valid_unit_header =
    "# Units: fm fm fm fm GeV GeV GeV GeV GeV none none e none fm none "
    "none none fm none none none none\n";

constexpr std::string_view valid_version_header = "# SMASH-3.2.1\n";

constexpr std::string_view valid_particle_row =
    "1.25 -2.5 3.75 -4 0.5 6.25 -7.5 8.75 -9 211 42 1 3 10.5 "
    "1.25 13 5 11.75 223 -211 0 -1\n";

/**
 * @brief Creates a unique temporary directory for reader fixtures.
 * @return Path to the newly created temporary directory.
 */
std::filesystem::path make_fixture_directory() {
    const auto stamp =
        std::chrono::steady_clock::now().time_since_epoch().count();
    const std::filesystem::path directory =
        std::filesystem::temp_directory_path() /
        ("afterburner_reader_test_" + std::to_string(stamp));

    std::filesystem::create_directories(directory);
    return directory;
}

/**
 * @brief Writes one complete text fixture into the temporary directory.
 * @param directory Directory that owns the fixture.
 * @param name Fixture filename.
 * @param content Complete file content.
 * @return Path to the created fixture file.
 */
std::filesystem::path write_fixture(
    const std::filesystem::path& directory,
    std::string_view name,
    std::string_view content) {
    const std::filesystem::path path = directory / std::string(name);
    std::ofstream output(path);

    if (!output.is_open()) {
        throw std::runtime_error("Unable to create test fixture.");
    }

    output << content;
    output.close();

    if (!output) {
        throw std::runtime_error("Unable to write test fixture.");
    }

    return path;
}

/**
 * @brief Builds a complete valid OSCAR prologue.
 * @param version_header Version line to append after schema and units.
 * @return Complete three-line prologue.
 */
std::string make_valid_prologue(
    std::string_view version_header = valid_version_header) {
    return std::string(valid_column_header) + std::string(valid_unit_header) +
           std::string(version_header);
}

/**
 * @brief Verifies that an action throws the requested exception type.
 * @tparam Exception Expected exception type.
 * @tparam Callable Callable action type.
 * @param action Operation expected to throw.
 * @param label Description reported on failure.
 * @return `true` if the expected exception type was thrown.
 */
template <typename Exception, typename Callable>
bool expect_exception(Callable&& action, std::string_view label) {
    try {
        action();
    } catch (const Exception&) {
        return true;
    } catch (const std::exception& error) {
        std::cerr << "Unexpected exception for " << label << ": "
                  << error.what() << '\n';
        return false;
    } catch (...) {
        std::cerr << "Unexpected non-standard exception for " << label
                  << '\n';
        return false;
    }

    std::cerr << "Expected exception was not thrown for " << label << '\n';
    return false;
}

/**
 * @brief Verifies constructor acceptance and rejection of file headers.
 * @param directory Temporary fixture directory.
 * @return `true` if every constructor contract case passes.
 */
bool verify_constructor_contract(const std::filesystem::path& directory) {
    bool success = true;

    const std::filesystem::path valid = write_fixture(
        directory,
        "valid_header.oscar",
        make_valid_prologue());

    try {
        const input::AfterburnerReader reader(valid);
    } catch (const std::exception& error) {
        std::cerr << "Valid header was rejected: " << error.what() << '\n';
        success = false;
    }

    const std::filesystem::path newer_version = write_fixture(
        directory,
        "newer_version.oscar",
        make_valid_prologue("# SMASH-9.4.0\n"));

    try {
        const input::AfterburnerReader reader(newer_version);
    } catch (const std::exception& error) {
        std::cerr << "Compatible SMASH version was rejected: "
                  << error.what() << '\n';
        success = false;
    }

    const std::filesystem::path flexible_whitespace = write_fixture(
        directory,
        "flexible_whitespace.oscar",
        "  #!OSCAR2013Extended   particle_lists t x y z mass p0 px py pz "
        "pdg ID charge ncoll form_time xsecfac proc_id_origin "
        "proc_type_origin time_last_coll pdg_mother1 pdg_mother2 "
        "baryon_number strangeness  \n"
        "#   Units: fm fm fm fm GeV GeV GeV GeV GeV none none e none fm "
        "none none none fm none none none none\n"
        "#   SMASH-3.2.1\n");

    try {
        const input::AfterburnerReader reader(flexible_whitespace);
    } catch (const std::exception& error) {
        std::cerr << "Whitespace-equivalent header was rejected: "
                  << error.what() << '\n';
        success = false;
    }

    const std::filesystem::path missing = directory / "missing.oscar";
    success &= expect_exception<std::runtime_error>(
        [&missing]() {
            const input::AfterburnerReader reader(missing);
        },
        "missing input file");

    const std::filesystem::path bad_columns = write_fixture(
        directory,
        "bad_columns.oscar",
        "#!OSCAR2013Extended particle_lists t x y z\n");
    success &= expect_exception<std::runtime_error>(
        [&bad_columns]() {
            const input::AfterburnerReader reader(bad_columns);
        },
        "incompatible column header");

    const std::filesystem::path bad_units = write_fixture(
        directory,
        "bad_units.oscar",
        std::string(valid_column_header) +
            "# Units: incompatible\n" +
            std::string(valid_version_header));
    success &= expect_exception<std::runtime_error>(
        [&bad_units]() {
            const input::AfterburnerReader reader(bad_units);
        },
        "incompatible unit header");

    const std::filesystem::path bad_version = write_fixture(
        directory,
        "bad_version.oscar",
        make_valid_prologue("# OTHER-3.2.1\n"));
    success &= expect_exception<std::runtime_error>(
        [&bad_version]() {
            const input::AfterburnerReader reader(bad_version);
        },
        "malformed SMASH version header");

    return success;
}

/**
 * @brief Verifies parsing and state rules for subevent opening markers.
 * @param directory Temporary fixture directory.
 * @return `true` if every subevent-opening case passes.
 */
bool verify_begin_next_subevent(const std::filesystem::path& directory) {
    bool success = true;

    const std::filesystem::path valid = write_fixture(
        directory,
        "valid_subevent.oscar",
        make_valid_prologue() + "# event 17 ensemble 3 out 424\n");
    input::AfterburnerReader reader(valid);
    const auto header = reader.begin_next_subevent();

    if (!header.has_value()) {
        std::cerr << "Valid subevent marker returned clean EOF.\n";
        success = false;
    } else if (header->subevent_id != 17 || header->ensemble_id != 3 ||
               header->particle_count != 424) {
        std::cerr << "Valid subevent metadata was parsed incorrectly.\n";
        success = false;
    }

    success &= expect_exception<std::logic_error>(
        [&reader]() {
            static_cast<void>(reader.begin_next_subevent());
        },
        "beginning a second subevent while one is open");

    const std::filesystem::path empty = write_fixture(
        directory,
        "clean_eof.oscar",
        make_valid_prologue());
    input::AfterburnerReader empty_reader(empty);

    if (empty_reader.begin_next_subevent().has_value()) {
        std::cerr << "Header-only file did not produce clean EOF.\n";
        success = false;
    }

    const std::filesystem::path malformed = write_fixture(
        directory,
        "malformed_marker.oscar",
        make_valid_prologue() + "# event 0 ensemble 0 output 2\n");
    success &= expect_exception<std::runtime_error>(
        [&malformed]() {
            input::AfterburnerReader malformed_reader(malformed);
            static_cast<void>(malformed_reader.begin_next_subevent());
        },
        "malformed subevent opening marker");

    const std::filesystem::path negative_id = write_fixture(
        directory,
        "negative_id.oscar",
        make_valid_prologue() + "# event -1 ensemble 0 out 2\n");
    success &= expect_exception<std::runtime_error>(
        [&negative_id]() {
            input::AfterburnerReader negative_reader(negative_id);
            static_cast<void>(negative_reader.begin_next_subevent());
        },
        "negative subevent identifier");

    const std::filesystem::path negative_count = write_fixture(
        directory,
        "negative_count.oscar",
        make_valid_prologue() + "# event 0 ensemble 0 out -1\n");
    success &= expect_exception<std::runtime_error>(
        [&negative_count]() {
            input::AfterburnerReader negative_reader(negative_count);
            static_cast<void>(negative_reader.begin_next_subevent());
        },
        "negative particle count");

    const std::filesystem::path trailing = write_fixture(
        directory,
        "trailing_token.oscar",
        make_valid_prologue() + "# event 0 ensemble 0 out 2 unexpected\n");
    success &= expect_exception<std::runtime_error>(
        [&trailing]() {
            input::AfterburnerReader trailing_reader(trailing);
            static_cast<void>(trailing_reader.begin_next_subevent());
        },
        "trailing opening-marker token");

    return success;
}

/**
 * @brief Verifies parsing and state rules for raw particle rows.
 * @param directory Temporary fixture directory.
 * @return `true` if every particle-reading case passes.
 */
bool verify_read_particle(const std::filesystem::path& directory) {
    bool success = true;

    const std::filesystem::path valid = write_fixture(
        directory,
        "valid_particle.oscar",
        make_valid_prologue() +
            "# event 0 ensemble 0 out 1\n" +
            std::string(valid_particle_row) +
            "# event 0 ensemble 0 end 0 impact -1.000 "
            "scattering_projectile_target no\n");
    input::AfterburnerReader reader(valid);

    success &= expect_exception<std::logic_error>(
        [&reader]() {
            static_cast<void>(reader.read_particle());
        },
        "reading a particle without an open subevent");

    static_cast<void>(reader.begin_next_subevent());
    const input::AfterburnerParticleRecord record = reader.read_particle();

    if (record.position.x0 != 1.25 || record.position.x1 != -2.5 ||
        record.position.x2 != 3.75 || record.position.x3 != -4.0 ||
        record.mass != 0.5 || record.momentum.x0 != 6.25 ||
        record.momentum.x1 != -7.5 || record.momentum.x2 != 8.75 ||
        record.momentum.x3 != -9.0 || record.pdg != 211 || record.id != 42 ||
        record.charge != 1 || record.ncoll != 3 ||
        record.form_time != 10.5 || record.xsecfac != 1.25 ||
        record.proc_id_origin != 13 || record.proc_type_origin != 5 ||
        record.time_last_coll != 11.75 || record.pdg_mother1 != 223 ||
        record.pdg_mother2 != -211 || record.baryon_number != 0 ||
        record.strangeness != -1) {
        std::cerr << "Valid particle row was parsed incorrectly.\n";
        success = false;
    }

    success &= expect_exception<std::logic_error>(
        [&reader]() {
            static_cast<void>(reader.read_particle());
        },
        "reading beyond the declared particle count");

    const std::filesystem::path short_row = write_fixture(
        directory,
        "short_particle.oscar",
        make_valid_prologue() +
            "# event 0 ensemble 0 out 1\n"
            "1 2 3\n");
    success &= expect_exception<std::runtime_error>(
        [&short_row]() {
            input::AfterburnerReader short_reader(short_row);
            static_cast<void>(short_reader.begin_next_subevent());
            static_cast<void>(short_reader.read_particle());
        },
        "particle row with too few fields");

    const std::filesystem::path extra_field = write_fixture(
        directory,
        "extra_particle_field.oscar",
        make_valid_prologue() +
            "# event 0 ensemble 0 out 1\n" +
            std::string(valid_particle_row).substr(
                0,
                valid_particle_row.size() - 1) +
            " unexpected\n");
    success &= expect_exception<std::runtime_error>(
        [&extra_field]() {
            input::AfterburnerReader extra_reader(extra_field);
            static_cast<void>(extra_reader.begin_next_subevent());
            static_cast<void>(extra_reader.read_particle());
        },
        "particle row with an extra field");

    const std::filesystem::path missing_row = write_fixture(
        directory,
        "missing_particle.oscar",
        make_valid_prologue() + "# event 0 ensemble 0 out 1\n");
    success &= expect_exception<std::runtime_error>(
        [&missing_row]() {
            input::AfterburnerReader missing_reader(missing_row);
            static_cast<void>(missing_reader.begin_next_subevent());
            static_cast<void>(missing_reader.read_particle());
        },
        "unexpected EOF while reading a particle row");

    return success;
}

/**
 * @brief Verifies closing-marker validation and subevent state transitions.
 * @param directory Temporary fixture directory.
 * @return `true` if every subevent-closing case passes.
 */
bool verify_finish_subevent(const std::filesystem::path& directory) {
    bool success = true;

    const std::filesystem::path valid = write_fixture(
        directory,
        "valid_finish.oscar",
        make_valid_prologue() +
            "# event 0 ensemble 0 out 0\n"
            "# event 0 ensemble 0 end 0 impact -1.000 "
            "scattering_projectile_target no\n"
            "# event 1 ensemble 0 out 0\n"
            "# event 1 ensemble 0 end 0 impact -1.000 "
            "scattering_projectile_target no\n");
    input::AfterburnerReader reader(valid);

    success &= expect_exception<std::logic_error>(
        [&reader]() {
            reader.finish_subevent();
        },
        "finishing without an open subevent");

    const auto first = reader.begin_next_subevent();
    if (!first.has_value() || first->subevent_id != 0) {
        std::cerr << "First zero-particle subevent was not opened correctly.\n";
        success = false;
    }

    reader.finish_subevent();

    const auto second = reader.begin_next_subevent();
    if (!second.has_value() || second->subevent_id != 1) {
        std::cerr << "Reader did not advance after finishing a subevent.\n";
        success = false;
    }

    reader.finish_subevent();

    if (reader.begin_next_subevent().has_value()) {
        std::cerr << "Reader did not reach clean EOF after final subevent.\n";
        success = false;
    }

    const std::filesystem::path unread = write_fixture(
        directory,
        "unread_particle.oscar",
        make_valid_prologue() +
            "# event 0 ensemble 0 out 1\n" +
            std::string(valid_particle_row));
    input::AfterburnerReader unread_reader(unread);
    static_cast<void>(unread_reader.begin_next_subevent());
    success &= expect_exception<std::logic_error>(
        [&unread_reader]() {
            unread_reader.finish_subevent();
        },
        "finishing with unread particle rows");

    const std::filesystem::path wrong_subevent = write_fixture(
        directory,
        "wrong_closing_subevent.oscar",
        make_valid_prologue() +
            "# event 2 ensemble 4 out 0\n"
            "# event 3 ensemble 4 end 0 impact -1.000 "
            "scattering_projectile_target no\n");
    success &= expect_exception<std::runtime_error>(
        [&wrong_subevent]() {
            input::AfterburnerReader wrong_reader(wrong_subevent);
            static_cast<void>(wrong_reader.begin_next_subevent());
            wrong_reader.finish_subevent();
        },
        "mismatched closing subevent identifier");

    const std::filesystem::path wrong_ensemble = write_fixture(
        directory,
        "wrong_closing_ensemble.oscar",
        make_valid_prologue() +
            "# event 2 ensemble 4 out 0\n"
            "# event 2 ensemble 5 end 0 impact -1.000 "
            "scattering_projectile_target no\n");
    success &= expect_exception<std::runtime_error>(
        [&wrong_ensemble]() {
            input::AfterburnerReader wrong_reader(wrong_ensemble);
            static_cast<void>(wrong_reader.begin_next_subevent());
            wrong_reader.finish_subevent();
        },
        "mismatched closing ensemble identifier");

    const std::filesystem::path malformed = write_fixture(
        directory,
        "malformed_closing_marker.oscar",
        make_valid_prologue() +
            "# event 0 ensemble 0 out 0\n"
            "# event 0 ensemble 0 stop 0 impact -1.000\n");
    success &= expect_exception<std::runtime_error>(
        [&malformed]() {
            input::AfterburnerReader malformed_reader(malformed);
            static_cast<void>(malformed_reader.begin_next_subevent());
            malformed_reader.finish_subevent();
        },
        "malformed subevent closing marker");

    const std::filesystem::path truncated = write_fixture(
        directory,
        "truncated_closing_marker.oscar",
        make_valid_prologue() +
            "# event 0 ensemble 0 out 0\n"
            "# event 0 ensemble 0 end\n");
    success &= expect_exception<std::runtime_error>(
        [&truncated]() {
            input::AfterburnerReader truncated_reader(truncated);
            static_cast<void>(truncated_reader.begin_next_subevent());
            truncated_reader.finish_subevent();
        },
        "truncated subevent closing marker");

    const std::filesystem::path missing = write_fixture(
        directory,
        "missing_closing_marker.oscar",
        make_valid_prologue() + "# event 0 ensemble 0 out 0\n");
    success &= expect_exception<std::runtime_error>(
        [&missing]() {
            input::AfterburnerReader missing_reader(missing);
            static_cast<void>(missing_reader.begin_next_subevent());
            missing_reader.finish_subevent();
        },
        "unexpected EOF while reading a closing marker");

    return success;
}

}  // namespace

/**
 * @brief Runs the Afterburner streaming-reader unit tests.
 * @return `EXIT_SUCCESS` when all checks pass, otherwise `EXIT_FAILURE`.
 */
int main() {
    const std::filesystem::path directory = make_fixture_directory();
    bool success = true;

    try {
        success &= verify_constructor_contract(directory);
        success &= verify_begin_next_subevent(directory);
        success &= verify_read_particle(directory);
        success &= verify_finish_subevent(directory);
    } catch (const std::exception& error) {
        std::cerr << "Unexpected test setup failure: " << error.what() << '\n';
        success = false;
    }

    std::error_code cleanup_error;
    std::filesystem::remove_all(directory, cleanup_error);

    if (cleanup_error) {
        std::cerr << "Unable to remove temporary fixture directory: "
                  << cleanup_error.message() << '\n';
        success = false;
    }

    return success ? EXIT_SUCCESS : EXIT_FAILURE;
}
