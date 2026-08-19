/**
 * @file sampler_reader_test.cpp
 * @brief Unit tests for indexed SMASH Sampler OSCAR input loading.
 */

#include "input/sampler_reader.h"

#include <cstdlib>
#include <exception>
#include <filesystem>
#include <fstream>
#include <iostream>
#include <optional>
#include <stdexcept>
#include <string>
#include <system_error>

namespace {

/**
 * @brief Returns the temporary directory used by SamplerReader tests.
 * @return Relative path to the test-data directory.
 */
std::filesystem::path test_directory() {
    return "sampler_reader_test_data";
}

/**
 * @brief Prepares an empty temporary directory for SamplerReader tests.
 * @return `true` when the directory is ready, otherwise `false`.
 */
bool prepare_test_directory() {
    const std::filesystem::path directory = test_directory();
    std::error_code error;

    static_cast<void>(std::filesystem::remove_all(directory, error));
    error.clear();
    static_cast<void>(std::filesystem::create_directories(directory, error));

    if (error) {
        std::cerr
            << "sampler_reader_test: failed to prepare test directory.\n";
        return false;
    }

    return true;
}

/**
 * @brief Removes the temporary test directory without throwing.
 */
void cleanup_test_directory() noexcept {
    std::error_code error;
    static_cast<void>(std::filesystem::remove_all(test_directory(), error));
}

/**
 * @brief Writes complete text to one temporary test file.
 * @param name Filename relative to the test directory.
 * @param content Complete file content to write.
 * @return Path of the created file.
 * @throws std::runtime_error If the file cannot be written.
 */
std::filesystem::path write_test_file(
    const std::string& name,
    const std::string& content
) {
    const std::filesystem::path path = test_directory() / name;
    std::ofstream output(path);

    if (!output) {
        throw std::runtime_error(
            "sampler_reader_test: failed to create test file."
        );
    }

    output << content;

    if (!output) {
        throw std::runtime_error(
            "sampler_reader_test: failed to write test file."
        );
    }

    return path;
}

/**
 * @brief Returns the canonical supported three-line Sampler file header.
 * @return Complete header text including trailing newlines.
 */
std::string supported_header() {
    return
        "#!OSCAR2013 particle_lists t x y z mass p0 px py pz pdg ID charge\n"
        "# Units: fm fm fm fm GeV GeV GeV GeV GeV none none e\n"
        "# SMASH-3.2.1\n";
}

/**
 * @brief Returns one complete valid Sampler particle row.
 * @param t Time component of the stored position.
 * @param x First spatial component of the stored position.
 * @param y Second spatial component of the stored position.
 * @param z Third spatial component of the stored position.
 * @param pdg Signed PDG particle code.
 * @param id Raw particle identifier.
 * @param charge Raw electric charge.
 * @return Complete 12-column row with a trailing newline.
 */
std::string particle_row(
    double t,
    double x,
    double y,
    double z,
    int pdg,
    int id,
    int charge
) {
    return std::to_string(t) + " " + std::to_string(x) + " " +
           std::to_string(y) + " " + std::to_string(z) +
           " 0.138 1.25 0.11 -0.22 1.20 " + std::to_string(pdg) +
           " " + std::to_string(id) + " " + std::to_string(charge) +
           "\n";
}

/**
 * @brief Returns one valid Sampler subevent opening marker.
 * @param subevent_id Subevent identifier to encode.
 * @param ensemble_id Ensemble identifier to encode.
 * @param particle_count Number of rows declared for the subevent.
 * @return Complete opening-marker line with trailing newline.
 */
std::string opening_marker(
    int subevent_id,
    int ensemble_id,
    std::size_t particle_count
) {
    return "# event " + std::to_string(subevent_id) + " ensemble " +
           std::to_string(ensemble_id) + " out " +
           std::to_string(particle_count) + "\n";
}

/**
 * @brief Returns one valid Sampler subevent closing marker.
 * @param subevent_id Subevent identifier to encode.
 * @param ensemble_id Ensemble identifier to encode.
 * @return Complete closing-marker line with trailing newline.
 */
std::string closing_marker(int subevent_id, int ensemble_id) {
    return "# event " + std::to_string(subevent_id) + " ensemble " +
           std::to_string(ensemble_id) +
           " end 0 impact 0.000 scattering_projectile_target yes\n";
}

/**
 * @brief Compares every component of one position with expected values.
 * @param value Position returned by SamplerReader.
 * @param t Expected time component.
 * @param x Expected first spatial component.
 * @param y Expected second spatial component.
 * @param z Expected third spatial component.
 * @return `true` only when all four components match exactly.
 */
bool position_equals(
    const common::FourVector& value,
    double t,
    double x,
    double y,
    double z
) {
    return value.x0 == t && value.x1 == x && value.x2 == y && value.x3 == z;
}

/**
 * @brief Verifies complete loading, counts, lookup, and key separation.
 * @return `true` when a valid two-subevent file is indexed correctly.
 */
bool verify_valid_index_and_lookup() {
    std::string content = supported_header();
    content += opening_marker(0, 0, 3);
    content += particle_row(1.0, 2.0, 3.0, 4.0, 211, 7, 1);
    content += particle_row(5.0, 6.0, 7.0, 8.0, -211, 7, -1);
    content += particle_row(9.0, 10.0, 11.0, 12.0, 211, 8, 1);
    content += closing_marker(0, 0);
    content += opening_marker(1, 4, 1);
    content += particle_row(13.0, 14.0, 15.0, 16.0, 211, 7, 1);
    content += closing_marker(1, 4);

    try {
        const input::SamplerReader reader(
            write_test_file("valid_lookup.oscar", content)
        );

        if (reader.subevent_count() != 2U || reader.particle_count() != 4U) {
            std::cerr
                << "sampler_reader_test: valid file produced wrong counts.\n";
            return false;
        }

        const std::optional<common::FourVector> first =
            reader.find_position(0, 7, 211);
        const std::optional<common::FourVector> different_pdg =
            reader.find_position(0, 7, -211);
        const std::optional<common::FourVector> different_id =
            reader.find_position(0, 8, 211);
        const std::optional<common::FourVector> different_subevent =
            reader.find_position(1, 7, 211);

        if (!first.has_value() ||
            !position_equals(first.value(), 1.0, 2.0, 3.0, 4.0)) {
            std::cerr
                << "sampler_reader_test: primary lookup returned wrong "
                << "position.\n";
            return false;
        }

        if (!different_pdg.has_value() ||
            !position_equals(different_pdg.value(), 5.0, 6.0, 7.0, 8.0)) {
            std::cerr
                << "sampler_reader_test: PDG was not part of the key.\n";
            return false;
        }

        if (!different_id.has_value() ||
            !position_equals(
                different_id.value(),
                9.0,
                10.0,
                11.0,
                12.0
            )) {
            std::cerr
                << "sampler_reader_test: ID was not part of the key.\n";
            return false;
        }

        if (!different_subevent.has_value() ||
            !position_equals(
                different_subevent.value(),
                13.0,
                14.0,
                15.0,
                16.0
            )) {
            std::cerr
                << "sampler_reader_test: subevent was not part of the key.\n";
            return false;
        }

        if (reader.find_position(0, 99, 211).has_value()) {
            std::cerr
                << "sampler_reader_test: absent key unexpectedly matched.\n";
            return false;
        }

        return true;
    } catch (const std::exception& exception) {
        std::cerr
            << "sampler_reader_test: valid lookup case unexpectedly threw: "
            << exception.what()
            << ".\n";
        return false;
    }
}

/**
 * @brief Verifies a supported header can be followed immediately by EOF.
 * @return `true` when an empty Sampler index is accepted.
 */
bool verify_clean_eof_after_header() {
    try {
        const input::SamplerReader reader(
            write_test_file("empty.oscar", supported_header())
        );

        if (reader.subevent_count() != 0U || reader.particle_count() != 0U) {
            std::cerr
                << "sampler_reader_test: empty file produced nonzero counts.\n";
            return false;
        }

        return !reader.find_position(0, 0, 211).has_value();
    } catch (const std::exception& exception) {
        std::cerr
            << "sampler_reader_test: clean EOF unexpectedly threw: "
            << exception.what()
            << ".\n";
        return false;
    }
}

/**
 * @brief Verifies a zero-particle subevent is closed and counted normally.
 * @return `true` when the empty subevent is accepted.
 */
bool verify_zero_particle_subevent() {
    const std::string content =
        supported_header() + opening_marker(3, 2, 0) + closing_marker(3, 2);

    try {
        const input::SamplerReader reader(
            write_test_file("zero_particle.oscar", content)
        );

        return reader.subevent_count() == 1U && reader.particle_count() == 0U;
    } catch (const std::exception& exception) {
        std::cerr
            << "sampler_reader_test: zero-particle subevent unexpectedly "
            << "threw: "
            << exception.what()
            << ".\n";
        return false;
    }
}

/**
 * @brief Verifies token-equivalent whitespace and another SMASH version.
 * @return `true` when supported schema tokens remain accepted.
 */
bool verify_header_whitespace_and_version_flexibility() {
    const std::string content =
        "  #!OSCAR2013   particle_lists\t t x y z mass p0 px py pz pdg ID "
        "charge  \n"
        "# Units:  fm fm fm fm GeV GeV GeV GeV GeV none none e\n"
        "# SMASH-9.9.9\n";

    try {
        const input::SamplerReader reader(
            write_test_file("header_whitespace.oscar", content)
        );

        return reader.subevent_count() == 0U && reader.particle_count() == 0U;
    } catch (const std::exception& exception) {
        std::cerr
            << "sampler_reader_test: token-equivalent header unexpectedly "
            << "threw: "
            << exception.what()
            << ".\n";
        return false;
    }
}

/**
 * @brief Expects SamplerReader construction to reject one file.
 * @param name Unique filename and diagnostic label for the invalid case.
 * @param content Complete invalid Sampler file content.
 * @return `true` only when construction throws `std::runtime_error`.
 */
bool expect_runtime_error(
    const std::string& name,
    const std::string& content
) {
    try {
        const input::SamplerReader reader(write_test_file(name, content));
        static_cast<void>(reader);
    } catch (const std::runtime_error&) {
        return true;
    } catch (const std::exception& exception) {
        std::cerr
            << "sampler_reader_test: "
            << name
            << " threw the wrong exception type: "
            << exception.what()
            << ".\n";
        return false;
    }

    std::cerr
        << "sampler_reader_test: "
        << name
        << " was unexpectedly accepted.\n";
    return false;
}

/**
 * @brief Verifies failure when the requested Sampler file does not exist.
 * @return `true` only when the constructor reports a runtime error.
 */
bool verify_missing_file_is_rejected() {
    try {
        const input::SamplerReader reader(
            test_directory() / "does_not_exist.oscar"
        );
        static_cast<void>(reader);
    } catch (const std::runtime_error&) {
        return true;
    } catch (const std::exception& exception) {
        std::cerr
            << "sampler_reader_test: missing file threw wrong exception: "
            << exception.what()
            << ".\n";
        return false;
    }

    std::cerr
        << "sampler_reader_test: missing file was unexpectedly accepted.\n";
    return false;
}

/**
 * @brief Verifies malformed or incompatible file headers are rejected.
 * @return `true` when every invalid header case fails loading.
 */
bool verify_invalid_headers_are_rejected() {
    bool success = true;

    success = expect_runtime_error(
        "bad_columns.oscar",
        "#!OSCAR2013 particle_lists t x y z mass p0 px py pz pdg ID\n"
        "# Units: fm fm fm fm GeV GeV GeV GeV GeV none none e\n"
        "# SMASH-3.2.1\n"
    ) && success;

    success = expect_runtime_error(
        "bad_units.oscar",
        "#!OSCAR2013 particle_lists t x y z mass p0 px py pz pdg ID charge\n"
        "# Units: fm fm fm fm MeV GeV GeV GeV GeV none none e\n"
        "# SMASH-3.2.1\n"
    ) && success;

    success = expect_runtime_error(
        "bad_version.oscar",
        "#!OSCAR2013 particle_lists t x y z mass p0 px py pz pdg ID charge\n"
        "# Units: fm fm fm fm GeV GeV GeV GeV GeV none none e\n"
        "SMASH-3.2.1\n"
    ) && success;

    success = expect_runtime_error(
        "missing_version.oscar",
        "#!OSCAR2013 particle_lists t x y z mass p0 px py pz pdg ID charge\n"
        "# Units: fm fm fm fm GeV GeV GeV GeV GeV none none e\n"
    ) && success;

    return success;
}

/**
 * @brief Verifies malformed subevent opening markers are rejected.
 * @return `true` when every invalid opening case fails loading.
 */
bool verify_invalid_opening_markers_are_rejected() {
    bool success = true;

    success = expect_runtime_error(
        "data_before_opening.oscar",
        supported_header() + particle_row(1.0, 2.0, 3.0, 4.0, 211, 7, 1)
    ) && success;

    success = expect_runtime_error(
        "bad_open_keyword.oscar",
        supported_header() + "# event 0 ensemble 0 output 0\n"
    ) && success;

    success = expect_runtime_error(
        "negative_subevent.oscar",
        supported_header() + "# event -1 ensemble 0 out 0\n"
    ) && success;

    success = expect_runtime_error(
        "negative_ensemble.oscar",
        supported_header() + "# event 0 ensemble -1 out 0\n"
    ) && success;

    success = expect_runtime_error(
        "negative_count.oscar",
        supported_header() + "# event 0 ensemble 0 out -1\n"
    ) && success;

    success = expect_runtime_error(
        "extra_open_token.oscar",
        supported_header() + "# event 0 ensemble 0 out 0 extra\n"
    ) && success;

    return success;
}

/**
 * @brief Verifies malformed particle rows and declared counts are rejected.
 * @return `true` when every invalid particle-layout case fails loading.
 */
bool verify_invalid_particle_rows_are_rejected() {
    bool success = true;

    success = expect_runtime_error(
        "short_particle.oscar",
        supported_header() + opening_marker(0, 0, 1) +
            "1 2 3 4 0.138 1.25 0.11 -0.22 1.20 211 7\n" +
            closing_marker(0, 0)
    ) && success;

    success = expect_runtime_error(
        "extra_particle.oscar",
        supported_header() + opening_marker(0, 0, 1) +
            particle_row(1.0, 2.0, 3.0, 4.0, 211, 7, 1).substr(
                0,
                particle_row(1.0, 2.0, 3.0, 4.0, 211, 7, 1).size() - 1
            ) +
            " extra\n" + closing_marker(0, 0)
    ) && success;

    success = expect_runtime_error(
        "unexpected_particle_eof.oscar",
        supported_header() + opening_marker(0, 0, 1)
    ) && success;

    success = expect_runtime_error(
        "declared_too_many.oscar",
        supported_header() + opening_marker(0, 0, 2) +
            particle_row(1.0, 2.0, 3.0, 4.0, 211, 7, 1) +
            closing_marker(0, 0)
    ) && success;

    success = expect_runtime_error(
        "declared_too_few.oscar",
        supported_header() + opening_marker(0, 0, 1) +
            particle_row(1.0, 2.0, 3.0, 4.0, 211, 7, 1) +
            particle_row(5.0, 6.0, 7.0, 8.0, 321, 8, 1) +
            closing_marker(0, 0)
    ) && success;

    return success;
}

/**
 * @brief Verifies duplicate legacy lookup keys are fatal input errors.
 * @return `true` when an ambiguous duplicate key is rejected.
 */
bool verify_duplicate_key_is_rejected() {
    const std::string content =
        supported_header() + opening_marker(0, 0, 2) +
        particle_row(1.0, 2.0, 3.0, 4.0, 211, 7, 1) +
        particle_row(5.0, 6.0, 7.0, 8.0, 211, 7, 1) +
        closing_marker(0, 0);

    return expect_runtime_error("duplicate_key.oscar", content);
}

/**
 * @brief Verifies malformed, missing, and mismatched closures are rejected.
 * @return `true` when every invalid closing-marker case fails loading.
 */
bool verify_invalid_closing_markers_are_rejected() {
    bool success = true;

    success = expect_runtime_error(
        "mismatched_close_subevent.oscar",
        supported_header() + opening_marker(0, 0, 0) + closing_marker(1, 0)
    ) && success;

    success = expect_runtime_error(
        "mismatched_close_ensemble.oscar",
        supported_header() + opening_marker(0, 2, 0) + closing_marker(0, 3)
    ) && success;

    success = expect_runtime_error(
        "bad_close_keyword.oscar",
        supported_header() + opening_marker(0, 0, 0) +
            "# event 0 ensemble 0 stop 0 impact 0.0 "
            "scattering_projectile_target yes\n"
    ) && success;

    success = expect_runtime_error(
        "truncated_close.oscar",
        supported_header() + opening_marker(0, 0, 0) +
            "# event 0 ensemble 0 end 0 impact 0.0\n"
    ) && success;

    success = expect_runtime_error(
        "extra_close_token.oscar",
        supported_header() + opening_marker(0, 0, 0) +
            "# event 0 ensemble 0 end 0 impact 0.0 "
            "scattering_projectile_target yes extra\n"
    ) && success;

    success = expect_runtime_error(
        "missing_close.oscar",
        supported_header() + opening_marker(0, 0, 0)
    ) && success;

    return success;
}

}  // namespace

/**
 * @brief Runs the complete SamplerReader unit-test collection.
 * @return `EXIT_SUCCESS` when every test passes, otherwise `EXIT_FAILURE`.
 */
int main() {
    if (!prepare_test_directory()) {
        return EXIT_FAILURE;
    }

    bool success = true;

    success = verify_valid_index_and_lookup() && success;
    success = verify_clean_eof_after_header() && success;
    success = verify_zero_particle_subevent() && success;
    success = verify_header_whitespace_and_version_flexibility() && success;
    success = verify_missing_file_is_rejected() && success;
    success = verify_invalid_headers_are_rejected() && success;
    success = verify_invalid_opening_markers_are_rejected() && success;
    success = verify_invalid_particle_rows_are_rejected() && success;
    success = verify_duplicate_key_is_rejected() && success;
    success = verify_invalid_closing_markers_are_rejected() && success;

    cleanup_test_directory();

    return success ? EXIT_SUCCESS : EXIT_FAILURE;
}
