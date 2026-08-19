/**
 * @file emission_point_resolver_test.cpp
 * @brief Unit tests for HBT emission-position resolution.
 */

#include "input/emission_point_resolver.h"

#include <cmath>
#include <cstdlib>
#include <filesystem>
#include <fstream>
#include <iostream>
#include <limits>
#include <stdexcept>
#include <string>
#include <system_error>

namespace {

/**
 * @brief Return the temporary directory used by resolver tests.
 * @return Relative test-data directory path.
 */
std::filesystem::path test_directory() {
    return "emission_point_resolver_test_data";
}

/**
 * @brief Prepare an empty temporary directory for resolver fixtures.
 * @return `true` when the directory is ready, otherwise `false`.
 */
bool prepare_test_directory() {
    std::error_code error;
    static_cast<void>(
        std::filesystem::remove_all(test_directory(), error)
    );
    error.clear();
    static_cast<void>(
        std::filesystem::create_directories(test_directory(), error)
    );

    if (error) {
        std::cerr
            << "emission_point_resolver_test: failed to prepare directory.\n";
        return false;
    }

    return true;
}

/**
 * @brief Remove the temporary resolver-test directory without throwing.
 */
void cleanup_test_directory() noexcept {
    std::error_code error;
    static_cast<void>(
        std::filesystem::remove_all(test_directory(), error)
    );
}

/**
 * @brief Write one complete temporary test file.
 * @param name Filename relative to the test directory.
 * @param content Complete file content.
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
            "emission_point_resolver_test: failed to create fixture."
        );
    }

    output << content;

    if (!output) {
        throw std::runtime_error(
            "emission_point_resolver_test: failed to write fixture."
        );
    }

    return path;
}

/**
 * @brief Build one valid Sampler file with a single indexed particle.
 * @return Complete Sampler OSCAR text.
 */
std::string sampler_fixture() {
    return
        "#!OSCAR2013 particle_lists t x y z mass p0 px py pz pdg ID charge\n"
        "# Units: fm fm fm fm GeV GeV GeV GeV GeV none none e\n"
        "# SMASH-3.2.1\n"
        "# event 3 ensemble 0 out 1\n"
        "1 2 3 4 0.138 9 0.1 0.2 0.3 211 42 1\n"
        "# event 3 ensemble 0 end 0 impact 0.000 "
        "scattering_projectile_target yes\n";
}

/**
 * @brief Construct one raw Afterburner record for resolver tests.
 * @param position Raw Afterburner position in fm.
 * @param momentum Raw Afterburner momentum in GeV.
 * @param ncoll Raw collision count.
 * @param time_last_coll Raw time of the last interaction in fm.
 * @param id Raw particle identifier used by Sampler lookup.
 * @param pdg Raw signed PDG code used by Sampler lookup.
 * @return Complete raw record with irrelevant fields set to fixed values.
 */
input::AfterburnerParticleRecord make_record(
    const common::FourVector& position,
    const common::FourVector& momentum,
    int ncoll,
    double time_last_coll,
    int id = 42,
    int pdg = 211
) {
    return {
        position,
        0.138,
        momentum,
        pdg,
        id,
        1,
        ncoll,
        0.0,
        1.0,
        0,
        0,
        time_last_coll,
        0,
        0,
        0,
        0
    };
}

/**
 * @brief Compare all four components of two positions exactly.
 * @param actual Position produced by the resolver.
 * @param expected Expected position.
 * @return `true` only when all components are equal.
 */
bool positions_equal(
    const common::FourVector& actual,
    const common::FourVector& expected
) {
    return actual.x0 == expected.x0 && actual.x1 == expected.x1 &&
           actual.x2 == expected.x2 && actual.x3 == expected.x3;
}

/**
 * @brief Verify one successful resolution result.
 * @param result Resolver result to inspect.
 * @param source Expected resolution source.
 * @param position Expected resolved position.
 * @param label Diagnostic label printed on failure.
 * @return true when status, source, and all position components are correct.
 */
bool verify_resolution(
    const input::EmissionPointResolutionResult& result,
    input::EmissionPointSource source,
    const common::FourVector& position,
    const char* label
) {
    if (result.status != input::EmissionPointResolutionStatus::Resolved ||
        !result.candidate.has_value()) {
        std::cerr
            << "emission_point_resolver_test: missing resolution for "
            << label
            << ".\n";
        return false;
    }

    if (result.candidate->source != source ||
        !positions_equal(result.candidate->position, position)) {
        std::cerr
            << "emission_point_resolver_test: wrong resolution for "
            << label
            << ".\n";
        return false;
    }

    return true;
}

/**
 * @brief Verify that Sampler has strict precedence for primordial ncoll zero.
 * @param sampler_reader Sampler index containing the matching legacy key.
 * @return `true` when the Sampler position is used exactly.
 */
bool verify_sampler_precedence(
    const input::SamplerReader& sampler_reader
) {
    const input::AfterburnerParticleRecord particle = make_record(
        {10.0, 5.0, 7.0, 9.0},
        {2.0, 1.0, -0.5, 0.25},
        0,
        std::numeric_limits<double>::quiet_NaN()
    );

    return verify_resolution(
        input::resolve_emission_point(3, particle, true, sampler_reader),
        input::EmissionPointSource::Sampler,
        {1.0, 2.0, 3.0, 4.0},
        "Sampler precedence"
    );
}

/**
 * @brief Verify that a required missing Sampler key has no fallback.
 * @param sampler_reader Sampler index that intentionally lacks the lookup key.
 * @return true when the resolver reports a missing mandatory Sampler key.
 */
bool verify_missing_sampler_match_has_no_fallback(
    const input::SamplerReader& sampler_reader
) {
    const input::AfterburnerParticleRecord particle = make_record(
        {10.0, 5.0, 7.0, 9.0},
        {2.0, 1.0, -0.5, 0.25},
        0,
        4.0,
        99,
        211
    );

    const input::EmissionPointResolutionResult result =
        input::resolve_emission_point(3, particle, true, sampler_reader);

    if (result.status !=
            input::EmissionPointResolutionStatus::MissingMandatorySampler ||
        result.candidate.has_value()) {
        std::cerr
            << "emission_point_resolver_test: missing mandatory Sampler "
            << "match did not preserve the missing-key status.\n";
        return false;
    }

    return true;
}

/**
 * @brief Verify primordial particles with collisions skip the Sampler branch.
 * @param sampler_reader Sampler index containing a matching key.
 * @return `true` when propagation is used instead of Sampler.
 */
bool verify_primordial_with_collisions_propagates(
    const input::SamplerReader& sampler_reader
) {
    const input::AfterburnerParticleRecord particle = make_record(
        {10.0, 5.0, 7.0, 9.0},
        {2.0, 1.0, -0.5, 0.25},
        1,
        4.0
    );

    return verify_resolution(
        input::resolve_emission_point(3, particle, true, sampler_reader),
        input::EmissionPointSource::Propagation,
        {4.0, 2.0, 8.5, 8.25},
        "primordial ncoll non-zero"
    );
}

/**
 * @brief Verify a non-primordial particle never uses the Sampler branch.
 * @param sampler_reader Sampler index containing a matching key.
 * @return `true` when propagation wins despite the matching Sampler key.
 */
bool verify_nonprimordial_particle_propagates(
    const input::SamplerReader& sampler_reader
) {
    const input::AfterburnerParticleRecord particle = make_record(
        {10.0, 5.0, 7.0, 9.0},
        {2.0, 1.0, -0.5, 0.25},
        0,
        4.0
    );

    return verify_resolution(
        input::resolve_emission_point(3, particle, false, sampler_reader),
        input::EmissionPointSource::Propagation,
        {4.0, 2.0, 8.5, 8.25},
        "non-primordial matching key"
    );
}

/**
 * @brief Verify both inclusive time boundaries permit propagation.
 * @param sampler_reader Sampler index unused by these non-primordial cases.
 * @return `true` when times zero and t_afterburner both propagate.
 */
bool verify_propagation_time_boundaries(
    const input::SamplerReader& sampler_reader
) {
    const input::AfterburnerParticleRecord zero_time = make_record(
        {10.0, 5.0, 7.0, 9.0},
        {2.0, 1.0, -0.5, 0.25},
        2,
        0.0
    );

    if (!verify_resolution(
            input::resolve_emission_point(
                3,
                zero_time,
                false,
                sampler_reader
            ),
            input::EmissionPointSource::Propagation,
            {0.0, 0.0, 9.5, 7.75},
            "time_last_coll zero boundary")) {
        return false;
    }

    const input::AfterburnerParticleRecord output_time = make_record(
        {10.0, 5.0, 7.0, 9.0},
        {2.0, 1.0, -0.5, 0.25},
        2,
        10.0
    );

    return verify_resolution(
        input::resolve_emission_point(
            3,
            output_time,
            false,
            sampler_reader
        ),
        input::EmissionPointSource::Propagation,
        {10.0, 5.0, 7.0, 9.0},
        "time_last_coll output-time boundary"
    );
}

/**
 * @brief Verify invalid propagation times use raw Afterburner positions.
 * @param sampler_reader Sampler index unused by these non-primordial cases.
 * @return `true` when every invalid propagation-time case falls back exactly.
 */
bool verify_propagation_time_fallbacks(
    const input::SamplerReader& sampler_reader
) {
    const common::FourVector position{10.0, 5.0, 7.0, 9.0};
    const common::FourVector momentum{2.0, 1.0, -0.5, 0.25};

    const input::AfterburnerParticleRecord negative_time = make_record(
        position,
        momentum,
        1,
        -0.5
    );
    if (!verify_resolution(
            input::resolve_emission_point(
                3,
                negative_time,
                false,
                sampler_reader
            ),
            input::EmissionPointSource::Afterburner,
            position,
            "negative time fallback")) {
        return false;
    }

    const input::AfterburnerParticleRecord late_time = make_record(
        position,
        momentum,
        1,
        10.5
    );
    if (!verify_resolution(
            input::resolve_emission_point(
                3,
                late_time,
                false,
                sampler_reader
            ),
            input::EmissionPointSource::Afterburner,
            position,
            "late time fallback")) {
        return false;
    }

    const double nan = std::numeric_limits<double>::quiet_NaN();
    const input::AfterburnerParticleRecord nan_time = make_record(
        position,
        momentum,
        1,
        nan
    );
    if (!verify_resolution(
            input::resolve_emission_point(
                3,
                nan_time,
                false,
                sampler_reader
            ),
            input::EmissionPointSource::Afterburner,
            position,
            "NaN time fallback")) {
        return false;
    }

    const double infinity = std::numeric_limits<double>::infinity();
    const input::AfterburnerParticleRecord positive_infinite_time = make_record(
        position,
        momentum,
        1,
        infinity
    );
    if (!verify_resolution(
            input::resolve_emission_point(
                3,
                positive_infinite_time,
                false,
                sampler_reader
            ),
            input::EmissionPointSource::Afterburner,
            position,
            "positive infinite time fallback")) {
        return false;
    }

    const input::AfterburnerParticleRecord negative_infinite_time = make_record(
        position,
        momentum,
        1,
        -infinity
    );

    return verify_resolution(
        input::resolve_emission_point(
            3,
            negative_infinite_time,
            false,
            sampler_reader
        ),
        input::EmissionPointSource::Afterburner,
        position,
        "negative infinite time fallback"
    );
}

/**
 * @brief Verify non-finite propagated positions are returned without fallback.
 * @param sampler_reader Sampler index unused by the non-primordial case.
 * @return true when the invalid propagation candidate and source are retained.
 */
bool verify_nonfinite_propagated_position_rejected(
    const input::SamplerReader& sampler_reader
) {
    const double maximum = std::numeric_limits<double>::max();
    const input::AfterburnerParticleRecord particle = make_record(
        {maximum, maximum, 0.0, 0.0},
        {1.0, -1.0, 0.0, 0.0},
        1,
        0.0
    );
    const input::EmissionPointResolutionResult result =
        input::resolve_emission_point(3, particle, false, sampler_reader);

    return result.status ==
               input::EmissionPointResolutionStatus::NonFinitePosition &&
           result.candidate.has_value() &&
           result.candidate->source ==
               input::EmissionPointSource::Propagation &&
           !std::isfinite(result.candidate->position.x1);
}

/**
 * @brief Verify non-finite raw Afterburner positions are returned unchanged.
 * @param sampler_reader Sampler index unused by the non-primordial case.
 * @return true when the invalid raw candidate is retained without fallback.
 */
bool verify_nonfinite_afterburner_position_rejected(
    const input::SamplerReader& sampler_reader
) {
    const double nan = std::numeric_limits<double>::quiet_NaN();
    const input::AfterburnerParticleRecord particle = make_record(
        {10.0, 5.0, nan, 9.0},
        {2.0, 1.0, -0.5, 0.25},
        1,
        -0.5
    );
    const input::EmissionPointResolutionResult result =
        input::resolve_emission_point(3, particle, false, sampler_reader);

    return result.status ==
               input::EmissionPointResolutionStatus::NonFinitePosition &&
           result.candidate.has_value() &&
           result.candidate->source ==
               input::EmissionPointSource::Afterburner &&
           std::isnan(result.candidate->position.x2);
}

}  // namespace

/**
 * @brief Run the complete emission-point-resolver unit-test collection.
 * @return `EXIT_SUCCESS` when every test passes, otherwise `EXIT_FAILURE`.
 */
int main() {
    if (!prepare_test_directory()) {
        return EXIT_FAILURE;
    }

    bool success = true;

    try {
        const input::SamplerReader sampler_reader(
            write_test_file("sampler.oscar", sampler_fixture())
        );

        success = verify_sampler_precedence(sampler_reader) && success;
        success = verify_missing_sampler_match_has_no_fallback(
            sampler_reader
        ) && success;
        success = verify_primordial_with_collisions_propagates(
            sampler_reader
        ) && success;
        success = verify_nonprimordial_particle_propagates(
            sampler_reader
        ) && success;
        success = verify_propagation_time_boundaries(sampler_reader) &&
                  success;
        success = verify_propagation_time_fallbacks(sampler_reader) && success;
        success = verify_nonfinite_propagated_position_rejected(
            sampler_reader
        ) && success;
        success = verify_nonfinite_afterburner_position_rejected(
            sampler_reader
        ) && success;
    } catch (const std::exception& error) {
        std::cerr
            << "emission_point_resolver_test: unexpected exception: "
            << error.what()
            << ".\n";
        success = false;
    }

    cleanup_test_directory();
    return success ? EXIT_SUCCESS : EXIT_FAILURE;
}
