cmake_minimum_required(VERSION 3.16)

project(smash-afterburner-analysis-modular LANGUAGES CXX)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)

include(CTest)
enable_testing()

find_package(Doxygen REQUIRED COMPONENTS dot)

add_custom_target(doxygen
    COMMAND
        ${CMAKE_COMMAND}
        "-DSOURCE_DIR=${PROJECT_SOURCE_DIR}"
        -P "${PROJECT_SOURCE_DIR}/cmake/run_doxygen.cmake"
    WORKING_DIRECTORY "${PROJECT_SOURCE_DIR}"
    COMMENT "Generating Doxygen documentation in docs/doxygen"
    VERBATIM
    USES_TERMINAL
)

add_test(
    NAME doxygen_contract
    COMMAND
        ${CMAKE_COMMAND}
        "-DSOURCE_DIR=${PROJECT_SOURCE_DIR}"
        -P "${PROJECT_SOURCE_DIR}/cmake/run_doxygen.cmake"
)

set(
    HBT_REAL_REGRESSION_MAIN_CONFIG
    ""
    CACHE FILEPATH
    "Optional main.yaml for the real-fixture event-preparation regression."
)

# External dependencies must be installed before configuration.
# This project never downloads third-party packages or source code.
find_package(yaml-cpp CONFIG REQUIRED)
find_package(Threads REQUIRED)

if(TARGET yaml-cpp::yaml-cpp)
    set(YAML_CPP_TARGET yaml-cpp::yaml-cpp)
elseif(TARGET yaml-cpp)
    set(YAML_CPP_TARGET yaml-cpp)
else()
    message(
        FATAL_ERROR
        "yaml-cpp was found, but no supported CMake target is available."
    )
endif()

find_package(ROOT REQUIRED COMPONENTS Minuit2)

if(NOT TARGET ROOT::Minuit2)
    message(
        FATAL_ERROR
        "ROOT was found, but ROOT::Minuit2 is unavailable."
    )
endif()

add_library(kinematics STATIC
    src/common/kinematics.cpp
)

target_include_directories(kinematics
    PUBLIC
        ${PROJECT_SOURCE_DIR}/src
)

target_compile_options(kinematics
    PRIVATE
        -Wall
        -Wextra
        -Wpedantic
        -Werror
)

add_library(kinematics_validation STATIC
    src/common/kinematics_validation.cpp
)

target_include_directories(kinematics_validation
    PUBLIC
        ${PROJECT_SOURCE_DIR}/src
)

target_compile_options(kinematics_validation
    PRIVATE
        -Wall
        -Wextra
        -Wpedantic
        -Werror
)

add_library(afterburner_reader STATIC
    src/input/afterburner_reader.cpp
)

target_include_directories(afterburner_reader
    PUBLIC
        ${PROJECT_SOURCE_DIR}/src
)

target_compile_options(afterburner_reader
    PRIVATE
        -Wall
        -Wextra
        -Wpedantic
        -Werror
)

add_library(sampler_reader STATIC
    src/input/sampler_reader.cpp
)

target_include_directories(sampler_reader
    PUBLIC
        ${PROJECT_SOURCE_DIR}/src
)

target_compile_options(sampler_reader
    PRIVATE
        -Wall
        -Wextra
        -Wpedantic
        -Werror
)

add_library(emission_point_resolver STATIC
    src/input/emission_point_resolver.cpp
)

target_include_directories(emission_point_resolver
    PUBLIC
        ${PROJECT_SOURCE_DIR}/src
)

target_link_libraries(emission_point_resolver
    PUBLIC
        afterburner_reader
        sampler_reader
)

target_compile_options(emission_point_resolver
    PRIVATE
        -Wall
        -Wextra
        -Wpedantic
        -Werror
)

add_library(hbt_species STATIC
    src/hbt/species/species.cpp
)

target_include_directories(hbt_species
    PUBLIC
        ${PROJECT_SOURCE_DIR}/src
)

target_compile_options(hbt_species
    PRIVATE
        -Wall
        -Wextra
        -Wpedantic
        -Werror
)

add_library(hbt_channels STATIC
    src/hbt/channels/channel_catalog.cpp
)

target_include_directories(hbt_channels
    PUBLIC
        ${PROJECT_SOURCE_DIR}/src
)

target_link_libraries(hbt_channels
    PUBLIC
        hbt_species
)

target_compile_options(hbt_channels
    PRIVATE
        -Wall
        -Wextra
        -Wpedantic
        -Werror
)

add_library(hbt_selection STATIC
    src/hbt/selection/required_channels.cpp
    src/hbt/selection/required_species.cpp
    src/hbt/selection/species_requirement.cpp
)

target_include_directories(hbt_selection
    PUBLIC
        ${PROJECT_SOURCE_DIR}/src
)

target_link_libraries(hbt_selection
    PUBLIC
        hbt_channels
)

target_compile_options(hbt_selection
    PRIVATE
        -Wall
        -Wextra
        -Wpedantic
        -Werror
)

add_library(hbt_config STATIC
    src/hbt/config/primitive_channel_token.cpp
    src/hbt/config/analysis_product_parser.cpp
    src/hbt/config/hbt_enabled_channels_parser.cpp
    src/hbt/config/hbt_config_loader.cpp
)

target_include_directories(hbt_config
    PUBLIC
        ${PROJECT_SOURCE_DIR}/src
)

target_link_libraries(hbt_config
    PUBLIC
        hbt_channels
        hbt_selection
    PRIVATE
        ${YAML_CPP_TARGET}
)

target_compile_options(hbt_config
    PRIVATE
        -Wall
        -Wextra
        -Wpedantic
        -Werror
)

add_library(hbt_particle_selector STATIC
    src/hbt/event/particle_selector.cpp
)

target_include_directories(hbt_particle_selector
    PUBLIC
        ${PROJECT_SOURCE_DIR}/src
)

target_link_libraries(hbt_particle_selector
    PUBLIC
        hbt_config
        hbt_species
    PRIVATE
        kinematics
        kinematics_validation
)

target_compile_options(hbt_particle_selector
    PRIVATE
        -Wall
        -Wextra
        -Wpedantic
        -Werror
)

add_library(hbt_event_buffers STATIC
    src/hbt/event/event_buffers.cpp
    src/hbt/event/subevent_shuffle.cpp
)

target_include_directories(hbt_event_buffers
    PUBLIC
        ${PROJECT_SOURCE_DIR}/src
)

target_compile_options(hbt_event_buffers
    PRIVATE
        -Wall
        -Wextra
        -Wpedantic
        -Werror
)

add_library(hbt_pair STATIC
    src/hbt/pair/pair_count_accumulator.cpp
    src/hbt/pair/pair_kinematics.cpp
    src/hbt/pair/pair_kinematics_validation.cpp
)

target_include_directories(hbt_pair
    PUBLIC
        ${PROJECT_SOURCE_DIR}/src
)

target_link_libraries(hbt_pair
    PUBLIC
        hbt_channels
        hbt_event_buffers
)

target_compile_options(hbt_pair
    PRIVATE
        -Wall
        -Wextra
        -Wpedantic
        -Werror
)

add_library(hbt_boosts STATIC
    src/hbt/boosts/lab_to_lcms.cpp
    src/hbt/boosts/lcms_to_prf.cpp
)

target_include_directories(hbt_boosts
    PUBLIC
        ${PROJECT_SOURCE_DIR}/src
)

target_link_libraries(hbt_boosts
    PUBLIC
        hbt_pair
)

target_compile_options(hbt_boosts
    PRIVATE
        -Wall
        -Wextra
        -Wpedantic
        -Werror
)

add_library(hbt_pair_frames STATIC
    src/hbt/pair/pair_frame_observables.cpp
    src/hbt/pair/pair_frame_observables_validation.cpp
)

target_include_directories(hbt_pair_frames
    PUBLIC
        ${PROJECT_SOURCE_DIR}/src
)

target_link_libraries(hbt_pair_frames
    PUBLIC
        hbt_boosts
        hbt_pair
)

target_compile_options(hbt_pair_frames
    PRIVATE
        -Wall
        -Wextra
        -Wpedantic
        -Werror
)

add_library(hbt_pair_origin_routing STATIC
    src/hbt/pair/pair_origin_routing.cpp
)

target_include_directories(hbt_pair_origin_routing
    PUBLIC
        ${PROJECT_SOURCE_DIR}/src
)

target_compile_options(hbt_pair_origin_routing
    PRIVATE
        -Wall
        -Wextra
        -Wpedantic
        -Werror
)

add_library(hbt_pair_slice_routing STATIC
    src/hbt/pair/pair_slice_routing.cpp
)

target_include_directories(hbt_pair_slice_routing
    PUBLIC
        ${PROJECT_SOURCE_DIR}/src
)

target_compile_options(hbt_pair_slice_routing
    PRIVATE
        -Wall
        -Wextra
        -Wpedantic
        -Werror
)

add_library(hbt_pair_slice_counts STATIC
    src/hbt/pair/pair_slice_count_accumulator.cpp
)

target_include_directories(hbt_pair_slice_counts
    PUBLIC
        ${PROJECT_SOURCE_DIR}/src
)

target_link_libraries(hbt_pair_slice_counts
    PUBLIC
        hbt_pair
        hbt_pair_origin_routing
        hbt_pair_slice_routing
)

target_compile_options(hbt_pair_slice_counts
    PRIVATE
        -Wall
        -Wextra
        -Wpedantic
        -Werror
)

add_library(hbt_rejected_pair_report STATIC
    src/hbt/reporting/rejected_pair_report.cpp
)

target_include_directories(hbt_rejected_pair_report
    PUBLIC
        ${PROJECT_SOURCE_DIR}/src
)

target_compile_options(hbt_rejected_pair_report
    PRIVATE
        -Wall
        -Wextra
        -Wpedantic
        -Werror
)

add_library(hbt_pair_processor STATIC
    src/hbt/pair/pair_processor.cpp
)

target_include_directories(hbt_pair_processor
    PUBLIC
        ${PROJECT_SOURCE_DIR}/src
)

target_link_libraries(hbt_pair_processor
    PUBLIC
        hbt_pair
        hbt_pair_frames
        hbt_pair_origin_routing
        hbt_pair_slice_counts
        hbt_rejected_pair_report
)

target_compile_options(hbt_pair_processor
    PRIVATE
        -Wall
        -Wextra
        -Wpedantic
        -Werror
)

add_library(hbt_rejected_particle_report STATIC
    src/hbt/reporting/rejected_particle_report.cpp
)

target_include_directories(hbt_rejected_particle_report
    PUBLIC
        ${PROJECT_SOURCE_DIR}/src
)

target_compile_options(hbt_rejected_particle_report
    PRIVATE
        -Wall
        -Wextra
        -Wpedantic
        -Werror
)

add_library(hbt_origin_selector STATIC
    src/hbt/event/origin_selector.cpp
)

target_include_directories(hbt_origin_selector
    PUBLIC
        ${PROJECT_SOURCE_DIR}/src
)

target_link_libraries(hbt_origin_selector
    PUBLIC
        hbt_config
)

target_compile_options(hbt_origin_selector
    PRIVATE
        -Wall
        -Wextra
        -Wpedantic
        -Werror
)

add_library(run_config STATIC
    src/config/run_config_loader.cpp
)

target_include_directories(run_config
    PUBLIC
        ${PROJECT_SOURCE_DIR}/src
)

target_link_libraries(run_config
    PRIVATE
        ${YAML_CPP_TARGET}
)

target_compile_options(run_config
    PRIVATE
        -Wall
        -Wextra
        -Wpedantic
        -Werror
)

add_library(hbt_config_gate STATIC
    src/app/hbt_config_gate.cpp
)

target_include_directories(hbt_config_gate
    PUBLIC
        ${PROJECT_SOURCE_DIR}/src
)

target_link_libraries(hbt_config_gate
    PUBLIC
        run_config
        hbt_config
)

target_compile_options(hbt_config_gate
    PRIVATE
        -Wall
        -Wextra
        -Wpedantic
        -Werror
)

add_library(hbt_startup STATIC
    src/hbt/startup/hbt_startup_builder.cpp
)

target_include_directories(hbt_startup
    PUBLIC
        ${PROJECT_SOURCE_DIR}/src
)

target_link_libraries(hbt_startup
    PUBLIC
        hbt_config
        hbt_selection
)

target_compile_options(hbt_startup
    PRIVATE
        -Wall
        -Wextra
        -Wpedantic
        -Werror
)

add_library(hbt_histograms STATIC
    src/hbt/histograms/product_fanout_plan.cpp
    src/hbt/histograms/raw_histograms.cpp
)

target_include_directories(hbt_histograms
    PUBLIC
        ${PROJECT_SOURCE_DIR}/src
)

target_link_libraries(hbt_histograms
    PUBLIC
        hbt_config
        hbt_startup
        hbt_pair_frames
        hbt_pair_origin_routing
        hbt_pair_slice_routing
)

target_compile_options(hbt_histograms
    PRIVATE
        -Wall
        -Wextra
        -Wpedantic
        -Werror
)

add_library(hbt_fit_statistics STATIC
    src/hbt/fits/fit_results.cpp
    src/hbt/fits/statistical_analysis.cpp
    src/hbt/fits/binned_models.cpp
)

target_include_directories(hbt_fit_statistics
    PUBLIC
        ${PROJECT_SOURCE_DIR}/src
)

target_link_libraries(hbt_fit_statistics
    PUBLIC
        hbt_config
)

target_compile_options(hbt_fit_statistics
    PRIVATE
        -Wall
        -Wextra
        -Wpedantic
        -Werror
)

add_library(hbt_minuit2_fits STATIC
    src/hbt/fits/minuit2_fitter.cpp
)

target_include_directories(hbt_minuit2_fits
    PUBLIC
        ${PROJECT_SOURCE_DIR}/src
)

target_link_libraries(hbt_minuit2_fits
    PUBLIC
        hbt_fit_statistics
        ROOT::Minuit2
)

target_compile_options(hbt_minuit2_fits
    PRIVATE
        -Wall
        -Wextra
        -Wpedantic
        -Werror
)

add_library(hbt_histogram_analysis STATIC
    src/hbt/fits/histogram_analysis.cpp
)

target_include_directories(hbt_histogram_analysis
    PUBLIC
        ${PROJECT_SOURCE_DIR}/src
)

target_link_libraries(hbt_histogram_analysis
    PUBLIC
        hbt_histograms
        hbt_minuit2_fits
)

target_compile_options(hbt_histogram_analysis
    PRIVATE
        -Wall
        -Wextra
        -Wpedantic
        -Werror
)

add_library(hbt_histograms_test_instrumented STATIC
    src/hbt/histograms/product_fanout_plan.cpp
    src/hbt/histograms/raw_histograms.cpp
)

target_include_directories(hbt_histograms_test_instrumented
    PUBLIC
        ${PROJECT_SOURCE_DIR}/src
)

target_link_libraries(hbt_histograms_test_instrumented
    PUBLIC
        hbt_config
        hbt_startup
        hbt_pair_frames
        hbt_pair_origin_routing
        hbt_pair_slice_routing
)

target_compile_definitions(hbt_histograms_test_instrumented
    PRIVATE
        HBT_RAW_HISTOGRAM_TEST_INSTRUMENTATION
)

target_compile_options(hbt_histograms_test_instrumented
    PRIVATE
        -Wall
        -Wextra
        -Wpedantic
        -Werror
)

add_library(analysis_runner STATIC
    src/app/analysis_runner.cpp
    src/app/console_progress_bar.cpp
)

target_include_directories(analysis_runner
    PUBLIC
        ${PROJECT_SOURCE_DIR}/src
)

target_link_libraries(analysis_runner
    PUBLIC
        run_config
        hbt_config_gate
        hbt_startup
        hbt_rejected_particle_report
        hbt_histograms
        hbt_histogram_analysis
    PRIVATE
        kinematics
        kinematics_validation
        afterburner_reader
        sampler_reader
        emission_point_resolver
        hbt_species
        hbt_selection
        hbt_particle_selector
        hbt_origin_selector
        hbt_pair_processor
        hbt_pair_slice_counts
        Threads::Threads
)

target_compile_options(analysis_runner
    PRIVATE
        -Wall
        -Wextra
        -Wpedantic
        -Werror
)

add_library(analysis_output_writer STATIC
    src/output/analysis_output_writer.cpp
    src/output/hbt_production_output.cpp
)

target_include_directories(analysis_output_writer
    PUBLIC
        ${PROJECT_SOURCE_DIR}/src
)

target_link_libraries(analysis_output_writer
    PUBLIC
        analysis_runner
    PRIVATE
        hbt_channels
        hbt_species
        hbt_rejected_pair_report
)

target_compile_options(analysis_output_writer
    PRIVATE
        -Wall
        -Wextra
        -Wpedantic
        -Werror
)

add_executable(smash-afterburner-analysis-modular
    src/main.cpp
)

target_link_libraries(smash-afterburner-analysis-modular
    PRIVATE
        analysis_output_writer
)

target_compile_options(smash-afterburner-analysis-modular
    PRIVATE
        -Wall
        -Wextra
        -Wpedantic
        -Werror
)

add_executable(species_test
    tests/unit/species_test.cpp
)

target_link_libraries(species_test
    PRIVATE
        hbt_species
)

target_compile_options(species_test
    PRIVATE
        -Wall
        -Wextra
        -Wpedantic
        -Werror
)

add_test(
    NAME species_test
    COMMAND species_test
)

add_executable(naming_catalog_test
    tests/unit/naming_catalog_test.cpp
)

target_link_libraries(naming_catalog_test
    PRIVATE
        hbt_channels
)

target_compile_options(naming_catalog_test
    PRIVATE
        -Wall
        -Wextra
        -Wpedantic
        -Werror
)

add_test(
    NAME naming_catalog_test
    COMMAND naming_catalog_test
)

add_executable(four_vector_test
    tests/unit/four_vector_test.cpp
)

target_include_directories(four_vector_test
    PRIVATE
        ${PROJECT_SOURCE_DIR}/src
)

target_compile_options(four_vector_test
    PRIVATE
        -Wall
        -Wextra
        -Wpedantic
        -Werror
)

add_test(
    NAME four_vector_test
    COMMAND four_vector_test
)

add_executable(kinematics_test
    tests/unit/kinematics_test.cpp
)

target_link_libraries(kinematics_test
    PRIVATE
        kinematics
)

target_compile_options(kinematics_test
    PRIVATE
        -Wall
        -Wextra
        -Wpedantic
        -Werror
)

add_test(
    NAME kinematics_test
    COMMAND kinematics_test
)

add_executable(kinematics_validation_test
    tests/unit/kinematics_validation_test.cpp
)

target_link_libraries(kinematics_validation_test
    PRIVATE
        kinematics_validation
)

target_compile_options(kinematics_validation_test
    PRIVATE
        -Wall
        -Wextra
        -Wpedantic
        -Werror
)

add_test(
    NAME kinematics_validation_test
    COMMAND kinematics_validation_test
)

add_executable(afterburner_reader_test
    tests/unit/afterburner_reader_test.cpp
)

target_link_libraries(afterburner_reader_test
    PRIVATE
        afterburner_reader
)

target_compile_options(afterburner_reader_test
    PRIVATE
        -Wall
        -Wextra
        -Wpedantic
        -Werror
)

add_test(
    NAME afterburner_reader_test
    COMMAND afterburner_reader_test
)

add_executable(sampler_reader_test
    tests/unit/sampler_reader_test.cpp
)

target_link_libraries(sampler_reader_test
    PRIVATE
        sampler_reader
)

target_compile_options(sampler_reader_test
    PRIVATE
        -Wall
        -Wextra
        -Wpedantic
        -Werror
)

add_test(
    NAME sampler_reader_test
    COMMAND sampler_reader_test
)

add_executable(emission_point_resolver_test
    tests/unit/emission_point_resolver_test.cpp
)

target_link_libraries(emission_point_resolver_test
    PRIVATE
        emission_point_resolver
)

target_compile_options(emission_point_resolver_test
    PRIVATE
        -Wall
        -Wextra
        -Wpedantic
        -Werror
)

add_test(
    NAME emission_point_resolver_test
    COMMAND emission_point_resolver_test
)

add_executable(particle_selector_test
    tests/unit/particle_selector_test.cpp
)

target_link_libraries(particle_selector_test
    PRIVATE
        hbt_particle_selector
)

target_compile_options(particle_selector_test
    PRIVATE
        -Wall
        -Wextra
        -Wpedantic
        -Werror
)

add_test(
    NAME particle_selector_test
    COMMAND particle_selector_test
)

add_executable(event_buffers_test
    tests/unit/event_buffers_test.cpp
)

target_link_libraries(event_buffers_test
    PRIVATE
        hbt_event_buffers
)

target_compile_options(event_buffers_test
    PRIVATE
        -Wall
        -Wextra
        -Wpedantic
        -Werror
)

add_test(
    NAME event_buffers_test
    COMMAND event_buffers_test
)

add_executable(subevent_shuffle_test
    tests/unit/subevent_shuffle_test.cpp
)

target_link_libraries(subevent_shuffle_test
    PRIVATE
        hbt_event_buffers
)

target_compile_options(subevent_shuffle_test
    PRIVATE
        -Wall
        -Wextra
        -Wpedantic
        -Werror
)

add_test(
    NAME subevent_shuffle_test
    COMMAND subevent_shuffle_test
)

add_executable(pair_iterator_test
    tests/unit/pair_iterator_test.cpp
)

target_link_libraries(pair_iterator_test
    PRIVATE
        hbt_pair
)

target_compile_options(pair_iterator_test
    PRIVATE
        -Wall
        -Wextra
        -Wpedantic
        -Werror
)

add_test(
    NAME pair_iterator_test
    COMMAND pair_iterator_test
)

add_executable(pair_count_accumulator_test
    tests/unit/pair_count_accumulator_test.cpp
)

target_link_libraries(pair_count_accumulator_test
    PRIVATE
        hbt_pair
)

target_compile_options(pair_count_accumulator_test
    PRIVATE
        -Wall
        -Wextra
        -Wpedantic
        -Werror
)

add_test(
    NAME pair_count_accumulator_test
    COMMAND pair_count_accumulator_test
)

add_executable(pair_kinematics_test
    tests/unit/pair_kinematics_test.cpp
)

target_link_libraries(pair_kinematics_test
    PRIVATE
        hbt_pair
)

target_compile_options(pair_kinematics_test
    PRIVATE
        -Wall
        -Wextra
        -Wpedantic
        -Werror
)

add_test(
    NAME pair_kinematics_test
    COMMAND pair_kinematics_test
)

add_executable(lab_to_lcms_test
    tests/unit/lab_to_lcms_test.cpp
)

target_link_libraries(lab_to_lcms_test
    PRIVATE
        hbt_boosts
)

target_compile_options(lab_to_lcms_test
    PRIVATE
        -Wall
        -Wextra
        -Wpedantic
        -Werror
)

add_test(
    NAME lab_to_lcms_test
    COMMAND lab_to_lcms_test
)

add_executable(lcms_to_prf_test
    tests/unit/lcms_to_prf_test.cpp
)

target_link_libraries(lcms_to_prf_test
    PRIVATE
        hbt_boosts
)

target_compile_options(lcms_to_prf_test
    PRIVATE
        -Wall
        -Wextra
        -Wpedantic
        -Werror
)

add_test(
    NAME lcms_to_prf_test
    COMMAND lcms_to_prf_test
)

add_executable(pair_frame_observables_test
    tests/unit/pair_frame_observables_test.cpp
)

target_link_libraries(pair_frame_observables_test
    PRIVATE
        hbt_pair_frames
)

target_compile_options(pair_frame_observables_test
    PRIVATE
        -Wall
        -Wextra
        -Wpedantic
        -Werror
)

add_test(
    NAME pair_frame_observables_test
    COMMAND pair_frame_observables_test
)

add_executable(pair_kinematics_validation_test
    tests/unit/pair_kinematics_validation_test.cpp
)

target_link_libraries(pair_kinematics_validation_test
    PRIVATE
        hbt_pair
)

target_compile_options(pair_kinematics_validation_test
    PRIVATE
        -Wall
        -Wextra
        -Wpedantic
        -Werror
)

add_test(
    NAME pair_kinematics_validation_test
    COMMAND pair_kinematics_validation_test
)

add_executable(rejected_pair_report_test
    tests/unit/rejected_pair_report_test.cpp
)

target_link_libraries(rejected_pair_report_test
    PRIVATE
        hbt_rejected_pair_report
)

target_compile_options(rejected_pair_report_test
    PRIVATE
        -Wall
        -Wextra
        -Wpedantic
        -Werror
)

add_test(
    NAME rejected_pair_report_test
    COMMAND rejected_pair_report_test
)

add_executable(pair_processor_test
    tests/unit/pair_processor_test.cpp
)

target_link_libraries(pair_processor_test
    PRIVATE
        hbt_pair_processor
)

target_compile_options(pair_processor_test
    PRIVATE
        -Wall
        -Wextra
        -Wpedantic
        -Werror
)

add_test(
    NAME pair_processor_test
    COMMAND pair_processor_test
)

# Compile pair_processor.cpp directly so this integration test can replace
# the frame calculator with a call-counting link double.
# Production behavior is unchanged.
add_executable(pair_frame_gate_integration_test
    tests/integration/pair_frame_gate_integration_test.cpp
    src/hbt/pair/pair_processor.cpp
    src/hbt/pair/pair_frame_observables_validation.cpp
)

target_include_directories(pair_frame_gate_integration_test
    PRIVATE
        ${PROJECT_SOURCE_DIR}/src
)

target_link_libraries(pair_frame_gate_integration_test
    PRIVATE
        hbt_pair_slice_counts
        hbt_rejected_pair_report
)

target_compile_options(pair_frame_gate_integration_test
    PRIVATE
        -Wall
        -Wextra
        -Wpedantic
        -Werror
)

add_test(
    NAME pair_frame_gate_integration_test
    COMMAND pair_frame_gate_integration_test
)

add_executable(pair_origin_routing_test
    tests/unit/pair_origin_routing_test.cpp
)

target_link_libraries(pair_origin_routing_test
    PRIVATE
        hbt_pair_origin_routing
)

target_compile_options(pair_origin_routing_test
    PRIVATE
        -Wall
        -Wextra
        -Wpedantic
        -Werror
)

add_test(
    NAME pair_origin_routing_test
    COMMAND pair_origin_routing_test
)

add_executable(pair_slice_routing_test
    tests/unit/pair_slice_routing_test.cpp
)

target_link_libraries(pair_slice_routing_test
    PRIVATE
        hbt_pair_slice_routing
)

target_compile_options(pair_slice_routing_test
    PRIVATE
        -Wall
        -Wextra
        -Wpedantic
        -Werror
)

add_test(
    NAME pair_slice_routing_test
    COMMAND pair_slice_routing_test
)

add_executable(pair_slice_count_accumulator_test
    tests/unit/pair_slice_count_accumulator_test.cpp
)

target_link_libraries(pair_slice_count_accumulator_test
    PRIVATE
        hbt_pair_slice_counts
)

target_compile_options(pair_slice_count_accumulator_test
    PRIVATE
        -Wall
        -Wextra
        -Wpedantic
        -Werror
)

add_test(
    NAME pair_slice_count_accumulator_test
    COMMAND pair_slice_count_accumulator_test
)

add_executable(rejected_particle_report_test
    tests/unit/rejected_particle_report_test.cpp
)

target_link_libraries(rejected_particle_report_test
    PRIVATE
        hbt_rejected_particle_report
)

target_compile_options(rejected_particle_report_test
    PRIVATE
        -Wall
        -Wextra
        -Wpedantic
        -Werror
)

add_test(
    NAME rejected_particle_report_test
    COMMAND rejected_particle_report_test
)

add_executable(origin_selector_test
    tests/unit/origin_selector_test.cpp
)

target_link_libraries(origin_selector_test
    PRIVATE
        hbt_origin_selector
)

target_compile_options(origin_selector_test
    PRIVATE
        -Wall
        -Wextra
        -Wpedantic
        -Werror
)

add_test(
    NAME origin_selector_test
    COMMAND origin_selector_test
)

add_executable(channel_catalog_test
    tests/unit/channel_catalog_test.cpp
)

target_link_libraries(channel_catalog_test
    PRIVATE
        hbt_channels
)

target_compile_options(channel_catalog_test
    PRIVATE
        -Wall
        -Wextra
        -Wpedantic
        -Werror
)

add_test(
    NAME channel_catalog_test
    COMMAND channel_catalog_test
)

add_executable(required_channels_test
    tests/unit/required_channels_test.cpp
)

target_link_libraries(required_channels_test
    PRIVATE
        hbt_selection
)

target_compile_options(required_channels_test
    PRIVATE
        -Wall
        -Wextra
        -Wpedantic
        -Werror
)

add_test(
    NAME required_channels_test
    COMMAND required_channels_test
)

add_executable(required_species_test
    tests/unit/required_species_test.cpp
)

target_link_libraries(required_species_test
    PRIVATE
        hbt_selection
)

target_compile_options(required_species_test
    PRIVATE
        -Wall
        -Wextra
        -Wpedantic
        -Werror
)

add_test(
    NAME required_species_test
    COMMAND required_species_test
)

add_executable(species_requirement_test
    tests/unit/species_requirement_test.cpp
)

target_link_libraries(species_requirement_test
    PRIVATE
        hbt_selection
)

target_compile_options(species_requirement_test
    PRIVATE
        -Wall
        -Wextra
        -Wpedantic
        -Werror
)

add_test(
    NAME species_requirement_test
    COMMAND species_requirement_test
)

add_executable(primitive_channel_token_test
    tests/unit/primitive_channel_token_test.cpp
)

target_link_libraries(primitive_channel_token_test
    PRIVATE
        hbt_config
)

target_compile_options(primitive_channel_token_test
    PRIVATE
        -Wall
        -Wextra
        -Wpedantic
        -Werror
)

add_test(
    NAME primitive_channel_token_test
    COMMAND primitive_channel_token_test
)

add_executable(analysis_product_parser_test
    tests/unit/analysis_product_parser_test.cpp
)

target_link_libraries(analysis_product_parser_test
    PRIVATE
        hbt_config
)

target_compile_options(analysis_product_parser_test
    PRIVATE
        -Wall
        -Wextra
        -Wpedantic
        -Werror
)

add_test(
    NAME analysis_product_parser_test
    COMMAND analysis_product_parser_test
)

add_executable(hbt_enabled_channels_parser_test
    tests/unit/hbt_enabled_channels_parser_test.cpp
)

target_link_libraries(hbt_enabled_channels_parser_test
    PRIVATE
        hbt_config
)

target_compile_options(hbt_enabled_channels_parser_test
    PRIVATE
        -Wall
        -Wextra
        -Wpedantic
        -Werror
)

add_test(
    NAME hbt_enabled_channels_parser_test
    COMMAND hbt_enabled_channels_parser_test
)

add_executable(hbt_config_loader_test
    tests/unit/hbt_config_loader_test.cpp
)

target_link_libraries(hbt_config_loader_test
    PRIVATE
        hbt_config
)

target_compile_options(hbt_config_loader_test
    PRIVATE
        -Wall
        -Wextra
        -Wpedantic
        -Werror
)

add_test(
    NAME hbt_config_loader_test
    COMMAND hbt_config_loader_test
)

add_executable(run_config_loader_test
    tests/unit/run_config_loader_test.cpp
)

target_link_libraries(run_config_loader_test
    PRIVATE
        run_config
)

target_compile_options(run_config_loader_test
    PRIVATE
        -Wall
        -Wextra
        -Wpedantic
        -Werror
)

add_test(
    NAME run_config_loader_test
    COMMAND run_config_loader_test
)

add_executable(hbt_config_gate_test
    tests/unit/hbt_config_gate_test.cpp
)

target_link_libraries(hbt_config_gate_test
    PRIVATE
        hbt_config_gate
)

target_compile_options(hbt_config_gate_test
    PRIVATE
        -Wall
        -Wextra
        -Wpedantic
        -Werror
)

add_test(
    NAME hbt_config_gate_test
    COMMAND hbt_config_gate_test
)

add_executable(hbt_startup_builder_test
    tests/unit/hbt_startup_builder_test.cpp
)

target_link_libraries(hbt_startup_builder_test
    PRIVATE
        hbt_startup
)

target_compile_options(hbt_startup_builder_test
    PRIVATE
        -Wall
        -Wextra
        -Wpedantic
        -Werror
)

add_test(
    NAME hbt_startup_builder_test
    COMMAND hbt_startup_builder_test
)

add_executable(event_preparation_integration_test
    tests/integration/event_preparation_integration_test.cpp
)

target_link_libraries(event_preparation_integration_test
    PRIVATE
        analysis_output_writer
)

target_compile_options(event_preparation_integration_test
    PRIVATE
        -Wall
        -Wextra
        -Wpedantic
        -Werror
)

add_test(
    NAME event_preparation_integration_test
    COMMAND event_preparation_integration_test
)

add_executable(event_preparation_real_fixture_test
    tests/regression/event_preparation_real_fixture_test.cpp
)

target_link_libraries(event_preparation_real_fixture_test
    PRIVATE
        analysis_runner
)

target_compile_options(event_preparation_real_fixture_test
    PRIVATE
        -Wall
        -Wextra
        -Wpedantic
        -Werror
)

if(HBT_REAL_REGRESSION_MAIN_CONFIG)
    get_filename_component(
        HBT_REAL_REGRESSION_MAIN_CONFIG_ABSOLUTE
        "${HBT_REAL_REGRESSION_MAIN_CONFIG}"
        ABSOLUTE
        BASE_DIR "${PROJECT_SOURCE_DIR}"
    )

    add_test(
        NAME event_preparation_real_fixture_test
        COMMAND
            event_preparation_real_fixture_test
            "${HBT_REAL_REGRESSION_MAIN_CONFIG_ABSOLUTE}"
    )
endif()

add_executable(analysis_runner_test
    tests/unit/analysis_runner_test.cpp
)

target_link_libraries(analysis_runner_test
    PRIVATE
        analysis_runner
)

target_compile_options(analysis_runner_test
    PRIVATE
        -Wall
        -Wextra
        -Wpedantic
        -Werror
)

add_test(
    NAME analysis_runner_test
    COMMAND analysis_runner_test
)

add_executable(analysis_output_writer_test
    tests/unit/analysis_output_writer_test.cpp
)

target_link_libraries(analysis_output_writer_test
    PRIVATE
        analysis_output_writer
)

target_compile_options(analysis_output_writer_test
    PRIVATE
        -Wall
        -Wextra
        -Wpedantic
        -Werror
)

add_test(
    NAME analysis_output_writer_test
    COMMAND analysis_output_writer_test
)

add_test(
    NAME analysis_executable_startup_test
    COMMAND
        smash-afterburner-analysis-modular
        ${PROJECT_SOURCE_DIR}/tests/fixtures/main_hbt_disabled.yaml
)

add_test(
    NAME analysis_executable_requires_run_config_path_test
    COMMAND smash-afterburner-analysis-modular
)

set_tests_properties(
    analysis_executable_requires_run_config_path_test
    PROPERTIES
        WILL_FAIL TRUE
)

add_executable(pair_frame_observables_validation_test
    tests/unit/pair_frame_observables_validation_test.cpp
)

target_link_libraries(pair_frame_observables_validation_test
    PRIVATE
        hbt_pair_frames
)

target_compile_options(pair_frame_observables_validation_test
    PRIVATE
        -Wall
        -Wextra
        -Wpedantic
        -Werror
)

add_test(
    NAME pair_frame_observables_validation_test
    COMMAND pair_frame_observables_validation_test
)

add_executable(product_fanout_plan_test
    tests/unit/product_fanout_plan_test.cpp
)

target_link_libraries(product_fanout_plan_test
    PRIVATE
        hbt_histograms
)

target_compile_options(product_fanout_plan_test
    PRIVATE
        -Wall
        -Wextra
        -Wpedantic
        -Werror
)

add_test(
    NAME product_fanout_plan_test
    COMMAND product_fanout_plan_test
)

add_executable(raw_histograms_test
    tests/unit/raw_histograms_test.cpp
)

target_link_libraries(raw_histograms_test
    PRIVATE
        hbt_histograms_test_instrumented
)

target_compile_options(raw_histograms_test
    PRIVATE
        -Wall
        -Wextra
        -Wpedantic
        -Werror
)

add_test(
    NAME raw_histograms_test
    COMMAND raw_histograms_test
)

add_executable(raw_histogram_accumulation_integration_test
    tests/integration/raw_histogram_accumulation_integration_test.cpp
)

target_link_libraries(raw_histogram_accumulation_integration_test
    PRIVATE
        hbt_histograms
        hbt_pair_processor
)

target_compile_options(raw_histogram_accumulation_integration_test
    PRIVATE
        -Wall
        -Wextra
        -Wpedantic
        -Werror
)

add_test(
    NAME raw_histogram_accumulation_integration_test
    COMMAND raw_histogram_accumulation_integration_test
)

add_executable(statistical_analysis_test
    tests/unit/statistical_analysis_test.cpp
)

target_link_libraries(statistical_analysis_test
    PRIVATE
        hbt_fit_statistics
)

target_compile_options(statistical_analysis_test
    PRIVATE
        -Wall
        -Wextra
        -Wpedantic
        -Werror
)

add_test(
    NAME statistical_analysis_test
    COMMAND statistical_analysis_test
)

add_executable(binned_models_test
    tests/unit/binned_models_test.cpp
)

target_link_libraries(binned_models_test
    PRIVATE
        hbt_fit_statistics
)

target_compile_options(binned_models_test
    PRIVATE
        -Wall
        -Wextra
        -Wpedantic
        -Werror
)

add_test(
    NAME binned_models_test
    COMMAND binned_models_test
)

add_executable(minuit2_fitter_test
    tests/unit/minuit2_fitter_test.cpp
)

target_link_libraries(minuit2_fitter_test
    PRIVATE
        hbt_minuit2_fits
)

target_compile_options(minuit2_fitter_test
    PRIVATE
        -Wall
        -Wextra
        -Wpedantic
        -Werror
)

add_test(
    NAME minuit2_fitter_test
    COMMAND minuit2_fitter_test
)

add_executable(hbt_production_output_integration_test
    tests/integration/hbt_production_output_integration_test.cpp
)

target_link_libraries(hbt_production_output_integration_test
    PRIVATE
        analysis_output_writer
)

target_compile_options(hbt_production_output_integration_test
    PRIVATE
        -Wall
        -Wextra
        -Wpedantic
        -Werror
)

add_test(
    NAME hbt_production_output_integration_test
    COMMAND hbt_production_output_integration_test
)
