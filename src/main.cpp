/**
 * @file main.cpp
 * @brief Command-line entry point for the modular analysis executable.
 *
 * The executable entry point validates only its command-line contract,
 * delegates orchestration to app::AnalysisRunner, and delegates result
 * serialization to the output writer.
 */

#include "app/analysis_runner.h"
#include "app/console_progress_bar.h"
#include "output/analysis_output_writer.h"

#include <cstdlib>
#include <exception>
#include <filesystem>
#include <iostream>

/**
 * @brief Start the modular analysis application.
 *
 * Exactly one argument is required: the path to the global run-configuration
 * YAML file. Configuration loading, analysis orchestration, and output
 * serialization are delegated to dedicated components.
 *
 * @param argc Number of command-line arguments including the executable name.
 * @param argv Command-line argument array.
 * @return EXIT_SUCCESS when the current modular run completes, otherwise
 *         EXIT_FAILURE.
 */
int main(
    int argc,
    char* argv[]
) {
    if (argc != 2) {
        std::cerr
            << "usage: smash-afterburner-analysis-modular <main.yaml>\n";
        return EXIT_FAILURE;
    }

    app::ConsoleProgressBar progress;

    try {
        const app::AnalysisRunner runner{
            std::filesystem::path(argv[1])
        };

        const app::AnalysisRunSummary result = runner.run(&progress);
        output::write_analysis_output(result, std::cout);

        if (result.startup.run_config.hbt_enabled) {
            progress.begin_output();
        }

        output::write_production_output(result);

        if (result.startup.run_config.hbt_enabled) {
            progress.finish();
        }
    } catch (const std::exception& exception) {
        progress.fail();
        std::cerr
            << "analysis run failed: "
            << exception.what()
            << '\n';
        return EXIT_FAILURE;
    }

    return EXIT_SUCCESS;
}
