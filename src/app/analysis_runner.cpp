/**
 * @file analysis_runner.cpp
 * @brief Application-level startup, preparation, and pair-processing
 *        orchestration.
 */

#include "app/analysis_runner.h"

#include "app/hbt_config_gate.h"
#include "common/kinematics.h"
#include "common/kinematics_validation.h"
#include "config/run_config_loader.h"
#include "hbt/event/event_buffers.h"
#include "hbt/event/origin_selector.h"
#include "hbt/event/particle_selector.h"
#include "hbt/event/subevent_shuffle.h"
#include "hbt/fits/histogram_analysis.h"
#include "hbt/histograms/raw_histograms.h"
#include "hbt/pair/pair_count_accumulator.h"
#include "hbt/pair/pair_processor.h"
#include "hbt/pair/pair_slice_count_accumulator.h"
#include "hbt/selection/species_requirement.h"
#include "hbt/species/species.h"
#include "hbt/startup/hbt_startup_builder.h"
#include "input/afterburner_reader.h"
#include "input/emission_point_resolver.h"
#include "input/sampler_reader.h"

#include <algorithm>
#include <atomic>
#include <condition_variable>
#include <exception>
#include <limits>
#include <mutex>
#include <optional>
#include <stdexcept>
#include <string>
#include <thread>
#include <utility>
#include <vector>

namespace app {
namespace {

/**
 * @brief Exception marking one explicitly recoverable subevent-local failure.
 *
 * Only failures that leave the current Afterburner stream recoverable are
 * converted to this type. Internal invariants, overflows, pair-processing
 * failures, and unclassified exceptions remain fatal.
 */
class RecoverableSubeventFailure final : public std::runtime_error {
public:
    using std::runtime_error::runtime_error;
};

/**
 * @brief Exception marking one explicitly recoverable outer-event failure.
 *
 * The outer event is transactionally discarded before this exception is
 * converted to SkippedDueToEventFailure. Failures outside the event-input
 * boundary are never converted to this type.
 */
class RecoverableEventFailure final : public std::runtime_error {
public:
    using std::runtime_error::runtime_error;
};

/**
 * @brief Execute one event-input operation and classify runtime input errors.
 * @tparam Operation Callable returning the wrapped operation result.
 * @param operation Event-local input operation to invoke.
 * @return The wrapped operation result, including void where applicable.
 * @throws RecoverableEventFailure If the input operation reports runtime
 *         failure.
 *
 * The wrapper is used only around Sampler/Afterburner input operations and
 * event-local cardinality validation. It is deliberately not used around
 * scientific processing, histogram accumulation, arithmetic, or invariants.
 */
template <typename Operation>
decltype(auto) event_input_operation(Operation&& operation) {
    try {
        return std::forward<Operation>(operation)();
    } catch (const std::runtime_error& error) {
        throw RecoverableEventFailure(error.what());
    }
}

/**
 * @brief Build the directory path for one one-based outer event.
 * @param events_path Directory containing all outer-event directories.
 * @param outer_event_number One-based outer-event number.
 * @return Path to the selected outer-event directory.
 */
std::filesystem::path outer_event_directory(
    const std::filesystem::path& events_path,
    std::size_t outer_event_number
) {
    return events_path / std::to_string(outer_event_number);
}

/**
 * @brief Build the Afterburner particle-list path for one outer event.
 * @param outer_event_dir Directory of the selected outer event.
 * @return Path to its Afterburner particle list.
 */
std::filesystem::path afterburner_path(
    const std::filesystem::path& outer_event_dir
) {
    return outer_event_dir / "Afterburner" / "particle_lists.oscar";
}

/**
 * @brief Build the Sampler particle-list path for one outer event.
 * @param outer_event_dir Directory of the selected outer event.
 * @return Path to its Sampler particle list.
 */
std::filesystem::path sampler_path(
    const std::filesystem::path& outer_event_dir
) {
    return outer_event_dir / "Sampler" / "particle_lists.oscar";
}

/**
 * @brief Build zeroed per-species counters for the required species.
 * @param required_species Startup-resolved species required by HBT.
 * @return Counters in the same stable order as required_species.
 */
std::vector<HBTSpeciesPreparationCounts> make_species_counts(
    const std::vector<hbt::SpeciesId>& required_species
) {
    std::vector<HBTSpeciesPreparationCounts> counts;
    counts.reserve(required_species.size());

    for (const hbt::SpeciesId species : required_species) {
        counts.push_back({species, 0U, 0U, 0U, 0U});
    }

    return counts;
}

/**
 * @brief Find mutable preparation counters for one required species.
 * @param counts Per-required-species preparation counters.
 * @param species Required species whose counters are requested.
 * @return Mutable counters for the species.
 * @throws std::logic_error If the required species has no counter entry.
 */
HBTSpeciesPreparationCounts& species_counts_for(
    std::vector<HBTSpeciesPreparationCounts>& counts,
    hbt::SpeciesId species
) {
    for (HBTSpeciesPreparationCounts& entry : counts) {
        if (entry.species == species) {
            return entry;
        }
    }

    throw std::logic_error(
        "AnalysisRunner: required species has no preparation counter"
    );
}

/**
 * @brief Increment the counter for one emission-position source.
 * @param counts Aggregate emission-position source counters.
 * @param source Source reported by EmissionPointResolver.
 * @throws std::invalid_argument If source is not a valid enum value.
 */
void count_emission_source(
    HBTEmissionPointCounts& counts,
    input::EmissionPointSource source
) {
    switch (source) {
        case input::EmissionPointSource::Sampler:
            ++counts.sampler;
            return;

        case input::EmissionPointSource::Propagation:
            ++counts.propagation;
            return;

        case input::EmissionPointSource::Afterburner:
            ++counts.afterburner;
            return;
    }

    throw std::invalid_argument(
        "AnalysisRunner: invalid EmissionPointSource"
    );
}

/**
 * @brief Count one accepted particle and all of its nested origin memberships.
 * @param counts Per-required-species preparation counters.
 * @param species Accepted particle species.
 * @param origin Nested origin memberships of the accepted particle.
 */
void count_accepted_particle(
    std::vector<HBTSpeciesPreparationCounts>& counts,
    hbt::SpeciesId species,
    const hbt::OriginFlags& origin
) {
    HBTSpeciesPreparationCounts& entry =
        species_counts_for(counts, species);

    ++entry.accepted;

    if (origin.primordial) {
        ++entry.primordial;
    }

    if (origin.primordial_rescattering) {
        ++entry.primordial_rescattering;
    }

    if (origin.primordial_rescattering_decay) {
        ++entry.primordial_rescattering_decay;
    }
}

/**
 * @brief Throw a diagnostic for one mandatory Sampler lookup miss.
 * @param outer_event_number One-based outer-event number.
 * @param subevent_id Current Afterburner subevent identifier.
 * @param particle Raw particle whose mandatory lookup was absent.
 * @throws RecoverableSubeventFailure Always.
 */
[[noreturn]] void throw_missing_sampler_match(
    std::size_t outer_event_number,
    int subevent_id,
    const input::AfterburnerParticleRecord& particle
) {
    throw RecoverableSubeventFailure(
        "AnalysisRunner: missing mandatory Sampler match for outer event " +
        std::to_string(outer_event_number) +
        ", subevent " + std::to_string(subevent_id) +
        ", ID " + std::to_string(particle.id) +
        ", PDG " + std::to_string(particle.pdg)
    );
}

/**
 * @brief Map one acceptance-stage numerical reason to report terminology.
 * @param reason Numerical rejection reason returned by particle acceptance.
 * @return Matching stable particle-report rejection reason.
 * @throws std::invalid_argument If reason is not a valid enum value.
 */
hbt::ParticleRejectionReason report_reason(
    hbt::ParticleAcceptanceNumericalReason reason
) {
    switch (reason) {
        case hbt::ParticleAcceptanceNumericalReason::NonFiniteMomentum:
            return hbt::ParticleRejectionReason::NonFiniteMomentum;

        case hbt::ParticleAcceptanceNumericalReason::NonPositiveEnergy:
            return hbt::ParticleRejectionReason::NonPositiveEnergy;

        case hbt::ParticleAcceptanceNumericalReason::
                NonFiniteTransverseMomentum:
            return hbt::ParticleRejectionReason::NonFiniteTransverseMomentum;

        case hbt::ParticleAcceptanceNumericalReason::InvalidRapidityInput:
            return hbt::ParticleRejectionReason::InvalidRapidityInput;

        case hbt::ParticleAcceptanceNumericalReason::NonFiniteRapidity:
            return hbt::ParticleRejectionReason::NonFiniteRapidity;

        case hbt::ParticleAcceptanceNumericalReason::
                InvalidPseudorapidityInput:
            return hbt::ParticleRejectionReason::InvalidPseudorapidityInput;

        case hbt::ParticleAcceptanceNumericalReason::
                NonFinitePseudorapidity:
            return hbt::ParticleRejectionReason::NonFinitePseudorapidity;
    }

    throw std::invalid_argument(
        "AnalysisRunner: invalid particle-acceptance numerical reason"
    );
}

/**
 * @brief Store one complete recoverable numerical particle rejection.
 * @param report Aggregate rejected-particle report to update.
 * @param outer_event_number One-based outer-event number.
 * @param subevent_id Current Afterburner subevent identifier.
 * @param particle Raw rejected Afterburner particle.
 * @param species Canonical identified HBT species.
 * @param reason Exact numerical rejection reason.
 * @param diagnostic_value Invalid calculated scalar when one exists.
 * @param diagnostic_position Invalid selected position when one exists.
 */
void record_numerical_rejection(
    hbt::RejectedParticleReport& report,
    std::size_t outer_event_number,
    int subevent_id,
    const input::AfterburnerParticleRecord& particle,
    hbt::SpeciesId species,
    hbt::ParticleRejectionReason reason,
    std::optional<double> diagnostic_value = std::nullopt,
    std::optional<common::FourVector> diagnostic_position = std::nullopt
) {
    report.add({
        outer_event_number,
        subevent_id,
        particle.id,
        particle.pdg,
        particle.charge,
        species,
        particle.momentum,
        particle.position,
        particle.mass,
        particle.ncoll,
        particle.time_last_coll,
        reason,
        diagnostic_value,
        diagnostic_position
    });
}

/**
 * @brief Map a non-finite emission branch to its rejection-report reason.
 * @param source Emission-position branch that produced the invalid candidate.
 * @return Stable rejection reason for that source.
 * @throws std::invalid_argument If source is not a valid enum value.
 */
hbt::ParticleRejectionReason emission_rejection_reason(
    input::EmissionPointSource source
) {
    switch (source) {
        case input::EmissionPointSource::Sampler:
            return hbt::ParticleRejectionReason::
                NonFiniteSamplerEmissionPosition;

        case input::EmissionPointSource::Propagation:
            return hbt::ParticleRejectionReason::
                NonFinitePropagationEmissionPosition;

        case input::EmissionPointSource::Afterburner:
            return hbt::ParticleRejectionReason::
                NonFiniteAfterburnerEmissionPosition;
    }

    throw std::invalid_argument(
        "AnalysisRunner: invalid emission source for numerical rejection"
    );
}

/**
 * @brief Process one raw particle through the ordered HBT preparation pipeline.
 * @param outer_event_number One-based outer-event number for diagnostics.
 * @param subevent_id Current Afterburner subevent identifier.
 * @param particle Raw Afterburner particle record.
 * @param config Resolved HBT scientific configuration.
 * @param startup Resolved HBT startup requirements.
 * @param sampler_reader Sampler index for the same outer event.
 * @param accepted_particles Current-subevent accepted-particle staging.
 *        Accepted Particle values remain here until the subevent is complete,
 *        when orchestration shuffles the full staging vector before species
 *        grouping and pair formation.
 * @param summary Aggregate HBT preparation diagnostics to update.
 * @return `true` when a final Particle was accepted into staging.
 * @throws RecoverableSubeventFailure If a mandatory Sampler match is absent.
 */
bool prepare_particle(
    std::size_t outer_event_number,
    int subevent_id,
    const input::AfterburnerParticleRecord& particle,
    const hbt::HBTConfig& config,
    const hbt::HBTStartupState& startup,
    const input::SamplerReader& sampler_reader,
    std::vector<hbt::Particle>& accepted_particles,
    HBTEventPreparationSummary& summary
) {
    ++summary.raw_particles;

    const std::optional<hbt::SpeciesId> species =
        hbt::identify_species(particle.charge, particle.pdg);

    if (!species.has_value()) {
        ++summary.unsupported_species;
        return false;
    }

    if (!hbt::is_species_required(
            species.value(),
            startup.required_species)) {
        ++summary.unrequired_species;
        return false;
    }

    const hbt::ParticleAcceptanceDecision acceptance =
        hbt::evaluate_particle_acceptance(
            particle.momentum,
            species.value(),
            config.particle_acceptance
        );

    if (acceptance.status == hbt::ParticleAcceptanceStatus::OutsideCuts) {
        ++summary.particle_acceptance_rejections;
        return false;
    }

    if (acceptance.status ==
        hbt::ParticleAcceptanceStatus::NumericalRejection) {
        if (!acceptance.numerical_reason.has_value()) {
            throw std::logic_error(
                "AnalysisRunner: numerical acceptance rejection has no reason"
            );
        }

        record_numerical_rejection(
            summary.numerical_rejections,
            outer_event_number,
            subevent_id,
            particle,
            species.value(),
            report_reason(acceptance.numerical_reason.value()),
            acceptance.diagnostic_value
        );
        return false;
    }

    const hbt::OriginFlags origin = hbt::classify_origin(
        particle.pdg_mother1,
        particle.pdg_mother2
    );

    if (!hbt::is_origin_eligible(origin, config.origin_mode)) {
        ++summary.origin_rejections;
        return false;
    }

    const double mass_squared =
        common::invariant_mass_squared(particle.momentum);

    if (!common::is_finite_kinematic_result(mass_squared)) {
        record_numerical_rejection(
            summary.numerical_rejections,
            outer_event_number,
            subevent_id,
            particle,
            species.value(),
            hbt::ParticleRejectionReason::NonFiniteInvariantMassSquared,
            mass_squared
        );
        return false;
    }

    if (!common::is_valid_invariant_mass_squared(mass_squared)) {
        record_numerical_rejection(
            summary.numerical_rejections,
            outer_event_number,
            subevent_id,
            particle,
            species.value(),
            hbt::ParticleRejectionReason::NonPositiveInvariantMassSquared,
            mass_squared
        );
        return false;
    }

    const double invariant_mass = common::invariant_mass(mass_squared);

    if (!common::is_finite_kinematic_result(invariant_mass)) {
        record_numerical_rejection(
            summary.numerical_rejections,
            outer_event_number,
            subevent_id,
            particle,
            species.value(),
            hbt::ParticleRejectionReason::NonFiniteInvariantMass,
            invariant_mass
        );
        return false;
    }

    const input::EmissionPointResolutionResult emission =
        input::resolve_emission_point(
            subevent_id,
            particle,
            origin.primordial,
            sampler_reader
        );

    if (emission.status ==
        input::EmissionPointResolutionStatus::MissingMandatorySampler) {
        throw_missing_sampler_match(
            outer_event_number,
            subevent_id,
            particle
        );
    }

    if (!emission.candidate.has_value()) {
        throw std::logic_error(
            "AnalysisRunner: emission result has no selected candidate"
        );
    }

    const input::ResolvedEmissionPoint& candidate =
        emission.candidate.value();

    if (emission.status ==
        input::EmissionPointResolutionStatus::NonFinitePosition) {
        record_numerical_rejection(
            summary.numerical_rejections,
            outer_event_number,
            subevent_id,
            particle,
            species.value(),
            emission_rejection_reason(candidate.source),
            std::nullopt,
            candidate.position
        );
        return false;
    }

    if (emission.status != input::EmissionPointResolutionStatus::Resolved) {
        throw std::logic_error(
            "AnalysisRunner: invalid emission resolution status"
        );
    }

    accepted_particles.push_back({
        species.value(),
        candidate.position,
        particle.momentum,
        invariant_mass,
        origin,
        particle.pdg,
        particle.charge
    });

    ++summary.accepted_particles;
    count_accepted_particle(summary.species, species.value(), origin);
    count_emission_source(summary.emission_points, candidate.source);
    return true;
}

/**
 * @brief Validate Sampler cardinality against the configured subevent count.
 * @param reader Loaded Sampler reader for one outer event.
 * @param expected Configured number of subevents per outer event.
 * @param outer_event_number One-based outer-event number for diagnostics.
 * @throws std::runtime_error If the counts differ.
 */
void require_sampler_subevent_count(
    const input::SamplerReader& reader,
    std::size_t expected,
    std::size_t outer_event_number
) {
    if (reader.subevent_count() != expected) {
        throw std::runtime_error(
            "AnalysisRunner: Sampler subevent count mismatch for outer event " +
            std::to_string(outer_event_number) +
            ": expected " + std::to_string(expected) +
            ", found " + std::to_string(reader.subevent_count())
        );
    }
}

/**
 * @brief Calculate the configured global subevent count without overflow.
 * @param run_config Resolved global run configuration.
 * @return number_of_events multiplied by number_of_subevents.
 * @throws std::overflow_error If the product cannot fit in std::size_t.
 */
std::size_t total_subevent_count(
    const config::RunConfig& run_config
) {
    if (
        run_config.number_of_subevents != 0U &&
        run_config.number_of_events >
            std::numeric_limits<std::size_t>::max() /
                run_config.number_of_subevents
    ) {
        throw std::overflow_error(
            "AnalysisRunner: total subevent count overflow"
        );
    }

    return run_config.number_of_events *
           run_config.number_of_subevents;
}

/**
 * @brief Build a zeroed preparation summary with stable species structure.
 * @param startup Resolved HBT startup requirements.
 * @return Empty summary ready for event-local or run-total accumulation.
 */
HBTEventPreparationSummary make_zero_preparation_summary(
    const hbt::HBTStartupState& startup
) {
    return {
        0U,
        0U,
        0U,
        0U,
        0U,
        0U,
        {},
        0U,
        0U,
        {0U, 0U, 0U},
        make_species_counts(startup.required_species),
        {},
        {},
        {},
        {}
    };
}

/**
 * @brief Build a zeroed pair summary with the configured stable layout.
 * @param config Resolved HBT scientific configuration.
 * @param startup Resolved HBT startup requirements.
 * @return Empty pair-processing summary.
 */
hbt::HBTPairProcessingSummary make_zero_pair_summary(
    const hbt::HBTConfig& config,
    const hbt::HBTStartupState& startup
) {
    const std::vector<hbt::PrimitiveChannelId>& required_channels =
        startup.required_primitive_channels;

    return {
        hbt::make_zero_pair_count_summary(required_channels),
        hbt::make_zero_pair_count_summary(required_channels),
        hbt::make_zero_pair_count_summary(required_channels),
        {
            config.origin_mode,
            hbt::make_zero_pair_count_summary(required_channels),
            hbt::make_zero_pair_count_summary(required_channels),
            hbt::make_zero_pair_count_summary(required_channels)
        },
        hbt::make_zero_pair_slice_count_summary(
            config.pair_slicing,
            config.origin_mode,
            required_channels
        ),
        {},
        {}
    };
}

/**
 * @brief Require one std::size_t reduction to fit without wraparound.
 * @param total Existing run-total value.
 * @param local Event-local value to add.
 * @param label Stable diagnostic label.
 * @throws std::overflow_error If total + local is not representable.
 */
void require_size_merge_fits(
    std::size_t total,
    std::size_t local,
    const char* label
) {
    if (local > std::numeric_limits<std::size_t>::max() - total) {
        throw std::overflow_error(
            std::string("AnalysisRunner: ") + label + " overflow"
        );
    }
}

/**
 * @brief Require one uint64_t reduction to fit without wraparound.
 * @param total Existing run-total value.
 * @param local Event-local value to add.
 * @param label Stable diagnostic label.
 * @throws std::overflow_error If total + local is not representable.
 */
void require_u64_merge_fits(
    std::uint64_t total,
    std::uint64_t local,
    const char* label
) {
    if (local > std::numeric_limits<std::uint64_t>::max() - total) {
        throw std::overflow_error(
            std::string("AnalysisRunner: ") + label + " overflow"
        );
    }
}

/**
 * @brief Increment one execution-status count with overflow checking.
 * @param count Mutable status counter.
 * @param label Stable diagnostic label.
 * @throws std::overflow_error If the counter is already maximal.
 */
void increment_status_count(
    std::uint64_t& count,
    const char* label
) {
    require_u64_merge_fits(count, 1U, label);
    ++count;
}

/**
 * @brief Validate all scalar preparation counters before one merge.
 * @param total Existing run-total preparation summary.
 * @param local Completed event-local preparation summary.
 * @throws std::overflow_error If any scalar reduction would overflow.
 */
void require_preparation_merge_fits(
    const HBTEventPreparationSummary& total,
    const HBTEventPreparationSummary& local
) {
    require_size_merge_fits(
        total.outer_events_processed,
        local.outer_events_processed,
        "outer-event preparation count"
    );
    require_size_merge_fits(
        total.subevents_processed,
        local.subevents_processed,
        "subevent preparation count"
    );
    require_size_merge_fits(
        total.raw_particles,
        local.raw_particles,
        "raw-particle count"
    );
    require_size_merge_fits(
        total.unsupported_species,
        local.unsupported_species,
        "unsupported-species count"
    );
    require_size_merge_fits(
        total.unrequired_species,
        local.unrequired_species,
        "unrequired-species count"
    );
    require_size_merge_fits(
        total.particle_acceptance_rejections,
        local.particle_acceptance_rejections,
        "particle-acceptance rejection count"
    );
    require_size_merge_fits(
        total.origin_rejections,
        local.origin_rejections,
        "origin-rejection count"
    );
    require_size_merge_fits(
        total.accepted_particles,
        local.accepted_particles,
        "accepted-particle count"
    );
    require_size_merge_fits(
        total.emission_points.sampler,
        local.emission_points.sampler,
        "Sampler emission count"
    );
    require_size_merge_fits(
        total.emission_points.propagation,
        local.emission_points.propagation,
        "propagation emission count"
    );
    require_size_merge_fits(
        total.emission_points.afterburner,
        local.emission_points.afterburner,
        "Afterburner emission count"
    );
    require_u64_merge_fits(
        total.event_status_counts.processed,
        local.event_status_counts.processed,
        "processed event-status count"
    );
    require_u64_merge_fits(
        total.event_status_counts.skipped_due_to_event_empty,
        local.event_status_counts.skipped_due_to_event_empty,
        "empty event-status count"
    );
    require_u64_merge_fits(
        total.event_status_counts.skipped_due_to_event_failure,
        local.event_status_counts.skipped_due_to_event_failure,
        "failed event-status count"
    );
    require_u64_merge_fits(
        total.subevent_status_counts.processed,
        local.subevent_status_counts.processed,
        "processed subevent-status count"
    );
    require_u64_merge_fits(
        total.subevent_status_counts.skipped_due_to_subevent_empty,
        local.subevent_status_counts.skipped_due_to_subevent_empty,
        "empty subevent-status count"
    );
    require_u64_merge_fits(
        total.subevent_status_counts.skipped_due_to_subevent_failure,
        local.subevent_status_counts.skipped_due_to_subevent_failure,
        "failed subevent-status count"
    );
}

/**
 * @brief Merge one event-local preparation summary into run totals.
 * @param total Existing run-total preparation summary.
 * @param local Completed event-local preparation summary to consume.
 * @throws std::logic_error If stable per-species structure differs.
 * @throws std::overflow_error If any run-total counter would overflow.
 *
 * The caller invokes this function only after all workers join and in
 * ascending outer-event order. Rejection records and ordered summaries are
 * therefore reconstructed deterministically and never reflect scheduler
 * completion order.
 */
void accumulate_preparation_summary(
    HBTEventPreparationSummary& total,
    HBTEventPreparationSummary local
) {
    if (total.species.size() != local.species.size()) {
        throw std::logic_error(
            "AnalysisRunner: preparation species-count structure mismatch"
        );
    }

    require_preparation_merge_fits(total, local);

    for (std::size_t index = 0U; index < total.species.size(); ++index) {
        const HBTSpeciesPreparationCounts& total_entry = total.species[index];
        const HBTSpeciesPreparationCounts& local_entry = local.species[index];
        if (total_entry.species != local_entry.species) {
            throw std::logic_error(
                "AnalysisRunner: preparation species-order mismatch"
            );
        }
        require_size_merge_fits(
            total_entry.accepted,
            local_entry.accepted,
            "per-species accepted count"
        );
        require_size_merge_fits(
            total_entry.primordial,
            local_entry.primordial,
            "per-species primordial count"
        );
        require_size_merge_fits(
            total_entry.primordial_rescattering,
            local_entry.primordial_rescattering,
            "per-species PR count"
        );
        require_size_merge_fits(
            total_entry.primordial_rescattering_decay,
            local_entry.primordial_rescattering_decay,
            "per-species PRD count"
        );
    }

    total.outer_events_processed += local.outer_events_processed;
    total.subevents_processed += local.subevents_processed;
    total.raw_particles += local.raw_particles;
    total.unsupported_species += local.unsupported_species;
    total.unrequired_species += local.unrequired_species;
    total.particle_acceptance_rejections +=
        local.particle_acceptance_rejections;
    total.origin_rejections += local.origin_rejections;
    total.accepted_particles += local.accepted_particles;
    total.emission_points.sampler += local.emission_points.sampler;
    total.emission_points.propagation += local.emission_points.propagation;
    total.emission_points.afterburner += local.emission_points.afterburner;
    total.event_status_counts.processed +=
        local.event_status_counts.processed;
    total.event_status_counts.skipped_due_to_event_empty +=
        local.event_status_counts.skipped_due_to_event_empty;
    total.event_status_counts.skipped_due_to_event_failure +=
        local.event_status_counts.skipped_due_to_event_failure;
    total.subevent_status_counts.processed +=
        local.subevent_status_counts.processed;
    total.subevent_status_counts.skipped_due_to_subevent_empty +=
        local.subevent_status_counts.skipped_due_to_subevent_empty;
    total.subevent_status_counts.skipped_due_to_subevent_failure +=
        local.subevent_status_counts.skipped_due_to_subevent_failure;

    for (std::size_t index = 0U; index < total.species.size(); ++index) {
        total.species[index].accepted += local.species[index].accepted;
        total.species[index].primordial += local.species[index].primordial;
        total.species[index].primordial_rescattering +=
            local.species[index].primordial_rescattering;
        total.species[index].primordial_rescattering_decay +=
            local.species[index].primordial_rescattering_decay;
    }

    for (const hbt::RejectedParticleRecord& record :
         local.numerical_rejections.records()) {
        total.numerical_rejections.add(record);
    }
    for (HBTEventExecutionSummary& event : local.events) {
        total.events.push_back(std::move(event));
    }
    for (HBTSubeventPreparationSummary& subevent : local.subevents) {
        total.subevents.push_back(std::move(subevent));
    }
}

/**
 * @brief Resolve the effective Phase-8 outer-event worker count.
 * @param run_config Resolved run configuration.
 * @return Worker count in [1, number_of_events].
 *
 * threads == 0 resolves std::thread::hardware_concurrency() with a fallback
 * of one. Explicit and automatic counts are capped by number_of_events so no
 * idle worker is created for a run with fewer events than requested threads.
 */
std::size_t effective_worker_count(
    const config::RunConfig& run_config
) {
    std::size_t requested = run_config.threads;
    if (requested == 0U) {
        const unsigned int hardware = std::thread::hardware_concurrency();
        requested = hardware == 0U
            ? 1U
            : static_cast<std::size_t>(hardware);
    }

    return std::min(requested, run_config.number_of_events);
}

/**
 * @brief Claim the next dynamically scheduled outer-event index.
 * @param next_event Atomic next-unclaimed zero-based event index.
 * @param event_count Configured number of outer events.
 * @return Claimed zero-based index, or std::nullopt after all claims.
 *
 * The compare/exchange loop avoids incrementing an already exhausted counter,
 * so the scheduler remains well-defined even when event_count is the maximum
 * representable std::size_t value. Relaxed ordering is sufficient because the
 * atomic controls unique index ownership only; worker-result visibility is
 * established by thread join before any merge.
 */
std::optional<std::size_t> claim_next_event(
    std::atomic<std::size_t>& next_event,
    std::size_t event_count
) {
    std::size_t current = next_event.load(std::memory_order_relaxed);
    while (current < event_count) {
        const std::size_t desired = current + 1U;
        if (next_event.compare_exchange_weak(
                current,
                desired,
                std::memory_order_relaxed,
                std::memory_order_relaxed)) {
            return current;
        }
    }
    return std::nullopt;
}

/**
 * @brief Shared progress coordination for one parallel event-processing pool.
 *
 * This object is the only mutable orchestration state intentionally shared by
 * event workers. Its mutex protects the completion count and cancellation
 * wake-up state. Workers never call AnalysisProgressObserver directly; the
 * main orchestration thread is the sole observer caller. No scientific state
 * is stored here.
 */
class ParallelProgressState {
public:
    /**
     * @brief Publish one fully completed subevent.
     * @param cancellation Global fatal-cancellation flag.
     * @return true when completion was recorded, false after cancellation.
     */
    bool publish(const std::atomic<bool>& cancellation) {
        const std::lock_guard<std::mutex> lock(mutex_);
        if (cancelled_ ||
            cancellation.load(std::memory_order_acquire)) {
            return false;
        }
        ++completed_subevents_;
        changed_.notify_one();
        return true;
    }

    /**
     * @brief Report that one worker has exited its worker loop.
     *
     * This method is called exactly once by each created worker.
     */
    void worker_finished() {
        const std::lock_guard<std::mutex> lock(mutex_);
        ++finished_workers_;
        changed_.notify_all();
    }

    /**
     * @brief Wake the observer loop after fatal cancellation begins.
     */
    void cancel() {
        const std::lock_guard<std::mutex> lock(mutex_);
        cancelled_ = true;
        changed_.notify_all();
    }

    /**
     * @brief Wait until progress, cancellation, or worker completion changes.
     * @param reported_subevents Number of completions already sent to observer.
     * @param worker_count Total number of created workers.
     * @return Snapshot used exclusively by the main orchestration thread.
     */
    struct Snapshot {
        std::size_t completed_subevents;  ///< Published completion count.
        std::size_t finished_workers;     ///< Workers that exited their loop.
        bool cancelled;                  ///< Fatal cancellation was signaled.
    };

    [[nodiscard]] Snapshot wait_for_change(
        std::size_t reported_subevents,
        std::size_t worker_count
    ) {
        std::unique_lock<std::mutex> lock(mutex_);
        changed_.wait(lock, [&] {
            return completed_subevents_ > reported_subevents ||
                   finished_workers_ == worker_count ||
                   cancelled_;
        });
        return {
            completed_subevents_,
            finished_workers_,
            cancelled_
        };
    }

private:
    std::mutex mutex_;  ///< Protects all fields below.
    std::condition_variable changed_;  ///< Main-thread wake-up condition.
    std::size_t completed_subevents_{0U};  ///< Completed work units.
    std::size_t finished_workers_{0U};     ///< Exited worker count.
    bool cancelled_{false};                ///< Fatal-cancellation wake-up.
};

/**
 * @brief Process one outer event and all configured independent subevents.
 * @tparam CompletionCallback Callable invoked once after each full subevent.
 * @param outer_event_number One-based outer-event number.
 * @param run_config Resolved global run configuration.
 * @param config Resolved HBT scientific configuration.
 * @param startup Resolved HBT startup requirements.
 * @param summary Preparation diagnostics owned by the caller.
 * @param pair_summary Pair-processing results owned by the caller.
 * @param frame_consumer Raw histogram consumer owned by the caller.
 * @param cancellation Optional fatal-cancellation flag; null in serial use.
 * @param completed Completion callback outside the pair hot path.
 * @return true after full event completion, false after external cancellation.
 * @throws std::overflow_error If pair accounting or status accounting would
 *         overflow.
 * @throws std::logic_error If internal pair-accounting invariants fail.
 *
 * Event-local Sampler/Afterburner input failures are classified as
 * SkippedDueToEventFailure and transactionally discard the complete event. A
 * missing mandatory Sampler match is classified as
 * SkippedDueToSubeventFailure after the remaining declared rows are consumed.
 * Unclassified exceptions remain fatal and propagate to the caller.
 *
 * The SamplerReader, AfterburnerReader, EventBuffers, accepted-particle
 * staging, and per-subevent shuffle RNG created here are thread-confined to
 * the calling worker. Accepted and required Particle values are shuffled once
 * after preparation of each complete subevent and before species grouping and
 * pair formation. Read-only configuration and startup metadata may be shared
 * by all workers. Fatal cancellation is observed only at subevent boundaries,
 * so an already-started subevent is never partially abandoned by orchestration
 * logic.
 */
template <typename CompletionCallback>
bool process_outer_event(
    std::size_t outer_event_number,
    const config::RunConfig& run_config,
    const hbt::HBTConfig& config,
    const hbt::HBTStartupState& startup,
    HBTEventPreparationSummary& summary,
    hbt::HBTPairProcessingSummary& pair_summary,
    hbt::PairFrameConsumer& frame_consumer,
    const std::atomic<bool>* cancellation,
    CompletionCallback&& completed
) {
    if (cancellation != nullptr &&
        cancellation->load(std::memory_order_acquire)) {
        return false;
    }

    std::size_t resolved_subevents = 0U;
    const auto mark_resolved = [&](std::size_t subevent_number) {
        ++resolved_subevents;
        completed(outer_event_number, subevent_number);
    };

    try {
        const std::filesystem::path event_directory =
            outer_event_directory(
                run_config.events_path,
                outer_event_number
            );

        const input::SamplerReader sampler_reader = event_input_operation(
            [&] {
                return input::SamplerReader(sampler_path(event_directory));
            }
        );

        event_input_operation(
            [&] {
                require_sampler_subevent_count(
                    sampler_reader,
                    run_config.number_of_subevents,
                    outer_event_number
                );
            }
        );

        input::AfterburnerReader afterburner_reader = event_input_operation(
            [&] {
                return input::AfterburnerReader(
                    afterburner_path(event_directory)
                );
            }
        );
        bool event_empty = true;

        for (std::size_t index = 0U;
             index < run_config.number_of_subevents;
             ++index) {
            if (cancellation != nullptr &&
                cancellation->load(std::memory_order_acquire)) {
                return false;
            }

            const std::optional<input::AfterburnerSubeventHeader> header =
                event_input_operation(
                    [&] {
                        return afterburner_reader.begin_next_subevent();
                    }
                );

            if (!header.has_value()) {
                throw RecoverableEventFailure(
                    "AnalysisRunner: Afterburner ended before configured "
                    "subevent count for outer event " +
                    std::to_string(outer_event_number)
                );
            }

            const bool subevent_empty = header->particle_count == 0U;
            event_empty = event_empty && subevent_empty;
            hbt::EventBuffers buffers;
            std::size_t accepted_in_subevent = 0U;
            std::optional<std::string> subevent_failure;
            HBTEventPreparationSummary subevent_preparation =
                make_zero_preparation_summary(startup);

            {
                std::vector<hbt::Particle> accepted_particles;

                for (std::size_t particle_index = 0U;
                     particle_index < header->particle_count;
                     ++particle_index) {
                    const input::AfterburnerParticleRecord particle =
                        event_input_operation(
                            [&] {
                                return afterburner_reader.read_particle();
                            }
                        );

                    if (subevent_failure.has_value()) {
                        continue;
                    }

                    try {
                        static_cast<void>(prepare_particle(
                            outer_event_number,
                            header->subevent_id,
                            particle,
                            config,
                            startup,
                            sampler_reader,
                            accepted_particles,
                            subevent_preparation
                        ));
                    } catch (const RecoverableSubeventFailure& error) {
                        subevent_failure = error.what();
                        accepted_particles.clear();
                    }
                }

                event_input_operation(
                    [&] {
                        afterburner_reader.finish_subevent();
                    }
                );

                if (!subevent_failure.has_value()) {
                    accepted_in_subevent = accepted_particles.size();
                    hbt::shuffle_subevent_particles(
                        accepted_particles,
                        outer_event_number,
                        header->subevent_id
                    );

                    for (hbt::Particle& accepted_particle :
                         accepted_particles) {
                        buffers.add(std::move(accepted_particle));
                    }
                }
            }

            if (subevent_failure.has_value()) {
                increment_status_count(
                    summary.subevent_status_counts
                        .skipped_due_to_subevent_failure,
                    "failed subevent-status count"
                );
                summary.subevents.push_back({
                    outer_event_number,
                    header->subevent_id,
                    0U,
                    SubeventStatus::SkippedDueToSubeventFailure,
                    subevent_failure.value()
                });
                ++summary.subevents_processed;
                mark_resolved(index + 1U);
                continue;
            }

            accumulate_preparation_summary(
                summary,
                std::move(subevent_preparation)
            );

            hbt::PairSubeventProcessingResult local_pair_result =
                hbt::process_subevent_pairs(
                    outer_event_number,
                    header->subevent_id,
                    buffers,
                    startup.required_primitive_channels,
                    config.origin_mode,
                    config.pair_slicing,
                    frame_consumer
                );

            hbt::accumulate_pair_processing_result(
                pair_summary,
                std::move(local_pair_result)
            );

            const SubeventStatus subevent_status = subevent_empty
                ? SubeventStatus::SkippedDueToSubeventEmpty
                : SubeventStatus::Processed;
            if (subevent_empty) {
                increment_status_count(
                    summary.subevent_status_counts
                        .skipped_due_to_subevent_empty,
                    "empty subevent-status count"
                );
            } else {
                increment_status_count(
                    summary.subevent_status_counts.processed,
                    "processed subevent-status count"
                );
            }

            summary.subevents.push_back({
                outer_event_number,
                header->subevent_id,
                accepted_in_subevent,
                subevent_status,
                {}
            });
            ++summary.subevents_processed;
            mark_resolved(index + 1U);
        }

        const bool has_extra_subevent = event_input_operation(
            [&] {
                return afterburner_reader.begin_next_subevent().has_value();
            }
        );
        if (has_extra_subevent) {
            throw RecoverableEventFailure(
                "AnalysisRunner: Afterburner contains more subevents than "
                "configured for outer event " +
                std::to_string(outer_event_number)
            );
        }

        const EventStatus event_status = event_empty
            ? EventStatus::SkippedDueToEventEmpty
            : EventStatus::Processed;
        ++summary.outer_events_processed;
        summary.events.push_back({outer_event_number, event_status, {}});
        if (event_empty) {
            increment_status_count(
                summary.event_status_counts.skipped_due_to_event_empty,
                "empty event-status count"
            );
        } else {
            increment_status_count(
                summary.event_status_counts.processed,
                "processed event-status count"
            );
        }
        return true;
    } catch (const RecoverableEventFailure& error) {
        summary = make_zero_preparation_summary(startup);
        pair_summary = make_zero_pair_summary(config, startup);
        ++summary.outer_events_processed;
        increment_status_count(
            summary.event_status_counts.skipped_due_to_event_failure,
            "failed event-status count"
        );
        summary.events.push_back({
            outer_event_number,
            EventStatus::SkippedDueToEventFailure,
            error.what()
        });

        while (resolved_subevents < run_config.number_of_subevents) {
            mark_resolved(resolved_subevents + 1U);
        }
        return true;
    }
}

/**
 * @brief Complete event-processing outputs before post-sample analysis.
 */
struct EventProcessingOutputs {
    HBTEventPreparationSummary preparation;  ///< Ordered preparation summary.
    hbt::HBTPairProcessingSummary pairs;     ///< Ordered pair summary.
    hbt::RawHistogramState histograms;       ///< Complete reduced raw state.
};

/**
 * @brief Process all outer events through the serial reference path.
 * @param run_config Resolved global run configuration.
 * @param config Resolved HBT scientific configuration.
 * @param startup Resolved HBT startup requirements.
 * @param progress Optional non-owning progress observer.
 * @return Complete serial event-processing outputs.
 */
EventProcessingOutputs process_events_serial(
    const config::RunConfig& run_config,
    const hbt::HBTConfig& config,
    const hbt::HBTStartupState& startup,
    AnalysisProgressObserver* progress
) {
    HBTEventPreparationSummary summary =
        make_zero_preparation_summary(startup);
    hbt::HBTPairProcessingSummary pair_summary =
        make_zero_pair_summary(config, startup);
    hbt::RawHistogramState raw_histograms =
        hbt::make_zero_raw_histogram_state(config);
    std::size_t resolved_subevents = 0U;

    for (std::size_t event_index = 0U;
         event_index < run_config.number_of_events;
         ++event_index) {
        HBTEventPreparationSummary local_preparation =
            make_zero_preparation_summary(startup);
        hbt::HBTPairProcessingSummary local_pairs =
            make_zero_pair_summary(config, startup);
        hbt::RawHistogramState local_histograms =
            hbt::make_zero_raw_histogram_state(config);
        hbt::RawHistogramAccumulator local_histogram_accumulator(
            config,
            startup,
            local_histograms
        );

        const bool event_completed = process_outer_event(
            event_index + 1U,
            run_config,
            config,
            startup,
            local_preparation,
            local_pairs,
            local_histogram_accumulator,
            nullptr,
            [&](std::size_t outer_event_number,
                std::size_t subevent_number) {
                ++resolved_subevents;
                if (progress != nullptr) {
                    progress->subevent_completed(
                        resolved_subevents,
                        outer_event_number,
                        subevent_number
                    );
                }
            }
        );
        if (!event_completed) {
            throw std::logic_error(
                "AnalysisRunner: serial event processing was cancelled"
            );
        }
        if (local_preparation.events.size() != 1U) {
            throw std::logic_error(
                "AnalysisRunner: serial event has invalid status cardinality"
            );
        }

        const bool event_failed =
            local_preparation.events.front().status ==
            EventStatus::SkippedDueToEventFailure;

        accumulate_preparation_summary(
            summary,
            std::move(local_preparation)
        );
        hbt::accumulate_pair_processing_summary(
            pair_summary,
            std::move(local_pairs)
        );
        if (!event_failed) {
            hbt::accumulate_raw_histogram_state(
                config,
                raw_histograms,
                local_histograms
            );
        }
    }

    return {
        std::move(summary),
        std::move(pair_summary),
        std::move(raw_histograms)
    };
}

/**
 * @brief Ordered result of exactly one completed outer event.
 *
 * One worker exclusively owns and mutates an EventResult while processing its
 * event. It then moves the result into the pre-sized event-result slot for
 * that event. No other worker accesses that slot. The main thread reads all
 * slots only after every worker has joined.
 */
struct EventResult {
    HBTEventPreparationSummary preparation;  ///< Event-local preparation data.
    hbt::HBTPairProcessingSummary pairs;     ///< Event-local pair data.
};

/**
 * @brief Thread-confined scientific accumulation state of one worker.
 *
 * Exactly one worker mutates this state. The contained histogram layout is
 * identical to the serial complete-run layout, but it receives only events
 * dynamically assigned to that worker. The main thread reads it only after
 * join, so no histogram counter is shared concurrently.
 */
struct WorkerState {
    hbt::RawHistogramState histograms;  ///< Private integer raw counts.
};

/**
 * @brief Join every joinable worker thread.
 * @param workers Worker threads to join.
 */
void join_workers(std::vector<std::thread>& workers) {
    for (std::thread& worker : workers) {
        if (worker.joinable()) {
            worker.join();
        }
    }
}

/**
 * @brief Convert one run-global completion count to canonical progress detail.
 * @param completed_subevents One-based completed-subevent count.
 * @param subevents_per_event Configured subevents per outer event.
 * @return One-based canonical outer-event and subevent positions.
 *
 * Parallel workers may finish out of order. Progress detail deliberately uses
 * the serial event-major logical position so scheduler order is not exposed by
 * the observer interface. The completed_subevents value itself always counts
 * real fully completed subevents.
 */
std::pair<std::size_t, std::size_t> canonical_progress_position(
    std::size_t completed_subevents,
    std::size_t subevents_per_event
) {
    if (completed_subevents == 0U || subevents_per_event == 0U) {
        throw std::logic_error(
            "AnalysisRunner: invalid canonical progress position"
        );
    }

    const std::size_t zero_based = completed_subevents - 1U;
    return {
        zero_based / subevents_per_event + 1U,
        zero_based % subevents_per_event + 1U
    };
}

/**
 * @brief Process outer events with dynamic Phase-8 worker scheduling.
 * @param run_config Resolved global run configuration.
 * @param config Resolved HBT scientific configuration.
 * @param startup Resolved HBT startup requirements.
 * @param worker_count Effective worker count, greater than one.
 * @param progress Optional non-owning progress observer.
 * @return Deterministically reduced event-processing outputs.
 * @throws Any fatal exception produced by event processing or the observer.
 *
 * Read-only configuration/startup objects are shared. Each worker owns one
 * WorkerState and writes only EventResult slots for events it obtained from
 * the atomic increasing scheduler. A fatal exception requests cancellation;
 * no new event is intentionally started after cancellation is observed. All
 * created workers are joined before an exception is propagated. Scientific
 * merge begins only after successful join of the complete pool.
 */
EventProcessingOutputs process_events_parallel(
    const config::RunConfig& run_config,
    const hbt::HBTConfig& config,
    const hbt::HBTStartupState& startup,
    std::size_t worker_count,
    AnalysisProgressObserver* progress
) {
    std::vector<WorkerState> worker_states;
    worker_states.reserve(worker_count);
    for (std::size_t worker_index = 0U;
         worker_index < worker_count;
         ++worker_index) {
        static_cast<void>(worker_index);
        worker_states.push_back({
            hbt::make_zero_raw_histogram_state(config)
        });
    }

    std::vector<std::optional<EventResult>> event_results(
        run_config.number_of_events
    );
    std::atomic<std::size_t> next_event{0U};
    std::atomic<bool> cancellation{false};
    std::mutex exception_mutex;
    std::exception_ptr first_exception;
    ParallelProgressState progress_state;

    std::vector<std::thread> workers;
    workers.reserve(worker_count);

    const auto worker_body = [&](std::size_t worker_index) {
        try {
            while (!cancellation.load(std::memory_order_acquire)) {
                const std::optional<std::size_t> event_index =
                    claim_next_event(
                        next_event,
                        run_config.number_of_events
                    );
                if (!event_index.has_value()) {
                    break;
                }
                if (cancellation.load(std::memory_order_acquire)) {
                    break;
                }

                EventResult local{
                    make_zero_preparation_summary(startup),
                    make_zero_pair_summary(config, startup)
                };
                local.preparation.events.reserve(1U);
                local.preparation.subevents.reserve(
                    run_config.number_of_subevents
                );
                local.pairs.subevents.reserve(
                    run_config.number_of_subevents
                );
                hbt::RawHistogramState event_histograms =
                    hbt::make_zero_raw_histogram_state(config);
                hbt::RawHistogramAccumulator event_histogram_accumulator(
                    config,
                    startup,
                    event_histograms
                );

                const bool event_completed = process_outer_event(
                    event_index.value() + 1U,
                    run_config,
                    config,
                    startup,
                    local.preparation,
                    local.pairs,
                    event_histogram_accumulator,
                    &cancellation,
                    [&](std::size_t, std::size_t) {
                        if (progress != nullptr) {
                            static_cast<void>(progress_state.publish(
                                cancellation
                            ));
                        }
                    }
                );

                if (!event_completed) {
                    break;
                }
                if (local.preparation.events.size() != 1U) {
                    throw std::logic_error(
                        "AnalysisRunner: parallel event has invalid status "
                        "cardinality"
                    );
                }

                const bool event_failed =
                    local.preparation.events.front().status ==
                    EventStatus::SkippedDueToEventFailure;
                if (!event_failed) {
                    hbt::accumulate_raw_histogram_state(
                        config,
                        worker_states[worker_index].histograms,
                        event_histograms
                    );
                }
                event_results[event_index.value()].emplace(std::move(local));
            }
        } catch (...) {
            {
                const std::lock_guard<std::mutex> lock(exception_mutex);
                if (first_exception == nullptr) {
                    first_exception = std::current_exception();
                }
            }
            cancellation.store(true, std::memory_order_release);
            progress_state.cancel();
        }

        progress_state.worker_finished();
    };

    try {
        for (std::size_t worker_index = 0U;
             worker_index < worker_count;
             ++worker_index) {
            workers.emplace_back(worker_body, worker_index);
        }
    } catch (...) {
        cancellation.store(true, std::memory_order_release);
        progress_state.cancel();
        join_workers(workers);
        throw;
    }

    try {
        if (progress != nullptr) {
            std::size_t reported_subevents = 0U;
            while (true) {
                const ParallelProgressState::Snapshot snapshot =
                    progress_state.wait_for_change(
                        reported_subevents,
                        worker_count
                    );

                while (reported_subevents <
                       snapshot.completed_subevents) {
                    ++reported_subevents;
                    const auto position = canonical_progress_position(
                        reported_subevents,
                        run_config.number_of_subevents
                    );
                    progress->subevent_completed(
                        reported_subevents,
                        position.first,
                        position.second
                    );
                }

                if (snapshot.cancelled ||
                    snapshot.finished_workers == worker_count) {
                    break;
                }
            }
        }
    } catch (...) {
        cancellation.store(true, std::memory_order_release);
        progress_state.cancel();
        join_workers(workers);
        throw;
    }

    join_workers(workers);

    if (first_exception != nullptr) {
        std::rethrow_exception(first_exception);
    }

    HBTEventPreparationSummary summary =
        make_zero_preparation_summary(startup);
    hbt::HBTPairProcessingSummary pair_summary =
        make_zero_pair_summary(config, startup);
    summary.events.reserve(run_config.number_of_events);
    summary.subevents.reserve(total_subevent_count(run_config));
    pair_summary.subevents.reserve(total_subevent_count(run_config));

    for (std::size_t event_index = 0U;
         event_index < event_results.size();
         ++event_index) {
        if (!event_results[event_index].has_value()) {
            throw std::logic_error(
                "AnalysisRunner: completed parallel run has "
                "missing event result"
            );
        }
        accumulate_preparation_summary(
            summary,
            std::move(event_results[event_index]->preparation)
        );
        hbt::accumulate_pair_processing_summary(
            pair_summary,
            std::move(event_results[event_index]->pairs)
        );
    }

    hbt::RawHistogramState raw_histograms =
        hbt::make_zero_raw_histogram_state(config);
    for (const WorkerState& worker : worker_states) {
        hbt::accumulate_raw_histogram_state(
            config,
            raw_histograms,
            worker.histograms
        );
    }

    return {
        std::move(summary),
        std::move(pair_summary),
        std::move(raw_histograms)
    };
}
}  // namespace

AnalysisRunner::AnalysisRunner(
    std::filesystem::path run_config_path
) : run_config_path_(std::move(run_config_path)) {
}

AnalysisStartupState AnalysisRunner::prepare_startup() const {
    config::RunConfig run_config =
        config::load_run_config(run_config_path_);

    std::optional<hbt::HBTConfig> hbt_config =
        load_hbt_config_if_enabled(run_config);

    std::optional<hbt::HBTStartupState> hbt_startup_state;

    if (hbt_config.has_value()) {
        hbt_startup_state =
            hbt::build_hbt_startup_state(hbt_config.value());
    }

    return AnalysisStartupState{
        std::move(run_config),
        std::move(hbt_config),
        std::move(hbt_startup_state)
    };
}

AnalysisRunSummary AnalysisRunner::run(
    AnalysisProgressObserver* progress
) const {
    AnalysisStartupState startup = prepare_startup();

    if (!startup.run_config.hbt_enabled) {
        return AnalysisRunSummary{
            std::move(startup),
            std::nullopt,
            std::nullopt,
            std::nullopt,
            std::nullopt
        };
    }

    if (!startup.hbt_config.has_value() ||
        !startup.hbt_startup_state.has_value()) {
        throw std::logic_error(
            "AnalysisRunner: enabled HBT startup state is incomplete"
        );
    }

    const std::size_t total_subevents =
        total_subevent_count(startup.run_config);

    if (progress != nullptr) {
        progress->begin(
            startup.run_config.number_of_events,
            total_subevents
        );
    }

    const std::size_t worker_count =
        effective_worker_count(startup.run_config);
    EventProcessingOutputs processing = worker_count == 1U
        ? process_events_serial(
              startup.run_config,
              startup.hbt_config.value(),
              startup.hbt_startup_state.value(),
              progress
          )
        : process_events_parallel(
              startup.run_config,
              startup.hbt_config.value(),
              startup.hbt_startup_state.value(),
              worker_count,
              progress
          );

    if (progress != nullptr) {
        progress->begin_postprocessing();
    }

    hbt::HistogramAnalysisState fit_results = hbt::analyze_histograms(
        startup.hbt_config.value(),
        processing.histograms
    );

    if (progress != nullptr) {
        progress->analysis_complete();
    }

    return AnalysisRunSummary{
        std::move(startup),
        std::move(processing.preparation),
        std::move(processing.pairs),
        std::move(processing.histograms),
        std::move(fit_results)
    };
}

}  // namespace app
