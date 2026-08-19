/**
 * @file raw_histograms.h
 * @brief Phase-6 raw HBT histogram state and accumulation consumer.
 */

#ifndef HBT_HISTOGRAMS_RAW_HISTOGRAMS_H
#define HBT_HISTOGRAMS_RAW_HISTOGRAMS_H

#include "hbt/config/hbt_config.h"
#include "hbt/histograms/product_fanout_plan.h"
#include "hbt/pair/pair_frame_consumer.h"
#include "hbt/startup/hbt_startup_state.h"

#include <array>
#include <cstddef>
#include <cstdint>
#include <vector>

namespace hbt {

/**
 * @brief Stable OSL histogram slots sharing one configured binning.
 */
enum class OSLHistogramSlot : std::size_t {
    ROutLcms,  ///< Absolute LCMS out component.
    ROutPrf,   ///< Absolute PRF out component.
    RSide,     ///< Absolute side component, invariant LCMS to PRF.
    RLong,     ///< Absolute long component, invariant LCMS to PRF.
    Count      ///< Number of OSL histogram slots.
};

/**
 * @brief Stable radial histogram slots sharing one configured binning.
 */
enum class RadialHistogramSlot : std::size_t {
    RadialLcms,  ///< LCMS radial source separation.
    RadialPrf,   ///< PRF radial source separation.
    Count        ///< Number of radial histogram slots.
};

/**
 * @brief Stable relative-time histogram slots sharing one configured binning.
 */
enum class DeltaTHistogramSlot : std::size_t {
    Lab,    ///< Lab-frame relative time.
    Lcms,   ///< LCMS relative time.
    Prf,    ///< PRF relative time.
    Count   ///< Number of relative-time histogram slots.
};

/**
 * @brief Physical origin identity of one stored histogram destination.
 */
enum class HistogramOrigin {
    Primordial,                     ///< P destination.
    PrimordialRescattering,         ///< PR destination.
    PrimordialRescatteringDecay     ///< PRD destination.
};

/**
 * @brief Mutable raw counts for the four OSL marginal histograms.
 *
 * bins is one contiguous block laid out as slot-major storage. Each slot owns
 * exactly histogram_config.osl.nbins consecutive counters. Underflow and
 * overflow diagnostics are kept separately for every marginal observable.
 */
struct OSLRawHistogramCounts {
    std::vector<std::uint64_t> bins;  ///< Four contiguous OSL count arrays.
    /// Underflow count for each OSL observable slot.
    std::array<std::uint64_t, 4U> underflow_counts{};
    /// Overflow count for each OSL observable slot.
    std::array<std::uint64_t, 4U> overflow_counts{};
};

/**
 * @brief Mutable raw counts for the two radial histograms.
 */
struct RadialRawHistogramCounts {
    std::vector<std::uint64_t> bins;  ///< Two contiguous radial count arrays.
    /// Underflow count for each radial observable slot.
    std::array<std::uint64_t, 2U> underflow_counts{};
    /// Overflow count for each radial observable slot.
    std::array<std::uint64_t, 2U> overflow_counts{};
};

/**
 * @brief Mutable raw counts for the three relative-time histograms.
 */
struct DeltaTRawHistogramCounts {
    std::vector<std::uint64_t> bins;  ///< Three contiguous delta-t arrays.
    /// Underflow count for each relative-time observable slot.
    std::array<std::uint64_t, 3U> underflow_counts{};
    /// Overflow count for each relative-time observable slot.
    std::array<std::uint64_t, 3U> overflow_counts{};
};

/**
 * @brief Complete nine-histogram raw state for one accumulation destination.
 */
struct RawHistogramSet {
    OSLRawHistogramCounts osl;        ///< Four OSL marginal histograms.
    RadialRawHistogramCounts radial;  ///< Two radial histograms.
    DeltaTRawHistogramCounts delta_t; ///< Three relative-time histograms.
};

/**
 * @brief Global and optional kinetic-slice histograms for one origin.
 */
struct OriginRawHistogramState {
    /// Global destination; with slicing, the union of configured slices.
    RawHistogramSet global;
    std::vector<RawHistogramSet> slices; ///< Direct flat-slice destinations.
};

/**
 * @brief Raw histogram states for all requested origins of one final product.
 *
 * origins contains exactly one entry for an individual OriginMode and exactly
 * three entries in P, PR, PRD order for OriginMode::All.
 */
struct ProductRawHistogramState {
    std::vector<OriginRawHistogramState> origins;  ///< Active origin states.
};

/**
 * @brief Raw histogram state accumulated over the complete processed sample.
 *
 * products is aligned exactly with HBTSelection::products. No state exists for
 * primitive channels unless that primitive channel is itself requested as a
 * final product.
 */
struct RawHistogramState {
    std::vector<ProductRawHistogramState> products;  ///< Final-product state.
};

/**
 * @brief Return the number of stored origin states for one configured mode.
 * @param mode Validated requested origin mode.
 * @return One for an individual mode and three for OriginMode::All.
 * @throws std::invalid_argument If @p mode is invalid.
 */
[[nodiscard]] std::size_t raw_histogram_origin_count(OriginMode mode);

/**
 * @brief Return the physical origin stored at one origin-state index.
 * @param mode Validated requested origin mode.
 * @param index Zero-based stored origin index.
 * @return Physical P, PR, or PRD identity at @p index.
 * @throws std::invalid_argument If @p mode or @p index is invalid.
 */
[[nodiscard]] HistogramOrigin raw_histogram_origin_at(
    OriginMode mode,
    std::size_t index
);

/**
 * @brief Return the number of allocated kinetic-slice destinations.
 * @param slicing Validated optional kT/mT slicing configuration.
 * @return Zero when slicing is disabled, otherwise configured slice count.
 * @throws std::invalid_argument If an enabled axis has fewer than two edges.
 * @throws std::overflow_error If the Cartesian slice count would overflow.
 */
[[nodiscard]] std::size_t raw_histogram_slice_count(
    const PairSlicingConfig& slicing
);

/**
 * @brief Allocate and zero all Phase-6 raw histogram state before events.
 * @param config Fully validated HBT scientific configuration.
 * @return Exact product/origin/global/slice raw histogram state.
 * @throws std::invalid_argument If origin or slicing structure is invalid.
 * @throws std::overflow_error If any required allocation size would overflow.
 *
 * All final sizes are resolved before allocation. No lazy histogram creation
 * is performed by the returned state during event processing.
 */
[[nodiscard]] RawHistogramState make_zero_raw_histogram_state(
    const HBTConfig& config
);

/**
 * @brief Require raw histogram state to match its configured eager layout.
 * @param config Fully validated HBT scientific configuration.
 * @param state Complete-sample raw histogram state to inspect.
 * @throws std::invalid_argument If origin or slicing structure is invalid.
 * @throws std::logic_error If any product, origin, slice, or family dimension
 *         differs from the exact configured layout.
 * @throws std::overflow_error If configured size arithmetic is not
 *         representable.
 *
 * This check performs no allocation and does not inspect or modify counter
 * values.
 */
void require_raw_histogram_state_layout(
    const HBTConfig& config,
    const RawHistogramState& state
);

/**
 * @brief Add one worker-private raw histogram state into a run-total state.
 * @param config Validated HBT configuration defining the common layout.
 * @param total Existing run-total raw histogram state to update.
 * @param worker Worker-private raw histogram state to reduce.
 * @throws std::logic_error If either state differs from the configured layout.
 * @throws std::overflow_error If any uint64_t counter addition would overflow.
 *
 * The complete merge is checked before mutation. A failed merge therefore
 * leaves @p total unchanged. Corresponding product, origin, slice, observable,
 * and kinetic-slice bins are added only to their exact homologous counters.
 * No floating-point reduction is performed. The function is synchronous and
 * not internally synchronized: callers must provide exclusive mutable access
 * to @p total and ensure @p worker is no longer being mutated. Phase 8 calls it
 * only from the orchestration thread after worker join.
 */
void accumulate_raw_histogram_state(
    const HBTConfig& config,
    RawHistogramState& total,
    const RawHistogramState& worker
);

/**
 * @brief Phase-6 synchronous accumulator for already calculated pair results.
 *
 * The accumulator owns the immutable startup fan-out plan and borrows the
 * validated configuration plus complete-run raw histogram state. It performs
 * no pair formation, frame transformation, kinematic calculation, slicing,
 * output, normalization, fit, or statistical analysis.
 *
 * The borrowed configuration and state must outlive this accumulator. consume()
 * is synchronous and never retains pointers or references from its arguments.
 * The accumulator is not thread-safe for shared mutation and is intended to be
 * thread-confined; Phase 8 constructs one per worker against that worker's
 * private RawHistogramState.
 */
class RawHistogramAccumulator final : public PairFrameConsumer {
public:
    /**
     * @brief Construct one complete-run raw histogram accumulator.
     * @param config Validated HBT configuration that defines histogram layout.
     * @param startup Validated startup selection and required-channel order.
     * @param state Preallocated zeroed state that receives all raw counts.
     * @throws std::invalid_argument If startup/config dimensions disagree.
     * @throws std::logic_error If preallocated state dimensions disagree.
     */
    RawHistogramAccumulator(
        const HBTConfig& config,
        const HBTStartupState& startup,
        RawHistogramState& state
    );

    /**
     * @brief Accumulate one pair into every required final destination.
     * @param context Reused channel, origin, and optional slice routing.
     * @param kinematics Existing pair kinematics, not recalculated here.
     * @param observables Existing finite frame observables.
     * @throws std::logic_error If routing or state is inconsistent.
     * @throws std::overflow_error If any raw or diagnostic counter overflows.
     *
     * Absolute OSL values and all nine bin classifications are calculated once
     * per consumed pair and reused across products, origins, global state, and
     * the already resolved flat slice.
     */
    void consume(
        const PairFrameRouteContext& context,
        const PairKinematics& kinematics,
        const PairFrameObservables& observables
    ) override;

private:
    const HBTConfig& config_;       ///< Borrowed validated configuration.
    RawHistogramState& state_;      ///< Borrowed complete-run mutable state.
    ProductFanoutPlan fanout_;      ///< Immutable startup channel fan-out.
    std::size_t slice_count_;        ///< Startup-resolved slice cardinality.
};

}  // namespace hbt

#endif  // HBT_HISTOGRAMS_RAW_HISTOGRAMS_H
