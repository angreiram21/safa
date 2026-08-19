/**
 * @file raw_histograms.cpp
 * @brief Phase-6 raw HBT histogram allocation and accumulation.
 */

#include "hbt/histograms/raw_histograms.h"

#include "hbt/selection/required_channels.h"

#include <array>
#include <cmath>
#include <limits>
#include <stdexcept>
#include <utility>

namespace hbt {

#ifdef HBT_RAW_HISTOGRAM_TEST_INSTRUMENTATION
namespace test {
namespace {

std::size_t raw_histogram_bin_classification_count_value = 0U;

}  // namespace

void reset_raw_histogram_bin_classification_count() noexcept {
    raw_histogram_bin_classification_count_value = 0U;
}

std::size_t raw_histogram_bin_classification_count() noexcept {
    return raw_histogram_bin_classification_count_value;
}

void record_raw_histogram_bin_classification() noexcept {
    ++raw_histogram_bin_classification_count_value;
}

}  // namespace test
#endif

namespace {

/// Number of OSL marginal histogram slots.
constexpr std::size_t osl_slot_count = 4U;
/// Number of radial histogram slots.
constexpr std::size_t radial_slot_count = 2U;
/// Number of relative-time histogram slots.
constexpr std::size_t delta_t_slot_count = 3U;

/**
 * @brief Classification of one finite observable against configured bins.
 */
enum class BinLocationKind {
    InRange,   ///< Value belongs to one ordinary histogram bin.
    Underflow, ///< Value is below the configured minimum.
    Overflow   ///< Value is at or above the configured maximum.
};

/**
 * @brief Reusable bin classification for one finite observable.
 */
struct BinLocation {
    BinLocationKind kind;      ///< Range classification.
    std::size_t slot;          ///< Logical observable slot.
    std::size_t flat_index;    ///< Valid only when kind is InRange.
};

/**
 * @brief All nine bin classifications calculated once for one pair.
 */
struct PairHistogramLocations {
    std::array<BinLocation, osl_slot_count> osl;        ///< OSL locations.
    std::array<BinLocation, radial_slot_count> radial;  ///< Radial locations.
    /// Relative-time locations.
    std::array<BinLocation, delta_t_slot_count> delta_t;
};

/**
 * @brief Multiply two sizes with explicit overflow detection.
 * @param lhs First factor.
 * @param rhs Second factor.
 * @param context Diagnostic context for an overflow error.
 * @return Exact product when representable.
 * @throws std::overflow_error If the product exceeds std::size_t.
 */
std::size_t checked_size_product(
    std::size_t lhs,
    std::size_t rhs,
    const char* context
) {
    if (lhs != 0U && rhs > std::numeric_limits<std::size_t>::max() / lhs) {
        throw std::overflow_error(context);
    }
    return lhs * rhs;
}

/**
 * @brief Add two sizes with explicit overflow detection.
 * @param lhs First addend.
 * @param rhs Second addend.
 * @param context Diagnostic context for an overflow error.
 * @return Exact sum when representable.
 * @throws std::overflow_error If the sum exceeds std::size_t.
 */
std::size_t checked_size_sum(
    std::size_t lhs,
    std::size_t rhs,
    const char* context
) {
    if (rhs > std::numeric_limits<std::size_t>::max() - lhs) {
        throw std::overflow_error(context);
    }
    return lhs + rhs;
}

/**
 * @brief Increment one raw histogram counter without silent wraparound.
 * @param counter Counter to increment.
 * @throws std::overflow_error If @p counter is already uint64_t maximum.
 */
void increment_histogram_counter(std::uint64_t& counter) {
    if (counter == std::numeric_limits<std::uint64_t>::max()) {
        throw std::overflow_error("raw histogram counter overflow");
    }
    ++counter;
}

/**
 * @brief Calculate the flat raw-count block size for one family.
 * @param slot_count Number of logical observables in the family.
 * @param nbins Number of bins per logical observable.
 * @return Exact flattened counter count.
 * @throws std::overflow_error If the flattened size would overflow.
 */
std::size_t family_counter_count(
    std::size_t slot_count,
    std::size_t nbins
) {
    return checked_size_product(
        slot_count,
        nbins,
        "raw histogram family allocation size overflow"
    );
}

/**
 * @brief Allocate one zeroed nine-histogram set.
 * @param config Validated histogram binning configuration.
 * @return One completely zeroed raw histogram destination.
 */
RawHistogramSet make_zero_histogram_set(
    const HBTHistogramConfig& config
) {
    return {
        {
            std::vector<std::uint64_t>(
                family_counter_count(osl_slot_count, config.osl.nbins),
                0U
            ),
            {},
            {}
        },
        {
            std::vector<std::uint64_t>(
                family_counter_count(
                    radial_slot_count,
                    config.radial.nbins
                ),
                0U
            ),
            {},
            {}
        },
        {
            std::vector<std::uint64_t>(
                family_counter_count(
                    delta_t_slot_count,
                    config.delta_t.nbins
                ),
                0U
            ),
            {},
            {}
        }
    };
}

/**
 * @brief Allocate one zeroed origin destination including all slices.
 * @param config Validated histogram binning configuration.
 * @param slice_count Number of flat kinetic-slice destinations.
 * @return Zeroed global state plus exactly @p slice_count slice states.
 */
OriginRawHistogramState make_zero_origin_state(
    const HBTHistogramConfig& config,
    std::size_t slice_count
) {
    OriginRawHistogramState state{make_zero_histogram_set(config), {}};
    state.slices.reserve(slice_count);
    for (std::size_t index = 0U; index < slice_count; ++index) {
        state.slices.push_back(make_zero_histogram_set(config));
    }
    return state;
}

/**
 * @brief Preflight the complete raw-count allocation arithmetic.
 * @param config Complete validated HBT configuration.
 * @param origin_count Number of stored origin states per product.
 * @param slice_count Number of stored slice states per origin.
 * @throws std::overflow_error If any total-size arithmetic would overflow.
 */
void require_allocation_size_fits(
    const HBTConfig& config,
    std::size_t origin_count,
    std::size_t slice_count
) {
    const std::size_t osl_counts =
        family_counter_count(osl_slot_count, config.histogram_config.osl.nbins);
    const std::size_t radial_counts = family_counter_count(
        radial_slot_count,
        config.histogram_config.radial.nbins
    );
    const std::size_t delta_t_counts = family_counter_count(
        delta_t_slot_count,
        config.histogram_config.delta_t.nbins
    );
    const std::size_t counts_per_set = checked_size_sum(
        checked_size_sum(
            osl_counts,
            radial_counts,
            "raw histogram counters-per-set overflow"
        ),
        delta_t_counts,
        "raw histogram counters-per-set overflow"
    );
    const std::size_t sets_per_origin = checked_size_sum(
        1U,
        slice_count,
        "raw histogram global-plus-slice count overflow"
    );
    const std::size_t product_origin_count = checked_size_product(
        config.selection.products.size(),
        origin_count,
        "raw histogram product-origin count overflow"
    );
    const std::size_t set_count = checked_size_product(
        product_origin_count,
        sets_per_origin,
        "raw histogram destination count overflow"
    );
    const std::size_t total_counts = checked_size_product(
        set_count,
        counts_per_set,
        "raw histogram total counter count overflow"
    );
    static_cast<void>(checked_size_product(
        total_counts,
        sizeof(std::uint64_t),
        "raw histogram total byte size overflow"
    ));
}

/**
 * @brief Classify one finite observable against one validated binning.
 * @param value Finite observable value.
 * @param binning Validated uniform histogram binning.
 * @return Reusable in-range, underflow, or overflow classification.
 * @throws std::logic_error If in-range arithmetic yields an impossible bin.
 */
BinLocation locate_histogram_bin(
    double value,
    const HistogramBinningConfig& binning,
    std::size_t slot
) {
#ifdef HBT_RAW_HISTOGRAM_TEST_INSTRUMENTATION
    test::record_raw_histogram_bin_classification();
#endif
    if (value < binning.minimum) {
        return {BinLocationKind::Underflow, slot, 0U};
    }
    if (value >= binning.maximum) {
        return {BinLocationKind::Overflow, slot, 0U};
    }

    const double scaled =
        (value - binning.minimum) * binning.inverse_bin_width;
    const std::size_t bin = static_cast<std::size_t>(scaled);
    if (bin >= binning.nbins) {
        throw std::logic_error(
            "raw histogram: in-range value produced invalid bin index"
        );
    }
    const std::size_t offset = checked_size_product(
        slot,
        binning.nbins,
        "raw histogram flattened offset overflow"
    );
    const std::size_t flat_index = checked_size_sum(
        offset,
        bin,
        "raw histogram flattened index overflow"
    );
    return {BinLocationKind::InRange, slot, flat_index};
}

/**
 * @brief Calculate all histogram transformations and bins once for one pair.
 * @param observables Finite frame observables supplied by Phase 5.
 * @param config Validated histogram binning configuration.
 * @return Reusable locations for all nine logical histograms.
 */
PairHistogramLocations locate_pair_histograms(
    const PairFrameObservables& observables,
    const HBTHistogramConfig& config
) {
    return {
        {
            locate_histogram_bin(
                std::abs(observables.r_out_lcms_fm),
                config.osl,
                0U
            ),
            locate_histogram_bin(
                std::abs(observables.r_out_prf_fm),
                config.osl,
                1U
            ),
            locate_histogram_bin(
                std::abs(observables.r_side_fm),
                config.osl,
                2U
            ),
            locate_histogram_bin(
                std::abs(observables.r_long_fm),
                config.osl,
                3U
            )
        },
        {
            locate_histogram_bin(
                observables.r_radial_lcms_fm,
                config.radial,
                0U
            ),
            locate_histogram_bin(
                observables.r_radial_prf_fm,
                config.radial,
                1U
            )
        },
        {
            locate_histogram_bin(
                observables.delta_t_lab_fm,
                config.delta_t,
                0U
            ),
            locate_histogram_bin(
                observables.delta_t_lcms_fm,
                config.delta_t,
                1U
            ),
            locate_histogram_bin(
                observables.delta_t_prf_fm,
                config.delta_t,
                2U
            )
        }
    };
}

/**
 * @brief Apply one reusable location to one flattened histogram family.
 * @param bins Flattened slot-major raw counts.
 * @param underflows Underflow diagnostics indexed by logical slot.
 * @param overflows Overflow diagnostics indexed by logical slot.
 * @param location Reusable classification and flat index for the pair.
 * @throws std::logic_error If flattened storage is inconsistent.
 * @throws std::overflow_error If the selected counter would overflow.
 */
template <std::size_t SlotCount>
void apply_location(
    std::vector<std::uint64_t>& bins,
    std::array<std::uint64_t, SlotCount>& underflows,
    std::array<std::uint64_t, SlotCount>& overflows,
    const BinLocation& location
) {
    if (location.slot >= SlotCount) {
        throw std::logic_error("raw histogram: invalid observable slot");
    }
    if (location.kind == BinLocationKind::Underflow) {
        increment_histogram_counter(underflows[location.slot]);
        return;
    }
    if (location.kind == BinLocationKind::Overflow) {
        increment_histogram_counter(overflows[location.slot]);
        return;
    }
    if (location.flat_index >= bins.size()) {
        throw std::logic_error("raw histogram: flattened storage mismatch");
    }
    increment_histogram_counter(bins[location.flat_index]);
}

/**
 * @brief Fill one nine-histogram destination from reusable locations.
 * @param set Mutable histogram destination.
 * @param locations Pair locations calculated once before any fan-out.
 */
void fill_histogram_set(
    RawHistogramSet& set,
    const PairHistogramLocations& locations
) {
    for (std::size_t slot = 0U; slot < osl_slot_count; ++slot) {
        apply_location(
            set.osl.bins,
            set.osl.underflow_counts,
            set.osl.overflow_counts,
            locations.osl[slot]
        );
    }
    for (std::size_t slot = 0U; slot < radial_slot_count; ++slot) {
        apply_location(
            set.radial.bins,
            set.radial.underflow_counts,
            set.radial.overflow_counts,
            locations.radial[slot]
        );
    }
    for (std::size_t slot = 0U; slot < delta_t_slot_count; ++slot) {
        apply_location(
            set.delta_t.bins,
            set.delta_t.underflow_counts,
            set.delta_t.overflow_counts,
            locations.delta_t[slot]
        );
    }
}

/**
 * @brief Fill one origin's global and optional already routed slice states.
 * @param origin Mutable origin state for one final product.
 * @param locations Pair locations calculated once before fan-out.
 * @param slice_route Already resolved slice, or nullptr when slicing is off.
 * @throws std::logic_error If slice presence or index disagrees with layout.
 */
void fill_origin_state(
    OriginRawHistogramState& origin,
    const PairHistogramLocations& locations,
    const PairSliceRoute* slice_route
) {
    fill_histogram_set(origin.global, locations);

    if (origin.slices.empty()) {
        if (slice_route != nullptr) {
            throw std::logic_error(
                "raw histogram: slice route present for unsliced state"
            );
        }
        return;
    }
    if (slice_route == nullptr ||
        slice_route->flat_slice_index >= origin.slices.size()) {
        throw std::logic_error(
            "raw histogram: missing or invalid flat slice index"
        );
    }
    fill_histogram_set(
        origin.slices[slice_route->flat_slice_index],
        locations
    );
}

/**
 * @brief Verify histogram-relevant startup state matches configuration.
 * @param config Complete validated HBT configuration.
 * @param startup Startup state derived from that configuration.
 * @throws std::invalid_argument If product or required-channel state differs.
 */
void require_matching_startup(
    const HBTConfig& config,
    const HBTStartupState& startup
) {
    if (config.selection.products.size() != startup.selection.products.size()) {
        throw std::invalid_argument(
            "raw histogram accumulator: startup product count mismatch"
        );
    }
    for (std::size_t index = 0U;
         index < config.selection.products.size();
         ++index) {
        if (config.selection.products[index].primitive_channels !=
            startup.selection.products[index].primitive_channels) {
            throw std::invalid_argument(
                "raw histogram accumulator: startup selection mismatch"
            );
        }
    }

    const std::vector<PrimitiveChannelId> expected_channels =
        required_primitive_channels(config.selection);
    if (startup.required_primitive_channels != expected_channels) {
        throw std::invalid_argument(
            "raw histogram accumulator: required-channel order mismatch"
        );
    }
}

/**
 * @brief Verify one preallocated destination has exact family block sizes.
 * @param config Validated histogram binning configuration.
 * @param set Preallocated destination to inspect.
 * @throws std::logic_error If any flattened family size differs.
 */
void require_matching_histogram_set_layout(
    const HBTHistogramConfig& config,
    const RawHistogramSet& set
) {
    if (set.osl.bins.size() !=
            family_counter_count(osl_slot_count, config.osl.nbins) ||
        set.radial.bins.size() !=
            family_counter_count(radial_slot_count, config.radial.nbins) ||
        set.delta_t.bins.size() !=
            family_counter_count(delta_t_slot_count, config.delta_t.nbins)) {
        throw std::logic_error(
            "raw histogram: family-state size mismatch"
        );
    }
}

/**
 * @brief Verify preallocated state dimensions before entering the hot path.
 * @param config Validated HBT configuration.
 * @param state Preallocated raw histogram state.
 * @throws std::logic_error If any product/origin/slice dimension differs.
 */
void require_matching_state_layout(
    const HBTConfig& config,
    const RawHistogramState& state
) {
    const std::size_t origin_count =
        raw_histogram_origin_count(config.origin_mode);
    const std::size_t slice_count =
        raw_histogram_slice_count(config.pair_slicing);

    if (state.products.size() != config.selection.products.size()) {
        throw std::logic_error(
            "raw histogram: product-state count mismatch"
        );
    }
    for (const ProductRawHistogramState& product : state.products) {
        if (product.origins.size() != origin_count) {
            throw std::logic_error(
                "raw histogram: origin-state count mismatch"
            );
        }
        for (const OriginRawHistogramState& origin : product.origins) {
            if (origin.slices.size() != slice_count) {
                throw std::logic_error(
                    "raw histogram: slice-state count mismatch"
                );
            }
            require_matching_histogram_set_layout(
                config.histogram_config,
                origin.global
            );
            for (const RawHistogramSet& slice : origin.slices) {
                require_matching_histogram_set_layout(
                    config.histogram_config,
                    slice
                );
            }
        }
    }
}


/**
 * @brief Require one uint64_t addition to be representable.
 * @param total Existing destination counter.
 * @param local Worker-private counter to add.
 * @throws std::overflow_error If total + local would overflow.
 */
void require_histogram_merge_fits(
    std::uint64_t total,
    std::uint64_t local
) {
    if (local > std::numeric_limits<std::uint64_t>::max() - total) {
        throw std::overflow_error("raw histogram merge: counter overflow");
    }
}

/**
 * @brief Check all homologous counters in one histogram family.
 * @tparam SlotCount Number of underflow/overflow slots.
 * @param total_bins Existing flattened destination bins.
 * @param total_underflows Existing destination underflows.
 * @param total_overflows Existing destination overflows.
 * @param local_bins Worker-private flattened bins.
 * @param local_underflows Worker-private underflows.
 * @param local_overflows Worker-private overflows.
 * @throws std::logic_error If flattened family sizes differ.
 * @throws std::overflow_error If any addition would overflow.
 */
template <std::size_t SlotCount>
void require_family_merge_fits(
    const std::vector<std::uint64_t>& total_bins,
    const std::array<std::uint64_t, SlotCount>& total_underflows,
    const std::array<std::uint64_t, SlotCount>& total_overflows,
    const std::vector<std::uint64_t>& local_bins,
    const std::array<std::uint64_t, SlotCount>& local_underflows,
    const std::array<std::uint64_t, SlotCount>& local_overflows
) {
    if (total_bins.size() != local_bins.size()) {
        throw std::logic_error("raw histogram merge: family size mismatch");
    }

    for (std::size_t index = 0U; index < total_bins.size(); ++index) {
        require_histogram_merge_fits(total_bins[index], local_bins[index]);
    }
    for (std::size_t slot = 0U; slot < SlotCount; ++slot) {
        require_histogram_merge_fits(
            total_underflows[slot],
            local_underflows[slot]
        );
        require_histogram_merge_fits(
            total_overflows[slot],
            local_overflows[slot]
        );
    }
}

/**
 * @brief Add homologous counters in one already-validated histogram family.
 * @tparam SlotCount Number of underflow/overflow slots.
 */
template <std::size_t SlotCount>
void merge_family(
    std::vector<std::uint64_t>& total_bins,
    std::array<std::uint64_t, SlotCount>& total_underflows,
    std::array<std::uint64_t, SlotCount>& total_overflows,
    const std::vector<std::uint64_t>& local_bins,
    const std::array<std::uint64_t, SlotCount>& local_underflows,
    const std::array<std::uint64_t, SlotCount>& local_overflows
) {
    for (std::size_t index = 0U; index < total_bins.size(); ++index) {
        total_bins[index] += local_bins[index];
    }
    for (std::size_t slot = 0U; slot < SlotCount; ++slot) {
        total_underflows[slot] += local_underflows[slot];
        total_overflows[slot] += local_overflows[slot];
    }
}

/** @brief Check one complete nine-histogram set before merge. */
void require_set_merge_fits(
    const RawHistogramSet& total,
    const RawHistogramSet& local
) {
    require_family_merge_fits(
        total.osl.bins,
        total.osl.underflow_counts,
        total.osl.overflow_counts,
        local.osl.bins,
        local.osl.underflow_counts,
        local.osl.overflow_counts
    );
    require_family_merge_fits(
        total.radial.bins,
        total.radial.underflow_counts,
        total.radial.overflow_counts,
        local.radial.bins,
        local.radial.underflow_counts,
        local.radial.overflow_counts
    );
    require_family_merge_fits(
        total.delta_t.bins,
        total.delta_t.underflow_counts,
        total.delta_t.overflow_counts,
        local.delta_t.bins,
        local.delta_t.underflow_counts,
        local.delta_t.overflow_counts
    );
}

/** @brief Merge one complete nine-histogram set after validation. */
void merge_set(
    RawHistogramSet& total,
    const RawHistogramSet& local
) {
    merge_family(
        total.osl.bins,
        total.osl.underflow_counts,
        total.osl.overflow_counts,
        local.osl.bins,
        local.osl.underflow_counts,
        local.osl.overflow_counts
    );
    merge_family(
        total.radial.bins,
        total.radial.underflow_counts,
        total.radial.overflow_counts,
        local.radial.bins,
        local.radial.underflow_counts,
        local.radial.overflow_counts
    );
    merge_family(
        total.delta_t.bins,
        total.delta_t.underflow_counts,
        total.delta_t.overflow_counts,
        local.delta_t.bins,
        local.delta_t.underflow_counts,
        local.delta_t.overflow_counts
    );
}

}  // namespace

std::size_t raw_histogram_origin_count(OriginMode mode) {
    switch (mode) {
        case OriginMode::Primordial:
        case OriginMode::PrimordialRescattering:
        case OriginMode::PrimordialRescatteringDecay:
            return 1U;
        case OriginMode::All:
            return 3U;
    }
    throw std::invalid_argument("raw histogram: invalid origin mode");
}

HistogramOrigin raw_histogram_origin_at(
    OriginMode mode,
    std::size_t index
) {
    switch (mode) {
        case OriginMode::Primordial:
            if (index == 0U) {
                return HistogramOrigin::Primordial;
            }
            break;
        case OriginMode::PrimordialRescattering:
            if (index == 0U) {
                return HistogramOrigin::PrimordialRescattering;
            }
            break;
        case OriginMode::PrimordialRescatteringDecay:
            if (index == 0U) {
                return HistogramOrigin::PrimordialRescatteringDecay;
            }
            break;
        case OriginMode::All:
            if (index == 0U) {
                return HistogramOrigin::Primordial;
            }
            if (index == 1U) {
                return HistogramOrigin::PrimordialRescattering;
            }
            if (index == 2U) {
                return HistogramOrigin::PrimordialRescatteringDecay;
            }
            break;
    }
    throw std::invalid_argument("raw histogram: invalid origin-state index");
}

std::size_t raw_histogram_slice_count(
    const PairSlicingConfig& slicing
) {
    const auto axis_count = [](const PairSlicingAxisConfig& axis) {
        if (!axis.enabled) {
            return std::size_t{0U};
        }
        if (axis.bin_edges_gev.size() < 2U) {
            throw std::invalid_argument(
                "raw histogram: enabled slicing axis has fewer than two edges"
            );
        }
        return axis.bin_edges_gev.size() - 1U;
    };

    const std::size_t kt_count = axis_count(slicing.kt);
    const std::size_t mt_count = axis_count(slicing.mt);
    if (kt_count == 0U) {
        return mt_count;
    }
    if (mt_count == 0U) {
        return kt_count;
    }
    return checked_size_product(
        kt_count,
        mt_count,
        "raw histogram Cartesian slice count overflow"
    );
}

void require_raw_histogram_state_layout(
    const HBTConfig& config,
    const RawHistogramState& state
) {
    require_matching_state_layout(config, state);
}

void accumulate_raw_histogram_state(
    const HBTConfig& config,
    RawHistogramState& total,
    const RawHistogramState& worker
) {
    require_raw_histogram_state_layout(config, total);
    require_raw_histogram_state_layout(config, worker);

    for (std::size_t product_index = 0U;
         product_index < total.products.size();
         ++product_index) {
        const ProductRawHistogramState& local_product =
            worker.products[product_index];
        const ProductRawHistogramState& total_product =
            total.products[product_index];

        for (std::size_t origin_index = 0U;
             origin_index < total_product.origins.size();
             ++origin_index) {
            const OriginRawHistogramState& total_origin =
                total_product.origins[origin_index];
            const OriginRawHistogramState& local_origin =
                local_product.origins[origin_index];

            require_set_merge_fits(total_origin.global, local_origin.global);
            for (std::size_t slice_index = 0U;
                 slice_index < total_origin.slices.size();
                 ++slice_index) {
                require_set_merge_fits(
                    total_origin.slices[slice_index],
                    local_origin.slices[slice_index]
                );
            }
        }
    }

    for (std::size_t product_index = 0U;
         product_index < total.products.size();
         ++product_index) {
        const ProductRawHistogramState& local_product =
            worker.products[product_index];
        ProductRawHistogramState& total_product = total.products[product_index];

        for (std::size_t origin_index = 0U;
             origin_index < total_product.origins.size();
             ++origin_index) {
            OriginRawHistogramState& total_origin =
                total_product.origins[origin_index];
            const OriginRawHistogramState& local_origin =
                local_product.origins[origin_index];

            merge_set(total_origin.global, local_origin.global);
            for (std::size_t slice_index = 0U;
                 slice_index < total_origin.slices.size();
                 ++slice_index) {
                merge_set(
                    total_origin.slices[slice_index],
                    local_origin.slices[slice_index]
                );
            }
        }
    }
}

RawHistogramState make_zero_raw_histogram_state(
    const HBTConfig& config
) {
    const std::size_t origin_count =
        raw_histogram_origin_count(config.origin_mode);
    const std::size_t slice_count =
        raw_histogram_slice_count(config.pair_slicing);
    require_allocation_size_fits(config, origin_count, slice_count);

    RawHistogramState state;
    state.products.reserve(config.selection.products.size());

    for (std::size_t product_index = 0U;
         product_index < config.selection.products.size();
         ++product_index) {
        ProductRawHistogramState product;
        product.origins.reserve(origin_count);
        for (std::size_t origin_index = 0U;
             origin_index < origin_count;
             ++origin_index) {
            static_cast<void>(origin_index);
            product.origins.push_back(
                make_zero_origin_state(config.histogram_config, slice_count)
            );
        }
        state.products.push_back(std::move(product));
    }
    return state;
}

RawHistogramAccumulator::RawHistogramAccumulator(
    const HBTConfig& config,
    const HBTStartupState& startup,
    RawHistogramState& state
) : config_(config),
    state_(state),
    fanout_(build_product_fanout_plan(
        startup.selection,
        startup.required_primitive_channels
    )),
    slice_count_(raw_histogram_slice_count(config.pair_slicing)) {
    require_matching_startup(config, startup);
    require_raw_histogram_state_layout(config, state);
    if (fanout_.offsets.size() !=
        startup.required_primitive_channels.size() + 1U) {
        throw std::logic_error(
            "raw histogram accumulator: fan-out channel count mismatch"
        );
    }
}

void RawHistogramAccumulator::consume(
    const PairFrameRouteContext& context,
    const PairKinematics& kinematics,
    const PairFrameObservables& observables
) {
    static_cast<void>(kinematics);

    if (fanout_.offsets.empty() ||
        context.channel_index >= fanout_.offsets.size() - 1U) {
        throw std::logic_error(
            "raw histogram accumulator: invalid channel_index"
        );
    }
    const bool slicing_active = slice_count_ != 0U;
    if (slicing_active != (context.pair_slice_route != nullptr)) {
        throw std::logic_error(
            "raw histogram accumulator: slice-routing state mismatch"
        );
    }
    if (context.pair_slice_route != nullptr &&
        context.pair_slice_route->flat_slice_index >= slice_count_) {
        throw std::logic_error(
            "raw histogram accumulator: invalid flat_slice_index"
        );
    }

    const PairHistogramLocations locations = locate_pair_histograms(
        observables,
        config_.histogram_config
    );
    const std::size_t begin = fanout_.offsets[context.channel_index];
    const std::size_t end = fanout_.offsets[context.channel_index + 1U];

    const auto fill_products = [&](std::size_t origin_index) {
        for (std::size_t fanout_index = begin;
             fanout_index < end;
             ++fanout_index) {
            const std::size_t product_index =
                fanout_.product_indices[fanout_index];
            if (product_index >= state_.products.size() ||
                origin_index >=
                    state_.products[product_index].origins.size()) {
                throw std::logic_error(
                    "raw histogram accumulator: fan-out state mismatch"
                );
            }
            fill_origin_state(
                state_.products[product_index].origins[origin_index],
                locations,
                context.pair_slice_route
            );
        }
    };

    switch (config_.origin_mode) {
        case OriginMode::Primordial:
            if (!context.origin_routes.primordial ||
                context.origin_routes.primordial_rescattering ||
                context.origin_routes.primordial_rescattering_decay) {
                throw std::logic_error(
                    "raw histogram accumulator: invalid P routing"
                );
            }
            fill_products(0U);
            return;
        case OriginMode::PrimordialRescattering:
            if (context.origin_routes.primordial ||
                !context.origin_routes.primordial_rescattering ||
                context.origin_routes.primordial_rescattering_decay) {
                throw std::logic_error(
                    "raw histogram accumulator: invalid PR routing"
                );
            }
            fill_products(0U);
            return;
        case OriginMode::PrimordialRescatteringDecay:
            if (context.origin_routes.primordial ||
                context.origin_routes.primordial_rescattering ||
                !context.origin_routes.primordial_rescattering_decay) {
                throw std::logic_error(
                    "raw histogram accumulator: invalid PRD routing"
                );
            }
            fill_products(0U);
            return;
        case OriginMode::All:
            if (!context.origin_routes.primordial_rescattering_decay ||
                (context.origin_routes.primordial &&
                 !context.origin_routes.primordial_rescattering)) {
                throw std::logic_error(
                    "raw histogram accumulator: invalid nested all routing"
                );
            }
            if (context.origin_routes.primordial) {
                fill_products(0U);
            }
            if (context.origin_routes.primordial_rescattering) {
                fill_products(1U);
            }
            fill_products(2U);
            return;
    }
    throw std::logic_error(
        "raw histogram accumulator: invalid configured origin mode"
    );
}

}  // namespace hbt
