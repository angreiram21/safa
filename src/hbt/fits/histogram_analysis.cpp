/**
 * @file histogram_analysis.cpp
 * @brief Complete post-sample HBT statistical orchestration.
 */

#include "hbt/fits/histogram_analysis.h"

#include "hbt/fits/minuit2_fitter.h"
#include "hbt/fits/statistical_analysis.h"

#include <array>
#include <cstddef>
#include <cstdint>
#include <optional>
#include <stdexcept>
#include <utility>
#include <vector>

namespace hbt {
namespace {

/**
 * @brief Return stable unattempted MIGRAD diagnostics for empty results.
 * @return Diagnostic with attempted == false and no numerical minimum.
 */
MigradDiagnostic unattempted_migrad() {
    return {false, false, false, false, false, false, 0, std::nullopt};
}

/**
 * @brief Return stable unattempted MINOS diagnostics for empty results.
 * @return Diagnostic with attempted == false and all side flags false.
 */
MinosDiagnostic unattempted_minos() {
    return {false, false, false, false, false, false, false, false, false};
}

/**
 * @brief Build an explicit invalid Gaussian fit result without fabricated data.
 * @param reason Stable reason explaining why no Gaussian result is published.
 * @return Invalid result carrying only diagnostics and @p reason.
 */
GaussianFitResult invalid_gaussian_result(FitFailureReason reason) {
    return {
        false,
        reason,
        unattempted_migrad(),
        unattempted_minos(),
        std::nullopt,
        std::nullopt,
        {}
    };
}

/**
 * @brief Build an explicit invalid mixed fit result without fabricated data.
 * @param reason Stable reason explaining why no mixed result is published.
 * @param estimator Objective identity retained even when no fit is attempted.
 * @return Invalid result carrying only diagnostics, @p reason, and @p estimator.
 */
MixedFitResult invalid_mixed_result(
    FitFailureReason reason,
    FitEstimator estimator
) {
    const MigradDiagnostic empty = unattempted_migrad();
    const std::array<MigradDiagnostic, MixedFitResult::kCoreStartCount> starts{
        empty, empty, empty, empty, empty
    };
    return {
        false,
        reason,
        estimator,
        starts,
        0U,
        0U,
        0U,
        std::nullopt,
        empty,
        unattempted_minos(),
        unattempted_minos(),
        unattempted_minos(),
        std::nullopt,
        std::nullopt,
        std::nullopt,
        std::nullopt,
        {}
    };
}

/**
 * @brief Build a derived shape placeholder for an intentionally skipped fit.
 * @return Empty, explicitly not-applicable Gaussian and mixed result.
 *
 * The approved production contract treats OSL as global-only. Raw kinetic
 * slice histograms may still exist because the accumulation layout is shared,
 * but no OSL slice likelihood is evaluated.
 */
ShapeHistogramResult not_applicable_shape_result() {
    return {
        std::nullopt,
        std::nullopt,
        0U,
        {},
        invalid_gaussian_result(FitFailureReason::NotApplicable),
        invalid_mixed_result(
            FitFailureReason::NotApplicable, FitEstimator::Poisson
        ),
        invalid_mixed_result(
            FitFailureReason::NotApplicable, FitEstimator::Neyman
        ),
        invalid_mixed_result(
            FitFailureReason::NotApplicable, FitEstimator::Pearson
        )
    };
}

/** Production quality cut validated for radial mT slices in the contract. */
constexpr std::uint64_t kMinimumRadialSliceSelectedCount = 10000U;

/**
 * @brief Analyze one OSL or radial logical histogram without normalization.
 * @param family OSL or radial physical model family.
 * @param bins Slot-major raw uint64_t histogram storage.
 * @param offset First raw counter for this logical histogram.
 * @param binning Validated uniform binning.
 * @param apply_radial_slice_quality_cut Whether the provisional radial kinetic
 *        slice threshold N_selected >= 10000 is required for this histogram.
 * @return Full/core regions plus one Gaussian result and three independent
 *         mixed results (Poisson, Neyman, Pearson).
 */
ShapeHistogramResult analyze_shape_histogram(
    FitObservableFamily family,
    const std::vector<std::uint64_t>& bins,
    std::size_t offset,
    const HistogramBinningConfig& binning,
    bool apply_radial_slice_quality_cut
) {
    const std::optional<StatisticalRegion> region =
        select_shape_region(bins, offset, binning);
    if (!region.has_value()) {
        return {
            std::nullopt,
            std::nullopt,
            0U,
            {},
            invalid_gaussian_result(FitFailureReason::EmptyHistogram),
            invalid_mixed_result(
                FitFailureReason::EmptyHistogram, FitEstimator::Poisson
            ),
            invalid_mixed_result(
                FitFailureReason::EmptyHistogram, FitEstimator::Neyman
            ),
            invalid_mixed_result(
                FitFailureReason::EmptyHistogram, FitEstimator::Pearson
            )
        };
    }

    const std::optional<StatisticalRegion> core_region =
        select_gaussian_core_region(
            family,
            bins,
            offset,
            binning,
            region.value()
        );
    if (apply_radial_slice_quality_cut &&
        region->selected_count < kMinimumRadialSliceSelectedCount) {
        return {
            region,
            core_region,
            region->selected_count,
            {},
            invalid_gaussian_result(FitFailureReason::InsufficientStatistics),
            invalid_mixed_result(
                FitFailureReason::InsufficientStatistics, FitEstimator::Poisson
            ),
            invalid_mixed_result(
                FitFailureReason::InsufficientStatistics, FitEstimator::Neyman
            ),
            invalid_mixed_result(
                FitFailureReason::InsufficientStatistics, FitEstimator::Pearson
            )
        };
    }
    if (!core_region.has_value()) {
        return {
            region,
            std::nullopt,
            region->selected_count,
            {},
            invalid_gaussian_result(FitFailureReason::InsufficientBins),
            invalid_mixed_result(
                FitFailureReason::InvalidGaussianCoreAnchor, FitEstimator::Poisson
            ),
            invalid_mixed_result(
                FitFailureReason::InvalidGaussianCoreAnchor, FitEstimator::Neyman
            ),
            invalid_mixed_result(
                FitFailureReason::InvalidGaussianCoreAnchor, FitEstimator::Pearson
            )
        };
    }

    GaussianFitResult gaussian = fit_gaussian_model(
        family,
        bins,
        offset,
        binning,
        core_region.value()
    );
    MixedFitResult mixed = fit_mixed_model(
        family,
        bins,
        offset,
        binning,
        region.value(),
        FitEstimator::Poisson,
        gaussian
    );
    MixedFitResult mixed_neyman = fit_mixed_model(
        family,
        bins,
        offset,
        binning,
        region.value(),
        FitEstimator::Neyman,
        gaussian
    );
    MixedFitResult mixed_pearson = fit_mixed_model(
        family,
        bins,
        offset,
        binning,
        region.value(),
        FitEstimator::Pearson,
        gaussian
    );
    return {
        region,
        core_region,
        region->selected_count,
        {},
        std::move(gaussian),
        std::move(mixed),
        std::move(mixed_neyman),
        std::move(mixed_pearson)
    };
}

/**
 * @brief Analyze one signed delta-t logical histogram without normalization.
 * @param bins Slot-major raw uint64_t histogram storage.
 * @param offset First raw counter for this logical histogram.
 * @param binning Validated signed delta-t uniform binning.
 * @return Selected region and required raw-count moments.
 */
DeltaTHistogramResult analyze_delta_t_histogram(
    const std::vector<std::uint64_t>& bins,
    std::size_t offset,
    const HistogramBinningConfig& binning
) {
    const std::optional<StatisticalRegion> region =
        select_delta_t_region(bins, offset, binning);
    if (!region.has_value()) {
        return {
            std::nullopt,
            0U,
            {},
            DeltaTStatisticsStatus::EmptyHistogram,
            std::nullopt,
            std::nullopt,
            std::nullopt
        };
    }
    return calculate_delta_t_statistics(
        bins,
        offset,
        binning,
        region.value()
    );
}

/**
 * @brief Analyze all nine logical histograms of one raw destination.
 * @param config Validated histogram binning metadata.
 * @param raw Immutable nine-histogram raw state.
 * @param kinetic_slice true only for a kT/mT slice destination. OSL fitting is
 *        skipped for such destinations because OSL is global-only.
 * @param apply_radial_mt_quality_cut Whether the provisional N_selected >=
 *        10000 quality cut applies to this one-dimensional radial mT slice.
 * @return Derived state with fits/moments complete and normalization empty.
 */
DerivedHistogramSet analyze_histogram_set(
    const HBTHistogramConfig& config,
    const RawHistogramSet& raw,
    bool kinetic_slice,
    bool apply_radial_mt_quality_cut
) {
    DerivedHistogramSet result{};
    for (std::size_t slot = 0U; slot < result.osl.size(); ++slot) {
        result.osl[slot] = kinetic_slice
            ? not_applicable_shape_result()
            : analyze_shape_histogram(
                FitObservableFamily::OSL,
                raw.osl.bins,
                slot * config.osl.nbins,
                config.osl,
                false
            );
    }
    for (std::size_t slot = 0U; slot < result.radial.size(); ++slot) {
        result.radial[slot] = analyze_shape_histogram(
            FitObservableFamily::Radial,
            raw.radial.bins,
            slot * config.radial.nbins,
            config.radial,
            apply_radial_mt_quality_cut
        );
    }
    for (std::size_t slot = 0U; slot < result.delta_t.size(); ++slot) {
        result.delta_t[slot] = analyze_delta_t_histogram(
            raw.delta_t.bins,
            slot * config.delta_t.nbins,
            config.delta_t
        );
    }
    return result;
}

/**
 * @brief Materialize final normalized distributions after analysis is complete.
 * @param config Validated histogram binning metadata.
 * @param raw Immutable nine-histogram raw state.
 * @param derived Matching derived result whose regions are already final.
 */
void normalize_histogram_set(
    const HBTHistogramConfig& config,
    const RawHistogramSet& raw,
    DerivedHistogramSet& derived
) {
    for (std::size_t slot = 0U; slot < derived.osl.size(); ++slot) {
        if (derived.osl[slot].region.has_value()) {
            derived.osl[slot].normalized_bins = normalize_region(
                raw.osl.bins,
                slot * config.osl.nbins,
                config.osl,
                derived.osl[slot].region.value()
            );
        }
    }
    for (std::size_t slot = 0U; slot < derived.radial.size(); ++slot) {
        if (derived.radial[slot].region.has_value()) {
            derived.radial[slot].normalized_bins = normalize_region(
                raw.radial.bins,
                slot * config.radial.nbins,
                config.radial,
                derived.radial[slot].region.value()
            );
        }
    }
    for (std::size_t slot = 0U; slot < derived.delta_t.size(); ++slot) {
        if (derived.delta_t[slot].region.has_value()) {
            derived.delta_t[slot].normalized_bins = normalize_region(
                raw.delta_t.bins,
                slot * config.delta_t.nbins,
                config.delta_t,
                derived.delta_t[slot].region.value()
            );
        }
    }
}

/**
 * @brief Return expected normalized-bin count for one optional region.
 * @param region Selected region or std::nullopt for an empty histogram.
 * @return Zero or inclusive selected region length.
 */
std::size_t region_bin_count(
    const std::optional<StatisticalRegion>& region
) {
    if (!region.has_value()) {
        return 0U;
    }
    return region->last_bin - region->first_bin + 1U;
}

/**
 * @brief Validate one shape result's internal cardinality invariants.
 * @param result Derived shape result to inspect.
 * @throws std::logic_error If region/count/series cardinality is inconsistent.
 */
void require_shape_result(const ShapeHistogramResult& result) {
    const std::size_t count = region_bin_count(result.region);
    const std::size_t core_count = region_bin_count(result.gaussian_core_region);
    if (!result.region.has_value()) {
        if (result.selected_count != 0U || !result.normalized_bins.empty() ||
            result.gaussian_core_region.has_value()) {
            throw std::logic_error(
                "HBT analysis state: empty shape result has derived bin data"
            );
        }
    } else if (result.selected_count != result.region->selected_count ||
               result.selected_count == 0U ||
               result.normalized_bins.size() != count) {
        throw std::logic_error(
            "HBT analysis state: shape region/count/normalization mismatch"
        );
    }

    if (result.gaussian.fully_valid) {
        if (!result.gaussian.radius.has_value() ||
            !result.gaussian.q_min.has_value() ||
            !result.gaussian_core_region.has_value() ||
            result.gaussian.fitted_pdf.size() != core_count) {
            throw std::logic_error(
                "HBT analysis state: valid Gaussian fit is incomplete"
            );
        }
    } else if (result.gaussian.radius.has_value() ||
               !result.gaussian.fitted_pdf.empty()) {
        throw std::logic_error(
            "HBT analysis state: invalid Gaussian fit has fabricated output"
        );
    }

    if (result.mixed.estimator != FitEstimator::Poisson ||
        result.mixed_neyman.estimator != FitEstimator::Neyman ||
        result.mixed_pearson.estimator != FitEstimator::Pearson) {
        throw std::logic_error(
            "HBT analysis state: mixed estimator identity mismatch"
        );
    }

    const std::array<const MixedFitResult*, 3U> mixed_results{
        &result.mixed,
        &result.mixed_neyman,
        &result.mixed_pearson
    };
    for (const MixedFitResult* mixed : mixed_results) {
        if (mixed->fully_valid) {
            if (!mixed->core_radius.has_value() ||
                !mixed->tail_radius.has_value() ||
                !mixed->core_fraction.has_value() ||
                !mixed->q_min.has_value() ||
                !mixed->selected_core_start.has_value() ||
                mixed->consensus_size < 4U ||
                mixed->fitted_pdf.size() != count) {
                throw std::logic_error(
                    "HBT analysis state: valid mixed fit is incomplete"
                );
            }
        } else if (mixed->core_radius.has_value() ||
                   mixed->tail_radius.has_value() ||
                   mixed->core_fraction.has_value() ||
                   !mixed->fitted_pdf.empty()) {
            throw std::logic_error(
                "HBT analysis state: invalid mixed fit has fabricated output"
            );
        }
    }
}

/**
 * @brief Validate one delta-t result's internal cardinality invariants.
 * @param result Derived delta-t result to inspect.
 * @throws std::logic_error If region/count/statistics are inconsistent.
 */
void require_delta_t_result(const DeltaTHistogramResult& result) {
    const std::size_t count = region_bin_count(result.region);
    if (!result.region.has_value()) {
        if (result.selected_count != 0U || !result.normalized_bins.empty() ||
            result.status != DeltaTStatisticsStatus::EmptyHistogram) {
            throw std::logic_error(
                "HBT analysis state: empty delta_t result is inconsistent"
            );
        }
        return;
    }

    if (result.selected_count != result.region->selected_count ||
        result.selected_count == 0U ||
        result.normalized_bins.size() != count) {
        throw std::logic_error(
            "HBT analysis state: delta_t region/count/normalization mismatch"
        );
    }
    if (result.status == DeltaTStatisticsStatus::Valid &&
        (!result.mean.has_value() || !result.sigma.has_value() ||
         !result.sigma_error.has_value())) {
        throw std::logic_error(
            "HBT analysis state: valid delta_t statistics are incomplete"
        );
    }
}

/**
 * @brief Validate all nine derived histograms of one destination.
 * @param derived Derived destination state to inspect.
 * @throws std::logic_error If any logical result is inconsistent.
 */
void require_derived_set(const DerivedHistogramSet& derived) {
    for (const ShapeHistogramResult& result : derived.osl) {
        require_shape_result(result);
    }
    for (const ShapeHistogramResult& result : derived.radial) {
        require_shape_result(result);
    }
    for (const DeltaTHistogramResult& result : derived.delta_t) {
        require_delta_t_result(result);
    }
}

}  // namespace

HistogramAnalysisState analyze_histograms(
    const HBTConfig& config,
    const RawHistogramState& raw
) {
    require_raw_histogram_state_layout(config, raw);

    HistogramAnalysisState derived;
    derived.products.resize(raw.products.size());

    for (std::size_t product = 0U; product < raw.products.size(); ++product) {
        const ProductRawHistogramState& raw_product = raw.products[product];
        ProductDerivedHistogramState& derived_product =
            derived.products[product];
        derived_product.origins.resize(raw_product.origins.size());

        for (std::size_t origin = 0U;
             origin < raw_product.origins.size();
             ++origin) {
            const OriginRawHistogramState& raw_origin =
                raw_product.origins[origin];
            OriginDerivedHistogramState& derived_origin =
                derived_product.origins[origin];
            derived_origin.global = analyze_histogram_set(
                config.histogram_config,
                raw_origin.global,
                false,
                false
            );
            derived_origin.slices.resize(raw_origin.slices.size());
            const bool apply_radial_mt_quality_cut =
                config.pair_slicing.mt.enabled &&
                !config.pair_slicing.kt.enabled;
            for (std::size_t slice = 0U;
                 slice < raw_origin.slices.size();
                 ++slice) {
                derived_origin.slices[slice] = analyze_histogram_set(
                    config.histogram_config,
                    raw_origin.slices[slice],
                    true,
                    apply_radial_mt_quality_cut
                );
            }
        }
    }

    for (std::size_t product = 0U; product < raw.products.size(); ++product) {
        for (std::size_t origin = 0U;
             origin < raw.products[product].origins.size();
             ++origin) {
            const OriginRawHistogramState& raw_origin =
                raw.products[product].origins[origin];
            OriginDerivedHistogramState& derived_origin =
                derived.products[product].origins[origin];
            normalize_histogram_set(
                config.histogram_config,
                raw_origin.global,
                derived_origin.global
            );
            for (std::size_t slice = 0U;
                 slice < raw_origin.slices.size();
                 ++slice) {
                normalize_histogram_set(
                    config.histogram_config,
                    raw_origin.slices[slice],
                    derived_origin.slices[slice]
                );
            }
        }
    }

    require_histogram_analysis_layout(config, raw, derived);
    return derived;
}

void require_histogram_analysis_layout(
    const HBTConfig& config,
    const RawHistogramState& raw,
    const HistogramAnalysisState& derived
) {
    require_raw_histogram_state_layout(config, raw);
    if (derived.products.size() != raw.products.size()) {
        throw std::logic_error(
            "HBT analysis state: product count differs from raw histogram state"
        );
    }

    for (std::size_t product = 0U; product < raw.products.size(); ++product) {
        const ProductRawHistogramState& raw_product = raw.products[product];
        const ProductDerivedHistogramState& derived_product =
            derived.products[product];
        if (derived_product.origins.size() != raw_product.origins.size()) {
            throw std::logic_error(
                "HBT analysis state: origin count differs from raw "
                "histogram state"
            );
        }
        for (std::size_t origin = 0U;
             origin < raw_product.origins.size();
             ++origin) {
            const OriginRawHistogramState& raw_origin =
                raw_product.origins[origin];
            const OriginDerivedHistogramState& derived_origin =
                derived_product.origins[origin];
            if (derived_origin.slices.size() != raw_origin.slices.size()) {
                throw std::logic_error(
                    "HBT analysis state: slice count differs from raw "
                    "histogram state"
                );
            }
            require_derived_set(derived_origin.global);
            for (const DerivedHistogramSet& slice : derived_origin.slices) {
                require_derived_set(slice);
            }
        }
    }
}

}  // namespace hbt
