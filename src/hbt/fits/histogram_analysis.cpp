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
 * @param estimator Objective identity retained even when no fit is attempted.
 * @return Invalid result carrying only diagnostics, @p reason, and @p estimator.
 */
GaussianFitResult invalid_gaussian_result(
    FitFailureReason reason,
    FitEstimator estimator
) {
    const MigradDiagnostic empty = unattempted_migrad();
    std::array<MigradDiagnostic, GaussianFitResult::kStartCount> starts{};
    starts.fill(empty);
    return {
        false,
        reason,
        estimator,
        starts,
        0U,
        0U,
        std::nullopt,
        empty,
        unattempted_minos(),
        unattempted_minos(),
        std::nullopt,
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
    std::array<MigradDiagnostic, MixedFitResult::kCoreStartCount> starts{};
    starts.fill(empty);
    std::array<
        MixedStartEndpointDiagnostic,
        MixedFitResult::kCoreStartCount
    > start_endpoints{};
    return {
        false,
        reason,
        estimator,
        starts,
        start_endpoints,
        0U,
        0U,
        0U,
        std::nullopt,
        empty,
        empty,
        unattempted_minos(),
        unattempted_minos(),
        unattempted_minos(),
        unattempted_minos(),
        std::nullopt,
        std::nullopt,
        std::nullopt,
        std::nullopt,
        std::nullopt,
        {}
    };
}

/**
 * @brief Return whether one estimator is enabled by the configured mode.
 */
bool fit_estimator_enabled(
    FitEstimatorMode mode,
    FitEstimator estimator
) {
    switch (mode) {
        case FitEstimatorMode::All:
            return true;
        case FitEstimatorMode::Poisson:
            return estimator == FitEstimator::Poisson;
        case FitEstimatorMode::Neyman:
            return estimator == FitEstimator::Neyman;
        case FitEstimatorMode::Pearson:
            return estimator == FitEstimator::Pearson;
    }
    throw std::logic_error("HBT analysis: unknown fit estimator mode");
}

/**
 * @brief Build a Gaussian failure or an explicit skipped-estimator result.
 */
GaussianFitResult configured_gaussian_failure(
    FitEstimatorMode mode,
    FitFailureReason reason,
    FitEstimator estimator
) {
    return invalid_gaussian_result(
        fit_estimator_enabled(mode, estimator)
            ? reason
            : FitFailureReason::NotApplicable,
        estimator
    );
}

/**
 * @brief Build a mixed failure or an explicit skipped-estimator result.
 */
MixedFitResult configured_mixed_failure(
    FitEstimatorMode mode,
    FitFailureReason reason,
    FitEstimator estimator
) {
    return invalid_mixed_result(
        fit_estimator_enabled(mode, estimator)
            ? reason
            : FitFailureReason::NotApplicable,
        estimator
    );
}

/**
 * @brief Build a derived shape placeholder for an intentionally skipped fit.
 * @return Empty, explicitly not-applicable Gaussian and mixed results.
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
        invalid_gaussian_result(
            FitFailureReason::NotApplicable, FitEstimator::Poisson
        ),
        invalid_gaussian_result(
            FitFailureReason::NotApplicable, FitEstimator::Neyman
        ),
        invalid_gaussian_result(
            FitFailureReason::NotApplicable, FitEstimator::Pearson
        ),
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

/** Radial mT-slice threshold; zero keeps the quality-cut machinery disabled. */
constexpr std::uint64_t kMinimumRadialSliceSelectedCount = 0U;

/**
 * @brief Analyze one OSL or radial logical histogram without normalization.
 * @param family OSL or radial physical model family.
 * @param bins Slot-major raw uint64_t histogram storage.
 * @param offset First raw counter for this logical histogram.
 * @param binning Validated uniform binning.
 * @param apply_radial_slice_quality_cut Whether the configured radial kinetic
 *        slice threshold is evaluated for this histogram.
 * @param estimator_mode Configured estimator set to execute. Disabled
 *        estimators are returned explicitly as NotApplicable and are never
 *        passed to MIGRAD or MINOS.
 * @return Full/core regions plus configured independent Gaussian and mixed
 *         fits, with stable NotApplicable placeholders for skipped estimators.
 *
 * N_selected and both statistical regions are established before any fit and
 * are never modified by estimator choice. The half-maximum radius seed is
 * derived once from the full selected histogram and shared numerically across
 * estimators. Each mixed estimator receives R_G only from its matching pure
 * Gaussian estimator.
 */
ShapeHistogramResult analyze_shape_histogram(
    FitObservableFamily family,
    const std::vector<std::uint64_t>& bins,
    std::size_t offset,
    const HistogramBinningConfig& binning,
    bool apply_radial_slice_quality_cut,
    FitEstimatorMode estimator_mode
) {
    const std::optional<StatisticalRegion> region =
        select_shape_region(bins, offset, binning);
    if (!region.has_value()) {
        return {
            std::nullopt,
            std::nullopt,
            0U,
            {},
            configured_gaussian_failure(
                estimator_mode,
                FitFailureReason::EmptyHistogram,
                FitEstimator::Poisson
            ),
            configured_gaussian_failure(
                estimator_mode,
                FitFailureReason::EmptyHistogram,
                FitEstimator::Neyman
            ),
            configured_gaussian_failure(
                estimator_mode,
                FitFailureReason::EmptyHistogram,
                FitEstimator::Pearson
            ),
            configured_mixed_failure(
                estimator_mode,
                FitFailureReason::EmptyHistogram,
                FitEstimator::Poisson
            ),
            configured_mixed_failure(
                estimator_mode,
                FitFailureReason::EmptyHistogram,
                FitEstimator::Neyman
            ),
            configured_mixed_failure(
                estimator_mode,
                FitFailureReason::EmptyHistogram,
                FitEstimator::Pearson
            )
        };
    }

    // Pure Gaussian fits use the complete shape region. Their free positive
    // amplitude lets the Gaussian vanish naturally outside the core without
    // an externally imposed 5%/10% threshold.
    const std::optional<StatisticalRegion> core_region = region;
    if (apply_radial_slice_quality_cut &&
        region->selected_count < kMinimumRadialSliceSelectedCount) {
        return {
            region,
            core_region,
            region->selected_count,
            {},
            configured_gaussian_failure(
                estimator_mode,
                FitFailureReason::InsufficientStatistics,
                FitEstimator::Poisson
            ),
            configured_gaussian_failure(
                estimator_mode,
                FitFailureReason::InsufficientStatistics,
                FitEstimator::Neyman
            ),
            configured_gaussian_failure(
                estimator_mode,
                FitFailureReason::InsufficientStatistics,
                FitEstimator::Pearson
            ),
            configured_mixed_failure(
                estimator_mode,
                FitFailureReason::InsufficientStatistics,
                FitEstimator::Poisson
            ),
            configured_mixed_failure(
                estimator_mode,
                FitFailureReason::InsufficientStatistics,
                FitEstimator::Neyman
            ),
            configured_mixed_failure(
                estimator_mode,
                FitFailureReason::InsufficientStatistics,
                FitEstimator::Pearson
            )
        };
    }
    if (!core_region.has_value()) {
        return {
            region,
            std::nullopt,
            region->selected_count,
            {},
            configured_gaussian_failure(
                estimator_mode,
                FitFailureReason::InsufficientBins,
                FitEstimator::Poisson
            ),
            configured_gaussian_failure(
                estimator_mode,
                FitFailureReason::InsufficientBins,
                FitEstimator::Neyman
            ),
            configured_gaussian_failure(
                estimator_mode,
                FitFailureReason::InsufficientBins,
                FitEstimator::Pearson
            ),
            configured_mixed_failure(
                estimator_mode,
                FitFailureReason::InvalidGaussianCoreAnchor,
                FitEstimator::Poisson
            ),
            configured_mixed_failure(
                estimator_mode,
                FitFailureReason::InvalidGaussianCoreAnchor,
                FitEstimator::Neyman
            ),
            configured_mixed_failure(
                estimator_mode,
                FitFailureReason::InvalidGaussianCoreAnchor,
                FitEstimator::Pearson
            )
        };
    }

    const std::optional<double> half_maximum_seed = half_maximum_radius_seed(
        family,
        bins,
        offset,
        binning,
        region.value()
    );
    if (!half_maximum_seed.has_value()) {
        return {
            region,
            core_region,
            region->selected_count,
            {},
            configured_gaussian_failure(
                estimator_mode,
                FitFailureReason::InvalidHalfMaximumSeed,
                FitEstimator::Poisson
            ),
            configured_gaussian_failure(
                estimator_mode,
                FitFailureReason::InvalidHalfMaximumSeed,
                FitEstimator::Neyman
            ),
            configured_gaussian_failure(
                estimator_mode,
                FitFailureReason::InvalidHalfMaximumSeed,
                FitEstimator::Pearson
            ),
            configured_mixed_failure(
                estimator_mode,
                FitFailureReason::InvalidHalfMaximumSeed,
                FitEstimator::Poisson
            ),
            configured_mixed_failure(
                estimator_mode,
                FitFailureReason::InvalidHalfMaximumSeed,
                FitEstimator::Neyman
            ),
            configured_mixed_failure(
                estimator_mode,
                FitFailureReason::InvalidHalfMaximumSeed,
                FitEstimator::Pearson
            )
        };
    }

    GaussianFitResult gaussian = fit_estimator_enabled(
        estimator_mode,
        FitEstimator::Poisson
    ) ? fit_gaussian_model(
        family,
        bins,
        offset,
        binning,
        core_region.value(),
        FitEstimator::Poisson,
        half_maximum_seed.value()
    ) : invalid_gaussian_result(
        FitFailureReason::NotApplicable,
        FitEstimator::Poisson
    );
    GaussianFitResult gaussian_neyman = fit_estimator_enabled(
        estimator_mode,
        FitEstimator::Neyman
    ) ? fit_gaussian_model(
        family,
        bins,
        offset,
        binning,
        core_region.value(),
        FitEstimator::Neyman,
        half_maximum_seed.value()
    ) : invalid_gaussian_result(
        FitFailureReason::NotApplicable,
        FitEstimator::Neyman
    );
    GaussianFitResult gaussian_pearson = fit_estimator_enabled(
        estimator_mode,
        FitEstimator::Pearson
    ) ? fit_gaussian_model(
        family,
        bins,
        offset,
        binning,
        core_region.value(),
        FitEstimator::Pearson,
        half_maximum_seed.value()
    ) : invalid_gaussian_result(
        FitFailureReason::NotApplicable,
        FitEstimator::Pearson
    );

    MixedFitResult mixed = fit_estimator_enabled(
        estimator_mode,
        FitEstimator::Poisson
    ) ? fit_mixed_model(
        family,
        bins,
        offset,
        binning,
        region.value(),
        FitEstimator::Poisson,
        gaussian,
        half_maximum_seed.value()
    ) : invalid_mixed_result(
        FitFailureReason::NotApplicable,
        FitEstimator::Poisson
    );
    MixedFitResult mixed_neyman = fit_estimator_enabled(
        estimator_mode,
        FitEstimator::Neyman
    ) ? fit_mixed_model(
        family,
        bins,
        offset,
        binning,
        region.value(),
        FitEstimator::Neyman,
        gaussian_neyman,
        half_maximum_seed.value()
    ) : invalid_mixed_result(
        FitFailureReason::NotApplicable,
        FitEstimator::Neyman
    );
    MixedFitResult mixed_pearson = fit_estimator_enabled(
        estimator_mode,
        FitEstimator::Pearson
    ) ? fit_mixed_model(
        family,
        bins,
        offset,
        binning,
        region.value(),
        FitEstimator::Pearson,
        gaussian_pearson,
        half_maximum_seed.value()
    ) : invalid_mixed_result(
        FitFailureReason::NotApplicable,
        FitEstimator::Pearson
    );
    return {
        region,
        core_region,
        region->selected_count,
        {},
        std::move(gaussian),
        std::move(gaussian_neyman),
        std::move(gaussian_pearson),
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
 * @param apply_radial_mt_quality_cut Whether the radial N_selected threshold
 *        machinery applies to this one-dimensional radial mT slice. The
 *        current threshold is zero, so no non-empty slice is vetoed by it.
 * @param estimator_mode Configured estimator set propagated to every OSL and
 *        radial fit in this destination.
 * @return Derived state with configured fits/moments complete and
 *         normalization empty.
 */
DerivedHistogramSet analyze_histogram_set(
    const HBTHistogramConfig& config,
    const RawHistogramSet& raw,
    bool kinetic_slice,
    bool apply_radial_mt_quality_cut,
    FitEstimatorMode estimator_mode
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
                false,
                estimator_mode
            );
    }
    for (std::size_t slot = 0U; slot < result.radial.size(); ++slot) {
        result.radial[slot] = analyze_shape_histogram(
            FitObservableFamily::Radial,
            raw.radial.bins,
            slot * config.radial.nbins,
            config.radial,
            apply_radial_mt_quality_cut,
            estimator_mode
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

    if (result.gaussian.estimator != FitEstimator::Poisson ||
        result.gaussian_neyman.estimator != FitEstimator::Neyman ||
        result.gaussian_pearson.estimator != FitEstimator::Pearson) {
        throw std::logic_error(
            "HBT analysis state: Gaussian estimator identity mismatch"
        );
    }
    const std::array<const GaussianFitResult*, 3U> gaussian_results{
        &result.gaussian,
        &result.gaussian_neyman,
        &result.gaussian_pearson
    };
    for (const GaussianFitResult* gaussian : gaussian_results) {
        if (gaussian->fully_valid) {
            if (!gaussian->radius.has_value() ||
                !gaussian->amplitude.has_value() ||
                !gaussian->q_min.has_value() ||
                !gaussian->selected_start.has_value() ||
                !result.gaussian_core_region.has_value() ||
                gaussian->fitted_pdf.size() != core_count) {
                throw std::logic_error(
                    "HBT analysis state: valid Gaussian fit is incomplete"
                );
            }
        } else if (gaussian->radius.has_value() ||
                   gaussian->amplitude.has_value() ||
                   !gaussian->fitted_pdf.empty()) {
            throw std::logic_error(
                "HBT analysis state: invalid Gaussian fit has fabricated output"
            );
        }
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
                !mixed->amplitude.has_value() ||
                !mixed->q_min.has_value() ||
                !mixed->selected_core_start.has_value() ||
                fit_failure_from_migrad(mixed->amplitude_profile_migrad) !=
                    FitFailureReason::None ||
                mixed->fitted_pdf.size() != count) {
                throw std::logic_error(
                    "HBT analysis state: valid mixed fit is incomplete"
                );
            }
        } else if (mixed->core_radius.has_value() ||
                   mixed->tail_radius.has_value() ||
                   mixed->core_fraction.has_value() ||
                   mixed->amplitude.has_value() ||
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
        (!result.mean.has_value() || !result.mean_error.has_value() ||
         !result.sigma.has_value() || !result.sigma_error.has_value())) {
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
                false,
                config.fit_estimator_mode
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
                    apply_radial_mt_quality_cut,
                    config.fit_estimator_mode
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
