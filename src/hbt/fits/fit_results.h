/**
 * @file fit_results.h
 * @brief Derived post-sample HBT statistical result types.
 */

#ifndef HBT_FITS_FIT_RESULTS_H
#define HBT_FITS_FIT_RESULTS_H

#include <array>
#include <cstddef>
#include <cstdint>
#include <optional>
#include <vector>

namespace hbt {

/**
 * @brief OSL-like or radial model family used by one fit.
 */
enum class FitObservableFamily {
    OSL,    ///< Absolute OSL marginal.
    Radial  ///< Three-dimensional radial separation.
};

/**
 * @brief Stable reason describing why a fit is not fully valid.
 */
enum class FitFailureReason {
    None,                    ///< Fit is fully valid.
    EmptyHistogram,          ///< No statistical region exists.
    InsufficientBins,        ///< K is smaller than P + 1.
    InvalidMomentSeed,       ///< Required data-derived seed is invalid.
    ObjectiveEvaluation,     ///< Objective could not be evaluated.
    MigradInvalid,           ///< MIGRAD returned an invalid minimum.
    MigradCallLimit,         ///< MIGRAD exhausted its call budget.
    MigradAboveMaxEdm,       ///< MIGRAD stopped above maximum EDM.
    NonFiniteMinimum,        ///< Minimum or physical parameter is non-finite.
    DegenerateCoreFraction,  ///< Mixed fit ended at f_core == 0 or 1.
    MinosLowerInvalid,       ///< A required lower MINOS side is invalid.
    MinosUpperInvalid,       ///< A required upper MINOS side is invalid.
    MinosLowerCallLimit,     ///< Lower MINOS side exhausted its call budget.
    MinosUpperCallLimit,     ///< Upper MINOS side exhausted its call budget.
    MinosLowerLimit,         ///< Lower MINOS side reached a parameter limit.
    MinosUpperLimit,         ///< Upper MINOS side reached a parameter limit.
    MinosLowerNewMinimum,    ///< Lower MINOS side found a new minimum.
    MinosUpperNewMinimum     ///< Upper MINOS side found a new minimum.
};

/**
 * @brief Contiguous selected range in one raw histogram.
 *
 * first_bin and last_bin are inclusive zero-based indices in the owning raw
 * histogram. selected_count is the exact sum of raw uint64_t counts over that
 * interval. The structure owns no raw counts and retains no references.
 */
struct StatisticalRegion {
    std::size_t first_bin;          ///< Inclusive first selected bin.
    std::size_t last_bin;           ///< Inclusive last selected bin.
    std::uint64_t selected_count;   ///< Sum of selected raw counts.
};

/**
 * @brief One normalized presentation bin produced after statistical analysis.
 */
struct NormalizedHistogramBin {
    std::size_t bin_index;     ///< Original zero-based raw-histogram index.
    double lower_edge;         ///< Exact configured lower edge.
    double upper_edge;         ///< Exact configured upper edge.
    double center;             ///< Arithmetic center of the uniform bin.
    double pdf;                ///< n_i / (N * Delta_x).
    double counting_error_pdf; ///< sqrt(n_i) / (N * Delta_x).
};

/**
 * @brief Diagnostic status returned by one independent MIGRAD start.
 */
struct MigradDiagnostic {
    bool attempted;               ///< Whether this start executed MIGRAD.
    bool valid;                   ///< Raw FunctionMinimum::IsValid state.
    bool reached_call_limit;      ///< MIGRAD call budget was exhausted.
    bool above_max_edm;           ///< Minimum remained above maximum EDM.
    bool objective_failure;       ///< Returned minimum cannot be evaluated.
    int function_calls;           ///< Function evaluations reported by Minuit.
    std::optional<double> q_min;  ///< Finite objective value when available.
};

/**
 * @brief Diagnostic status for both MINOS sides of one fitted parameter.
 */
struct MinosDiagnostic {
    /// Whether MINOS was invoked for this parameter.
    bool attempted;
    bool lower_valid;          ///< Lower MINOS crossing is valid.
    bool upper_valid;          ///< Upper MINOS crossing is valid.
    bool at_lower_limit;       ///< Lower side reached a parameter limit.
    bool at_upper_limit;       ///< Upper side reached a parameter limit.
    bool lower_call_limit;     ///< Lower side exhausted its call budget.
    bool upper_call_limit;     ///< Upper side exhausted its call budget.
    bool lower_new_minimum;    ///< Lower side found a new minimum.
    bool upper_new_minimum;    ///< Upper side found a new minimum.
};

/**
 * @brief Physical parameter value with asymmetric MINOS errors.
 *
 * lower_error and upper_error are non-negative distances from value in the
 * physical parameter space. For log-radius fits the transformation from the
 * MINOS interval is completed before this structure is created.
 */
struct FitParameterEstimate {
    double value;        ///< Physical fitted parameter value.
    double lower_error;  ///< Distance to the lower MINOS endpoint.
    double upper_error;  ///< Distance to the upper MINOS endpoint.
};

/**
 * @brief Classify the primary failure represented by one MIGRAD diagnostic.
 * @param diagnostic Completed stable MIGRAD diagnostic.
 * @return FitFailureReason::None when the minimum is fully valid.
 *
 * Objective evaluation failure has precedence over Minuit convergence states.
 * The classifier is pure and retains no Minuit object.
 */
[[nodiscard]] FitFailureReason fit_failure_from_migrad(
    const MigradDiagnostic& diagnostic
);

/**
 * @brief Classify the primary failure represented by one MINOS diagnostic.
 * @param diagnostic Completed stable two-sided MINOS diagnostic.
 * @param reject_limits Whether parameter-limit crossings invalidate the fit.
 * @return FitFailureReason::None when both required sides are valid.
 *
 * Limit rejection is enabled for the bounded core fraction and disabled for
 * log-radius parameters, which have no physical bounds.
 */
[[nodiscard]] FitFailureReason fit_failure_from_minos(
    const MinosDiagnostic& diagnostic,
    bool reject_limits
);

/**
 * @brief Classify the physical validity of a fitted mixed core fraction.
 * @param core_fraction Fitted physical fraction.
 * @return None for 0 < f_core < 1, DegenerateCoreFraction at either exact
 *         endpoint, or NonFiniteMinimum for a non-finite/out-of-domain value.
 */
[[nodiscard]] FitFailureReason mixed_core_fraction_failure(
    double core_fraction
);

/**
 * @brief Result of the one-parameter pure-Gaussian fit.
 */
struct GaussianFitResult {
    bool fully_valid;                 ///< MIGRAD and required MINOS are valid.
    FitFailureReason failure_reason;  ///< Primary invalidity cause, if any.
    MigradDiagnostic migrad;          ///< MIGRAD diagnostic for the only start.
    MinosDiagnostic minos_radius;     ///< MINOS diagnostic for log(R).
    std::optional<double> q_min;      ///< Selected minimum objective value.
    /// Physical R and asymmetric MINOS errors when fully valid.
    std::optional<FitParameterEstimate> radius;
    /// Valid fitted bin densities; empty when the fit is not fully valid.
    std::vector<double> fitted_pdf;
};

/**
 * @brief Result of the three-parameter Gaussian-plus-exponential fit.
 *
 * starts stores deterministic start A at index zero and optional start B at
 * index one. starts_attempted is one or two. selected_start is present only
 * when at least one valid MIGRAD minimum exists. exact_q_tie records equality
 * of two valid objective minima without asserting global optimality.
 */
struct MixedFitResult {
    bool fully_valid;                 ///< Selected fit and all MINOS are valid.
    FitFailureReason failure_reason;  ///< Primary invalidity cause, if any.
    std::array<MigradDiagnostic, 2U> starts; ///< Independent start diagnostics.
    std::size_t starts_attempted;     ///< Number of deterministic starts used.
    std::optional<std::size_t> selected_start; ///< Zero-based selected start.
    bool exact_q_tie;                 ///< Valid starts had exactly equal Q_min.
    bool tail_below_core;             ///< R_tail < R_core diagnostic.
    MigradDiagnostic selected_migrad; ///< Selected MIGRAD diagnostic.
    MinosDiagnostic minos_core_radius; ///< MINOS diagnostic for log(R_core).
    MinosDiagnostic minos_tail_radius; ///< MINOS diagnostic for log(R_tail).
    MinosDiagnostic minos_core_fraction; ///< MINOS diagnostic for f_core.
    std::optional<double> q_min;       ///< Best valid minimum found.
    /// Physical core radius and asymmetric MINOS errors when fully valid.
    std::optional<FitParameterEstimate> core_radius;
    /// Physical tail radius and asymmetric MINOS errors when fully valid.
    std::optional<FitParameterEstimate> tail_radius;
    /// Physical core fraction and asymmetric MINOS errors when fully valid.
    std::optional<FitParameterEstimate> core_fraction;
    /// Valid fitted bin densities; empty when the fit is not fully valid.
    std::vector<double> fitted_pdf;
};

/**
 * @brief Complete derived result for one logical OSL or radial histogram.
 */
struct ShapeHistogramResult {
    std::optional<StatisticalRegion> region; ///< Selected contiguous region.
    std::uint64_t selected_count;            ///< N over the selected region.
    /// Final normalized distribution, produced after both fit attempts.
    std::vector<NormalizedHistogramBin> normalized_bins;
    GaussianFitResult gaussian;              ///< Independent pure Gaussian fit.
    MixedFitResult mixed;                    ///< Independent mixed-model fit.
};

/**
 * @brief Diagnostic state for signed delta-t moment calculation.
 */
enum class DeltaTStatisticsStatus {
    Valid,             ///< Mean, sigma, and sigma error are valid.
    EmptyHistogram,    ///< No selected contiguous region exists.
    InsufficientCount, ///< N <= 1 prevents the required sigma error.
    InvalidVariance    ///< Computed variance is negative or non-finite.
};

/**
 * @brief Complete derived result for one signed delta-t histogram.
 */
struct DeltaTHistogramResult {
    std::optional<StatisticalRegion> region; ///< Selected peak-centered region.
    std::uint64_t selected_count;            ///< N over the selected region.
    /// Final normalized distribution, produced after moment calculation.
    std::vector<NormalizedHistogramBin> normalized_bins;
    /// Explicit moment validity state.
    DeltaTStatisticsStatus status;
    std::optional<double> mean;               ///< Weighted mean when defined.
    /// Population sigma when defined.
    std::optional<double> sigma;
    /// Required sigma error when valid.
    std::optional<double> sigma_error;
};

/**
 * @brief Derived nine-histogram state for one raw histogram destination.
 */
struct DerivedHistogramSet {
    std::array<ShapeHistogramResult, 4U> osl;     ///< Four OSL marginals.
    std::array<ShapeHistogramResult, 2U> radial;  ///< Two radial histograms.
    std::array<DeltaTHistogramResult, 3U> delta_t; ///< Three signed times.
};

/**
 * @brief Global and optional kinetic-slice post-sample state for one origin.
 */
struct OriginDerivedHistogramState {
    DerivedHistogramSet global;              ///< Global destination result.
    std::vector<DerivedHistogramSet> slices; ///< Flat-slice results.
};

/**
 * @brief post-sample origin results aligned with one configured final product.
 */
struct ProductDerivedHistogramState {
    /// Results aligned exactly with raw-histogram origin state order.
    std::vector<OriginDerivedHistogramState> origins;
};

/**
 * @brief Complete post-sample state aligned with raw products and origins.
 *
 * The state stores no raw uint64_t counts and no duplicate identity strings.
 * It is safe to retain after the synchronous post-sample analysis returns.
 */
struct HistogramAnalysisState {
    /// Final-product results aligned with RawHistogramState::products.
    std::vector<ProductDerivedHistogramState> products;
};

/**
 * @brief Return a stable ASCII token for one fit failure reason.
 * @param reason Failure reason to serialize or diagnose.
 * @return Static token whose lifetime is the complete program.
 * @throws std::invalid_argument If @p reason is not a valid enum value.
 */
[[nodiscard]] const char* fit_failure_reason_token(FitFailureReason reason);

/**
 * @brief Return a stable ASCII token for one delta-t status.
 * @param status Delta-t status to serialize or diagnose.
 * @return Static token whose lifetime is the complete program.
 * @throws std::invalid_argument If @p status is not a valid enum value.
 */
[[nodiscard]] const char* delta_t_status_token(
    DeltaTStatisticsStatus status
);

}  // namespace hbt

#endif  // HBT_FITS_FIT_RESULTS_H
