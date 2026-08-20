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
    NotApplicable,           ///< Observable is intentionally not analyzed here.
    EmptyHistogram,          ///< No statistical region exists.
    InsufficientBins,        ///< K is smaller than P + 1.
    InsufficientStatistics,  ///< Radial mT slice is below the production cut.
    InvalidMomentSeed,       ///< Required data-derived seed is invalid.
    InvalidGaussianCoreAnchor, ///< Required Gaussian-core anchor is unavailable.
    NoBasinConsensus,        ///< Fewer than four core starts agree on a basin.
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
    bool valid_covariance;        ///< Returned user state has covariance.
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
 * @brief Final external coordinates used to identify one mixed-fit basin.
 *
 * Log-radius coordinates make numerical equivalence scale independent. This
 * structure contains no physical ordering information and is used only to
 * determine whether independent starts converged to the same solution.
 */
struct MixedBasinPoint {
    double log_core_radius;  ///< Final log(R_core).
    double log_tail_radius;  ///< Final log(R_tail).
    double core_fraction;    ///< Final f_core.
};

/** Numerical same-basin tolerance for each log-radius coordinate. */
constexpr double kMixedBasinLogRadiusTolerance = 0.01;

/** Numerical same-basin absolute tolerance for f_core. */
constexpr double kMixedBasinCoreFractionTolerance = 0.01;

/**
 * @brief Test whether two mixed endpoints represent the same numerical basin.
 * @param lhs First converged mixed endpoint.
 * @param rhs Second converged mixed endpoint.
 * @return true when both log radii and f_core agree within the documented
 *         production convergence tolerances.
 *
 * The comparison is deliberately independent of R_tail/R_core ordering and of
 * the Gaussian-core anchor. The 0.01 tolerances identify repeated numerical
 * convergence and are much smaller than the separated basins observed in the
 * validation study.
 */
[[nodiscard]] bool same_mixed_basin(
    const MixedBasinPoint& lhs,
    const MixedBasinPoint& rhs
);

/**
 * @brief Return the largest connected same-basin group among valid starts.
 * @param endpoints Final endpoint for every attempted deterministic start.
 * @param valid_indices Start indices whose MIGRAD states are publishable.
 * @return Start indices in the largest connected numerical basin.
 * @throws std::out_of_range If a valid start index has no matching endpoint.
 *
 * Pairwise same-basin agreement defines an undirected graph. Connected
 * components make grouping deterministic even when tiny numerical differences
 * are not perfectly transitive. Ties in component size retain the component
 * reached from the lowest valid start index.
 */
[[nodiscard]] std::vector<std::size_t> largest_mixed_basin_group(
    const std::vector<MixedBasinPoint>& endpoints,
    const std::vector<std::size_t>& valid_indices
);

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
 * Five deterministic core-anchored starts are attempted when the Gaussian
 * core anchor and moment-derived tail seed are available. Every start uses
 * the same R_core seed from the independently fitted Gaussian core and varies
 * only R_tail and f_core. `consensus_size` records the largest numerically
 * equivalent solution group. A physical mixed result is published only when
 * at least four of the five starts agree on the same basin. Q is used only to
 * select the best numerical realization inside that consensus basin.
 */
struct MixedFitResult {
    /// Number of deterministic Gaussian-core-anchored MIGRAD starts.
    static constexpr std::size_t kCoreStartCount = 5U;

    bool fully_valid;                 ///< Consensus fit and all MINOS are valid.
    FitFailureReason failure_reason;  ///< Primary invalidity cause, if any.
    /// Diagnostics for the five deterministic Gaussian-core-anchored starts.
    std::array<MigradDiagnostic, kCoreStartCount> starts;
    std::size_t starts_attempted;     ///< Number of core starts actually run.
    std::size_t valid_starts;         ///< Starts with valid evaluable MIGRAD state.
    std::size_t consensus_size;       ///< Starts assigned to selected basin.
    /// Zero-based member of the consensus basin selected for MINOS.
    std::optional<std::size_t> selected_core_start;
    MigradDiagnostic selected_migrad; ///< Selected consensus MIGRAD diagnostic.
    MinosDiagnostic minos_core_radius; ///< MINOS diagnostic for log(R_core).
    MinosDiagnostic minos_tail_radius; ///< MINOS diagnostic for log(R_tail).
    MinosDiagnostic minos_core_fraction; ///< MINOS diagnostic for f_core.
    std::optional<double> q_min;       ///< Minimum Q inside consensus basin.
    /// Physical core radius and asymmetric MINOS errors when fully valid.
    std::optional<FitParameterEstimate> core_radius;
    /// Physical tail radius and asymmetric MINOS errors when fully valid.
    std::optional<FitParameterEstimate> tail_radius;
    /// Physical core fraction and asymmetric MINOS errors when fully valid.
    std::optional<FitParameterEstimate> core_fraction;
    /// Valid fitted bin densities over the full mixed region.
    std::vector<double> fitted_pdf;
};

/**
 * @brief Complete derived result for one logical OSL or radial histogram.
 */
struct ShapeHistogramResult {
    /// Full contiguous region retained for normalization and the mixed model.
    std::optional<StatisticalRegion> region;
    /// Compact region used only by the pure Gaussian core fit.
    std::optional<StatisticalRegion> gaussian_core_region;
    std::uint64_t selected_count;            ///< N over the full region.
    /// Final normalized distribution over the full statistical region.
    std::vector<NormalizedHistogramBin> normalized_bins;
    GaussianFitResult gaussian;              ///< Truncated pure Gaussian core fit.
    MixedFitResult mixed;                    ///< Full-range mixed-model fit.
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
