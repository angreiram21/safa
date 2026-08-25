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
 * @brief Statistical objective used to estimate Gaussian and mixed parameters.
 *
 * Poisson, Neyman, and Pearson are fitted independently. Each estimator keeps
 * its own deterministic starts, selected minimum, MIGRAD state, and MINOS
 * uncertainties; numerical results are never pooled across estimators.
 */
enum class FitEstimator {
    Poisson, ///< Binned Poisson deviance; production/default estimator.
    Neyman,  ///< Neyman chi-square, omitting bins with n_i == 0.
    Pearson  ///< Pearson chi-square with expected counts in the denominator.
};

/**
 * @brief Origin-dependent physical admissibility policy for mixed-fit basins.
 *
 * PRD requires both Gaussian core and exponential tail to remain appreciable.
 * P and PR also require a genuinely mixed solution and reject the degenerate
 * near-pure-Gaussian limit.
 */
enum class MixedCoreFractionPolicy {
    RequireCoreAndTail,   ///< PRD: require 0.1 < mean(f_core) < 0.9.
    RejectPureGaussian   ///< P/PR: require 0.1 < mean(f_core) < 0.99.
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
    InvalidMomentSeed,       ///< Required data-derived moment seed is invalid.
    InvalidHalfMaximumSeed,  ///< Required half-maximum radius seed is invalid.
    InvalidGaussianCoreAnchor, ///< Required Gaussian-core anchor is unavailable.
    NoBasinConsensus,        ///< Legacy diagnostic token; no longer an acceptance veto.
    ObjectiveEvaluation,     ///< Objective could not be evaluated.
    MigradInvalid,           ///< MIGRAD returned an invalid minimum.
    MigradCallLimit,         ///< MIGRAD exhausted its call budget.
    MigradAboveMaxEdm,       ///< MIGRAD stopped above maximum EDM.
    NonFiniteMinimum,        ///< Minimum or physical parameter is non-finite.
    DegenerateCoreFraction,  ///< No basin passes the origin-specific physical f_core policy.
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
    std::optional<double> q_min;  ///< Finite estimator objective when available.
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

/**
 * @brief Final physical endpoint retained for one deterministic mixed start.
 *
 * Values are recorded independently of whether the corresponding MIGRAD start
 * is ultimately accepted. A field is empty only when the terminal Minuit state
 * cannot be converted to a finite physical value. These diagnostics are used
 * solely to inspect the geometry of the 36-start basin search and do not alter
 * fit selection or validity.
 */
struct MixedStartEndpointDiagnostic {
    std::optional<double> core_radius;   ///< Final physical R_core when finite.
    std::optional<double> tail_radius;   ///< Final physical R_tail when finite.
    std::optional<double> core_fraction; ///< Final physical f_core when finite.
};

/** Numerical same-basin tolerance for each log-radius coordinate. */
constexpr double kMixedBasinLogRadiusTolerance = 0.01;

/** Numerical same-basin absolute tolerance for f_core. */
constexpr double kMixedBasinCoreFractionTolerance = 0.01;

/** Lower exclusive f_core bound for a physically non-degenerate mixed basin. */
constexpr double kMixedPhysicalCoreFractionMin = 0.1;

/** PRD-only upper exclusive f_core bound for a non-degenerate mixed basin. */
constexpr double kMixedPhysicalCoreFractionMax = 0.9;

/** P/PR upper exclusive f_core bound rejecting the Gaussian-limit degeneracy. */
constexpr double kMixedPPrCoreFractionMax = 0.99;

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
 * @brief Select the mixed start from the largest origin-admissible basin.
 * @param endpoints Final endpoint for every attempted deterministic start.
 * @param q_values Objective value associated with every attempted start.
 * @param valid_indices Start indices whose MIGRAD states are numerically valid.
 * @param core_fraction_policy Origin-dependent f_core admissibility policy.
 * @return Selected start index, or std::nullopt when every numerical basin is
 *         physically degenerate.
 * @throws std::invalid_argument If array sizes differ, no valid start is
 *         supplied, or a valid q is non-finite.
 * @throws std::out_of_range If a valid start index is outside the input arrays.
 *
 * Valid endpoints are partitioned into connected numerical basins using
 * same_mixed_basin(). Basin admissibility is determined first from the
 * arithmetic mean f_core. PRD retains the strict interval
 * 0.1 < mean(f_core) < 0.9. P and PR use
 * 0.1 < mean(f_core) < 0.99, rejecting the near-pure-Gaussian degeneracy.
 *
 * Among admissible basins, every origin selects the basin reached by the
 * largest number of converged deterministic starts. Equal-size basins are
 * ranked by their smallest finite q, then by the lowest start index. q is
 * never divided or otherwise normalized by basin multiplicity. Once the basin
 * is fixed, the lowest-q start in that basin is selected. R_HM remains part of
 * the deterministic seed set but does not rank final basins. No ordering
 * between R_core and R_tail is imposed.
 */
[[nodiscard]] std::optional<std::size_t> select_mixed_start_by_largest_basin(
    const std::vector<MixedBasinPoint>& endpoints,
    const std::vector<double>& q_values,
    const std::vector<std::size_t>& valid_indices,
    MixedCoreFractionPolicy core_fraction_policy
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
 * @brief Result of one independent two-parameter pure-Gaussian fit.
 *
 * Two deterministic radius starts are attempted when available: the moment
 * seed and the half-maximum-width seed converted to the model radius. The
 * numerically valid MIGRAD solution with the smallest estimator objective is
 * selected for MINOS. The Gaussian is fitted over the full contiguous shape
 * region. The pure Gaussian has a free positive
 * amplitude and therefore is not forced to carry unit area over that region.
 */
struct GaussianFitResult {
    /** Number of deterministic Gaussian MIGRAD starts. */
    static constexpr std::size_t kStartCount = 2U;

    bool fully_valid;                 ///< Selected MIGRAD and MINOS are valid.
    FitFailureReason failure_reason;  ///< Primary invalidity cause, if any.
    FitEstimator estimator;           ///< Objective minimized by this fit only.
    /// Diagnostics for moment and half-maximum starts, respectively.
    std::array<MigradDiagnostic, kStartCount> starts;
    std::size_t starts_attempted;      ///< Number of Gaussian starts actually run.
    std::size_t valid_starts;          ///< Starts with valid evaluable MIGRAD state.
    /// Zero-based start whose valid objective is the smallest found.
    std::optional<std::size_t> selected_start;
    MigradDiagnostic migrad;           ///< Selected minimum MIGRAD diagnostic.
    MinosDiagnostic minos_radius;      ///< MINOS diagnostic for log(R).
    MinosDiagnostic minos_amplitude;   ///< MINOS diagnostic for log(A_G).
    std::optional<double> q_min;        ///< Smallest valid objective found.
    /// Physical R and asymmetric MINOS errors when fully valid.
    std::optional<FitParameterEstimate> radius;
    /// Positive Gaussian amplitude A_G and asymmetric MINOS errors.
    std::optional<FitParameterEstimate> amplitude;
    /// Valid fitted bin densities; empty when the fit is not fully valid.
    std::vector<double> fitted_pdf;
};

/**
 * @brief Result of one independent three-parameter mixed fit.
 *
 * Thirty-six deterministic starts form the Cartesian product
 *
 * - R_core in {R_G, 0.5 R_HM, R_HM, 2 R_HM},
 * - R_tail in {0.5 R_tail,mom, R_tail,mom, 2 R_tail,mom}, and
 * - f_core in {0.25, 0.50, 0.75}.
 *
 * R_G always comes from the pure-Gaussian fit using the same estimator. Start
 * indices use core-major, then tail-major, then fraction-major ordering:
 * `index = (core_index * 3 + tail_index) * 3 + fraction_index`. Valid MIGRAD
 * endpoints are grouped into numerical basins. After the origin-specific
 * strict f_core admissibility filter, the basin reached by the largest number
 * of converged deterministic starts is selected. Equal-size basins are ranked
 * by their smallest q and then by lowest start index. The valid minimum with
 * the smallest q inside the selected basin is then selected for MINOS.
 * `consensus_size` records the selected basin multiplicity and is not an
 * independent acceptance veto. R_HM remains a deterministic seed only.
 */
struct MixedFitResult {
    /** Number of deterministic Cartesian-product mixed MIGRAD starts. */
    static constexpr std::size_t kCoreStartCount = 36U;

    bool fully_valid;                  ///< Selected fit and all MINOS are valid.
    FitFailureReason failure_reason;   ///< Primary invalidity cause, if any.
    FitEstimator estimator;            ///< Objective minimized by this fit only.
    /// Diagnostics for all 36 deterministic starts.
    std::array<MigradDiagnostic, kCoreStartCount> starts;
    /// Final physical endpoint of every start, retained for basin diagnostics.
    std::array<MixedStartEndpointDiagnostic, kCoreStartCount> start_endpoints;
    std::size_t starts_attempted;      ///< Number of starts actually run.
    std::size_t valid_starts;          ///< Starts with valid evaluable MIGRAD state.
    std::size_t consensus_size;        ///< Same-basin multiplicity of selected minimum.
    /// Zero-based lowest-q start inside the selected largest basin.
    std::optional<std::size_t> selected_core_start;
    MigradDiagnostic selected_migrad;  ///< Selected minimum MIGRAD diagnostic.
    MinosDiagnostic minos_core_radius; ///< MINOS diagnostic for log(R_core).
    MinosDiagnostic minos_tail_radius; ///< MINOS diagnostic for log(R_tail).
    MinosDiagnostic minos_core_fraction; ///< MINOS diagnostic for f_core.
    std::optional<double> q_min;        ///< Selected basin minimum objective value.
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
    /// Full region used by the free-amplitude pure Gaussian fit.
    /// Retained under the historical field name for output compatibility.
    std::optional<StatisticalRegion> gaussian_core_region;
    std::uint64_t selected_count;            ///< N over the full region.
    /// Final normalized distribution over the full statistical region.
    std::vector<NormalizedHistogramBin> normalized_bins;
    /// Full-range free-amplitude pure Gaussian using Poisson deviance.
    GaussianFitResult gaussian;
    /// Independent full-range free-amplitude Gaussian using Neyman chi-square.
    GaussianFitResult gaussian_neyman;
    /// Independent full-range free-amplitude Gaussian using Pearson chi-square.
    GaussianFitResult gaussian_pearson;
    /// Full-range mixed fit using the default binned Poisson deviance.
    MixedFitResult mixed;
    /// Full-range mixed fit using Neyman chi-square as an independent estimator.
    MixedFitResult mixed_neyman;
    /// Full-range mixed fit using Pearson chi-square as an independent estimator.
    MixedFitResult mixed_pearson;
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
 * @brief Return a stable ASCII token for one fit estimator.
 * @param estimator Statistical objective used by one independent fit.
 * @return Static token: `poisson`, `neyman`, or `pearson`.
 * @throws std::invalid_argument If @p estimator is not a valid enum value.
 */
[[nodiscard]] const char* fit_estimator_token(FitEstimator estimator);

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
