/**
 * @file hbt_config.h
 * @brief Scientific configuration data required by the HBT module.
 *
 * This file defines the domain-facing configuration data used by the HBT
 * module after configuration parsing and validation.
 *
 * File access, YAML syntax, configuration-key lookup, textual channel parsing,
 * and configuration validation are outside these data types' responsibility.
 */

#ifndef HBT_CONFIG_HBT_CONFIG_H
#define HBT_CONFIG_HBT_CONFIG_H

#include "hbt/selection/hbt_selection.h"

#include <cstddef>
#include <vector>

namespace hbt {

    /**
     * @brief Longitudinal kinematic variable used for particle acceptance.
     *
     * Rapidity selects the particle rapidity y.
     * Pseudorapidity selects the particle pseudorapidity eta.
     *
     * The configured longitudinal acceptance is applied to the absolute value
     * of the selected variable.
     */
    enum class LongitudinalVariable {
        Rapidity,       ///< Particle rapidity y.
        Pseudorapidity  ///< Particle pseudorapidity eta.
    };


    /**
     * @brief Requested statistical estimator set for HBT fits.
     *
     * A single estimator executes only its matching pure-Gaussian and mixed
     * fits. All preserves the historical behavior and executes Poisson,
     * Neyman, and Pearson independently.
     */
    enum class FitEstimatorMode {
        Poisson, ///< Execute only Poisson fits.
        Neyman,  ///< Execute only Neyman fits.
        Pearson, ///< Execute only Pearson fits.
        All      ///< Execute Poisson, Neyman, and Pearson fits.
    };

    /**
     * @brief Kinematic acceptance cuts shared by one particle-species group.
     *
     * longitudinal_abs_max is the positive upper bound applied to the absolute
     * value of the configured longitudinal variable.
     *
     * pt_min_gev and pt_max_gev are the lower and upper transverse-momentum
     * bounds in GeV.
     *
     * The eventual particle-acceptance operation uses open boundaries:
     *
     *     |longitudinal variable| < longitudinal_abs_max
     *     pt_min_gev < pT < pt_max_gev
     *
     * This data type stores already validated values. Configuration validation
     * is the responsibility of the HBT configuration loader.
     */
    struct ParticleAcceptanceCuts {
        /// Positive upper bound on the absolute longitudinal variable.
        double longitudinal_abs_max;
        /// Open lower transverse-momentum bound in GeV.
        double pt_min_gev;
        /// Open upper transverse-momentum bound in GeV.
        double pt_max_gev;
    };

    /**
     * @brief Particle-level kinematic acceptance configuration for HBT.
     *
     * The longitudinal variable is common to all HBT particle groups.
     *
     * Each physical species belongs to exactly one configured acceptance group:
     *
     * pions:
     *     PiPlus, PiMinus, PiZero
     *
     * kaons:
     *     KPlus, KMinus, KZero, KZeroBar
     *
     * nucleons:
     *     Proton, ProtonBar, Neutron, NeutronBar
     *
     * sigmas:
     *     SigmaPlus, SigmaBarMinus, SigmaZero
     *
     * lambdas:
     *     Lambda, LambdaBar
     *
     * These groups define particle-level acceptance only. They do not select
     * primitive channels, construct pairs, or define pair-level KT or mT
     * slices.
     */
    struct ParticleAcceptanceConfig {
        /// Longitudinal variable shared by all configured species groups.
        LongitudinalVariable longitudinal_variable;
        /// Validated pion acceptance cuts.
        ParticleAcceptanceCuts pions;
        /// Validated kaon acceptance cuts.
        ParticleAcceptanceCuts kaons;
        /// Validated nucleon acceptance cuts.
        ParticleAcceptanceCuts nucleons;
        /// Validated sigma acceptance cuts.
        ParticleAcceptanceCuts sigmas;
        /// Validated lambda acceptance cuts.
        ParticleAcceptanceCuts lambdas;
    };

    /**
     * @brief Pair-level slicing configuration for one kinematic variable.
     *
     * enabled determines whether this variable participates in pair routing.
     *
     * When bin_edges_gev is present, it contains at least two finite,
     * non-negative, strictly increasing boundaries in GeV. Consecutive entries
     * define half-open slices [edge_i, edge_(i+1)).
     *
     * When enabled is true, bin_edges_gev is required. When enabled is false,
     * bin_edges_gev may be empty or may retain a validated production binning
     * that does not participate in pair routing. No implicit or dummy slicing
     * interval is represented.
     *
     * This data type stores already validated configuration. Pair kinematics,
     * slice lookup, pair construction, and result accumulation are outside its
     * responsibility.
     */
    struct PairSlicingAxisConfig {
        /// Whether this axis participates in pair routing.
        bool enabled;
        /// Validated half-open slice boundaries in GeV.
        std::vector<double> bin_edges_gev;
    };

    /**
     * @brief Independent pair-level slicing configuration for kT and mT.
     *
     * kt configures optional pair transverse-momentum slicing.
     * mt configures optional pair transverse-mass slicing.
     *
     * Either axis may be enabled independently, both may be disabled, or both
     * may be enabled simultaneously. This structure does not define the
     * physical formulas used to calculate kT or mT.
     */
    struct PairSlicingConfig {
        /// Optional pair-transverse-momentum slicing axis.
        PairSlicingAxisConfig kt;
        /// Optional pair-transverse-mass slicing axis.
        PairSlicingAxisConfig mt;
    };

    /**
     * @brief Validated uniform binning for one histogram family.
     *
     * nbins is the exact number of contiguous half-open bins spanning
     * [minimum, maximum). minimum and maximum use the physical units of the
     * owning histogram family. inverse_bin_width is resolved once during
     * configuration loading and is reused by the accumulation hot path.
     *
     * All values are explicit scientific configuration or validated derived
     * state. No default, clamp, epsilon, or fallback is represented here.
     */
    struct HistogramBinningConfig {
        /// Number of uniform contiguous bins.
        std::size_t nbins;
        /// Inclusive lower histogram boundary.
        double minimum;
        /// Exclusive upper histogram boundary.
        double maximum;
        /// Precomputed reciprocal bin width.
        double inverse_bin_width;
    };

    /**
     * @brief Explicit raw-histogram binning for the three HBT families.
     *
     * osl is shared by |r_out_lcms|, |r_out_prf|, |r_side|, and |r_long|.
     * radial is shared by r_radial_lcms and r_radial_prf. delta_t is shared
     * by delta_t_lab, delta_t_lcms, and delta_t_prf.
     *
     * The three families are independent and may use different bin counts,
     * ranges, and widths. OSL and radial start exactly at zero. Delta-t spans
     * zero with a signed range.
     */
    struct HBTHistogramConfig {
        /// Absolute OSL-component histogram binning in fm.
        HistogramBinningConfig osl;
        /// Radial-source histogram binning in fm.
        HistogramBinningConfig radial;
        /// Relative-time histogram binning in fm/c.
        HistogramBinningConfig delta_t;
    };

    /**
     * @brief Requested nested HBT origin selection mode.
     *
     * Primordial requests the most restrictive origin selection.
     *
     * PrimordialRescattering requests the inclusive selection containing
     * primordial and rescattering contributions.
     *
     * PrimordialRescatteringDecay requests the fully inclusive configured
     * origin selection.
     *
     * All requests all three nested origin selections simultaneously in the
     * same analysis run.
     *
     * All is not a fourth physical origin slice. A physical pair is calculated
     * once and may later be routed to more than one compatible nested origin
     * slice.
     */
    enum class OriginMode {
        /// Primordial pairs only.
        Primordial,
        /// Primordial plus rescattering pairs.
        PrimordialRescattering,
        /// Primordial, rescattering, and decay pairs.
        PrimordialRescatteringDecay,
        /// Request all three nested origin selections simultaneously.
        All
    };

    /**
     * @brief Scientific configuration required to initialize the HBT module.
     *
     * This type contains only HBT-specific scientific configuration. Global run
     * controls, module activation, input/output paths, resource settings, and
     * configuration-file locations belong to the global run configuration.
     *
     * selection is the resolved representation of the configured
     * hbt_enabled_channels expression.
     *
     * particle_acceptance stores the validated particle-level kinematic
     * acceptance configuration.
     *
     * pair_slicing stores the validated optional kT and mT slicing axes.
     *
     * histogram_config stores the explicit raw-histogram binning.
     *
     * origin_mode stores the requested nested HBT origin-selection mode.
     *
     * fit_estimator_mode stores which independent statistical estimator fits
     * are executed.
     *
     * All six members are required scientific state. No implicit scientific
     * defaults are provided.
     */
    struct HBTConfig {
        /// Resolved configured analysis products and primitive channels.
        HBTSelection selection;
        /// Validated particle-level kinematic acceptance configuration.
        ParticleAcceptanceConfig particle_acceptance;
        /// Validated optional kT and mT pair-slicing configuration.
        PairSlicingConfig pair_slicing;
        /// Validated explicit raw-histogram binning.
        HBTHistogramConfig histogram_config;
        /// Requested nested origin-selection mode.
        OriginMode origin_mode;
        /// Requested statistical estimator set for Gaussian and mixed fits.
        FitEstimatorMode fit_estimator_mode;
    };

}  // namespace hbt

#endif  // HBT_CONFIG_HBT_CONFIG_H
