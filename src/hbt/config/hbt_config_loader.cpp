/**
 * @file hbt_config_loader.cpp
 * @brief Loading of the scientific HBT module configuration from a YAML file.
 */

#include "hbt/config/hbt_config_loader.h"

#include "hbt/config/hbt_enabled_channels_parser.h"

#include <yaml-cpp/yaml.h>

#include <algorithm>
#include <charconv>
#include <cmath>
#include <initializer_list>
#include <map>
#include <stdexcept>
#include <string>
#include <string_view>
#include <system_error>
#include <vector>

namespace hbt {

namespace {

    /**
     * @brief Load one YAML document from the HBT configuration file.
     *
     * YAML-library exceptions produced while opening or parsing the file are
     * translated into std::runtime_error so that YAML-specific exception types
     * do not escape through the public configuration interface.
     *
     * @param path Path to the HBT module configuration file.
     *
     * @return Parsed YAML root node.
     *
     * @throws std::runtime_error if the file cannot be read or the YAML
     *         document cannot be parsed.
     */
    YAML::Node load_yaml_document(
        const std::filesystem::path& path
    ) {
        try {
            return YAML::LoadFile(path.string());
        } catch (const YAML::Exception& exception) {
            throw std::runtime_error(
                "failed to load HBT configuration file '" +
                path.string() +
                "': " +
                exception.what()
            );
        }
    }

    /**
     * @brief Extract and validate one exact YAML mapping.
     *
     * The mapping must contain every expected key exactly once and no
     * additional keys.
     *
     * @param node YAML node to validate.
     * @param context Human-readable mapping description for error messages.
     * @param expected_keys Exact keys permitted in the mapping.
     *
     * @return Validated mapping from key text to YAML value nodes.
     *
     * @throws std::runtime_error if the node is not a mapping, a key is not a
     *         scalar, a key is unknown or duplicated, or a required key is
     *         missing.
     */
    std::map<std::string, YAML::Node> extract_exact_mapping(
        const YAML::Node& node,
        std::string_view context,
        std::initializer_list<std::string_view> expected_keys
    ) {
        if (!node.IsMap()) {
            throw std::runtime_error(
                std::string(context) +
                " must be a YAML mapping"
            );
        }

        std::map<std::string, YAML::Node> entries;

        for (const auto& entry : node) {
            if (!entry.first.IsScalar()) {
                throw std::runtime_error(
                    std::string(context) +
                    " keys must be YAML scalars"
                );
            }

            const std::string key = entry.first.Scalar();

            const auto expected_iterator = std::find(
                expected_keys.begin(),
                expected_keys.end(),
                std::string_view{key}
            );

            if (expected_iterator == expected_keys.end()) {
                throw std::runtime_error(
                    std::string(context) +
                    " contains unknown key '" +
                    key +
                    "'"
                );
            }

            if (entries.find(key) != entries.end()) {
                throw std::runtime_error(
                    std::string(context) +
                    " contains duplicate key '" +
                    key +
                    "'"
                );
            }

            entries.emplace(key, entry.second);
        }

        for (const std::string_view expected_key : expected_keys) {
            if (
                entries.find(std::string(expected_key)) ==
                entries.end()
            ) {
                throw std::runtime_error(
                    std::string(context) +
                    " is missing required key '" +
                    std::string(expected_key) +
                    "'"
                );
            }
        }

        return entries;
    }

    /**
     * @brief Extract one required YAML scalar as text.
     *
     * @param node YAML node expected to contain a scalar.
     * @param path Configuration path used in error messages.
     *
     * @return Scalar text.
     *
     * @throws std::runtime_error if @p node is not a YAML scalar.
     */
    std::string extract_scalar(
        const YAML::Node& node,
        std::string_view path
    ) {
        if (!node.IsScalar()) {
            throw std::runtime_error(
                "HBT configuration entry '" +
                std::string(path) +
                "' must be a YAML scalar"
            );
        }

        return node.Scalar();
    }

    /**
     * @brief Parse one finite floating-point scientific configuration value.
     *
     * @param node YAML scalar containing the numeric value.
     * @param path Configuration path used in error messages.
     *
     * @return Parsed finite double value.
     *
     * @throws std::runtime_error if @p node is not a YAML scalar.
     * @throws std::invalid_argument if the scalar is not a finite number.
     */
    double parse_finite_double(
        const YAML::Node& node,
        std::string_view path
    ) {
        static_cast<void>(
            extract_scalar(node, path)
        );

        double value = 0.0;

        try {
            value = node.as<double>();
        } catch (const YAML::Exception&) {
            throw std::invalid_argument(
                "HBT configuration value '" +
                std::string(path) +
                "' must be a finite number"
            );
        }

        if (!std::isfinite(value)) {
            throw std::invalid_argument(
                "HBT configuration value '" +
                std::string(path) +
                "' must be finite"
            );
        }

        return value;
    }

    /**
     * @brief Parse one strictly positive histogram bin count.
     *
     * Only unsigned decimal integer text is accepted. Signs, decimal points,
     * scientific notation, and values that do not fit std::size_t are
     * rejected rather than interpreted implicitly by the YAML library.
     *
     * @param node YAML scalar containing the bin count.
     * @param path Configuration path used in error messages.
     * @return Parsed positive bin count.
     * @throws std::runtime_error If @p node is not a YAML scalar.
     * @throws std::invalid_argument If the scalar is not a positive integer.
     */
    std::size_t parse_positive_bin_count(
        const YAML::Node& node,
        std::string_view path
    ) {
        const std::string text = extract_scalar(node, path);
        if (text.empty()) {
            throw std::invalid_argument(
                "HBT configuration value '" +
                std::string(path) +
                "' must be a positive integer"
            );
        }

        std::size_t value = 0U;
        const char* const begin = text.data();
        const char* const end = begin + text.size();
        const auto result = std::from_chars(begin, end, value, 10);

        if (result.ec != std::errc{} || result.ptr != end || value == 0U) {
            throw std::invalid_argument(
                "HBT configuration value '" +
                std::string(path) +
                "' must be a positive integer"
            );
        }

        return value;
    }

    /**
     * @brief Parse one YAML boolean scalar.
     *
     * @param node YAML node containing the boolean value.
     * @param path Configuration path used in error messages.
     *
     * @return Parsed boolean value.
     *
     * @throws std::runtime_error if @p node is not a YAML boolean scalar.
     */
    bool parse_boolean(
        const YAML::Node& node,
        std::string_view path
    ) {
        if (!node.IsScalar()) {
            throw std::runtime_error(
                "HBT configuration entry '" +
                std::string(path) +
                "' must be a YAML boolean scalar"
            );
        }

        try {
            return node.as<bool>();
        } catch (const YAML::Exception&) {
            throw std::runtime_error(
                "HBT configuration entry '" +
                std::string(path) +
                "' must be a YAML boolean scalar"
            );
        }
    }

    /**
     * @brief Parse validated consecutive bin edges for one pair-slicing axis.
     *
     * @param node YAML sequence containing bin boundaries in GeV.
     * @param path Configuration path used in error messages.
     *
     * @return Finite, non-negative, strictly increasing bin edges.
     *
     * @throws std::runtime_error if @p node is not a YAML sequence or one
     *         sequence element is not a scalar.
     * @throws std::invalid_argument if fewer than two edges are provided or
     *         an edge violates the scientific slicing contract.
     */
    std::vector<double> parse_pair_slicing_bin_edges(
        const YAML::Node& node,
        std::string_view path
    ) {
        if (!node.IsSequence()) {
            throw std::runtime_error(
                "HBT configuration entry '" +
                std::string(path) +
                "' must be a YAML sequence"
            );
        }

        if (node.size() < 2U) {
            throw std::invalid_argument(
                "HBT configuration entry '" +
                std::string(path) +
                "' must contain at least two bin edges"
            );
        }

        std::vector<double> edges;
        edges.reserve(node.size());

        for (std::size_t index = 0U; index < node.size(); ++index) {
            const std::string edge_path =
                std::string(path) +
                "[" +
                std::to_string(index) +
                "]";

            const double edge =
                parse_finite_double(node[index], edge_path);

            if (edge < 0.0) {
                throw std::invalid_argument(
                    "HBT configuration value '" +
                    edge_path +
                    "' must be >= 0"
                );
            }

            if (!edges.empty() && !(edge > edges.back())) {
                throw std::invalid_argument(
                    "HBT configuration entry '" +
                    std::string(path) +
                    "' must be strictly increasing"
                );
            }

            edges.push_back(edge);
        }

        return edges;
    }

    /**
     * @brief Parse one optional pair-slicing axis.
     *
     * The axis mapping always requires enabled. When enabled is true,
     * bin_edges_gev is also required. When enabled is false, bin_edges_gev is
     * optional; if present, it is validated and retained for production
     * configurations that toggle the axis without rewriting its binning.
     *
     * @param node YAML mapping containing one pair-slicing axis.
     * @param axis_name Canonical axis name, either kt or mt.
     *
     * @return Fully resolved pair-slicing axis configuration.
     *
     * @throws std::runtime_error if the YAML structure is invalid.
     * @throws std::invalid_argument if configured bin edges are invalid.
     */
    PairSlicingAxisConfig parse_pair_slicing_axis(
        const YAML::Node& node,
        std::string_view axis_name
    ) {
        const std::string context =
            "HBT pair-slicing axis '" +
            std::string(axis_name) +
            "'";

        if (!node.IsMap()) {
            throw std::runtime_error(
                context +
                " must be a YAML mapping"
            );
        }

        bool has_enabled = false;
        bool enabled = false;
        bool has_bin_edges = false;
        YAML::Node bin_edges;

        for (const auto& entry : node) {
            if (!entry.first.IsScalar()) {
                throw std::runtime_error(
                    context +
                    " keys must be YAML scalars"
                );
            }

            const std::string key = entry.first.Scalar();

            if (key == "enabled") {
                if (has_enabled) {
                    throw std::runtime_error(
                        context +
                        " contains duplicate key 'enabled'"
                    );
                }

                enabled = parse_boolean(
                    entry.second,
                    "hbt_pair_slicing." +
                    std::string(axis_name) +
                    ".enabled"
                );
                has_enabled = true;
                continue;
            }

            if (key == "bin_edges_gev") {
                if (has_bin_edges) {
                    throw std::runtime_error(
                        context +
                        " contains duplicate key 'bin_edges_gev'"
                    );
                }

                bin_edges = entry.second;
                has_bin_edges = true;
                continue;
            }

            throw std::runtime_error(
                context +
                " contains unknown key '" +
                key +
                "'"
            );
        }

        if (!has_enabled) {
            throw std::runtime_error(
                context +
                " is missing required key 'enabled'"
            );
        }

        if (!enabled) {
            if (!has_bin_edges) {
                return PairSlicingAxisConfig{false, {}};
            }

            return PairSlicingAxisConfig{
                false,
                parse_pair_slicing_bin_edges(
                    bin_edges,
                    "hbt_pair_slicing." +
                    std::string(axis_name) +
                    ".bin_edges_gev"
                )
            };
        }

        if (!has_bin_edges) {
            throw std::runtime_error(
                context +
                " is missing required key 'bin_edges_gev' when enabled"
            );
        }

        return PairSlicingAxisConfig{
            true,
            parse_pair_slicing_bin_edges(
                bin_edges,
                "hbt_pair_slicing." +
                std::string(axis_name) +
                ".bin_edges_gev"
            )
        };
    }

    /**
     * @brief Parse independent kT and mT pair-slicing configuration.
     *
     * @param node YAML mapping containing exactly kt and mt axes.
     *
     * @return Fully resolved pair-slicing configuration.
     *
     * @throws std::runtime_error if the YAML structure is invalid.
     * @throws std::invalid_argument if configured bin edges are invalid.
     */
    PairSlicingConfig parse_pair_slicing(
        const YAML::Node& node
    ) {
        const std::map<std::string, YAML::Node> entries =
            extract_exact_mapping(
                node,
                "HBT pair slicing",
                {
                    "kt",
                    "mt"
                }
            );

        return PairSlicingConfig{
            parse_pair_slicing_axis(entries.at("kt"), "kt"),
            parse_pair_slicing_axis(entries.at("mt"), "mt")
        };
    }

    /**
     * @brief Parse and validate one explicit histogram-family binning.
     * @param node YAML mapping containing nbins and physical boundaries.
     * @param family_name Stable family name used in diagnostics.
     * @param minimum_key Family-specific minimum-boundary key.
     * @param maximum_key Family-specific maximum-boundary key.
     * @param require_zero_minimum Whether the minimum must equal zero exactly.
     * @param require_signed_zero_span Whether the range must span zero.
     * @return Validated binning with reciprocal width resolved once.
     * @throws std::runtime_error If the YAML structure is invalid.
     * @throws std::invalid_argument If any scientific binning value is invalid.
     */
    HistogramBinningConfig parse_histogram_binning(
        const YAML::Node& node,
        std::string_view family_name,
        std::string_view minimum_key,
        std::string_view maximum_key,
        bool require_zero_minimum,
        bool require_signed_zero_span
    ) {
        const std::string context =
            "HBT histogram family '" +
            std::string(family_name) +
            "'";
        const std::map<std::string, YAML::Node> entries =
            extract_exact_mapping(
                node,
                context,
                {"nbins", minimum_key, maximum_key}
            );
        const std::string prefix =
            "hbt_histograms." +
            std::string(family_name) +
            ".";
        const std::size_t nbins = parse_positive_bin_count(
            entries.at("nbins"),
            prefix + "nbins"
        );
        const double minimum = parse_finite_double(
            entries.at(std::string(minimum_key)),
            prefix + std::string(minimum_key)
        );
        const double maximum = parse_finite_double(
            entries.at(std::string(maximum_key)),
            prefix + std::string(maximum_key)
        );

        if (!(maximum > minimum)) {
            throw std::invalid_argument(
                "HBT histogram family '" +
                std::string(family_name) +
                "' requires maximum > minimum"
            );
        }
        if (require_zero_minimum && minimum != 0.0) {
            throw std::invalid_argument(
                "HBT histogram family '" +
                std::string(family_name) +
                "' requires minimum == 0"
            );
        }
        if (require_signed_zero_span &&
            !(minimum < 0.0 && maximum > 0.0)) {
            throw std::invalid_argument(
                "HBT histogram family '" +
                std::string(family_name) +
                "' requires minimum < 0 < maximum"
            );
        }

        const double span = maximum - minimum;
        const double bin_count = static_cast<double>(nbins);
        const double bin_width = span / bin_count;
        const double inverse_bin_width = bin_count / span;

        if (!std::isfinite(span) || !(span > 0.0) ||
            !std::isfinite(bin_width) || !(bin_width > 0.0) ||
            !std::isfinite(inverse_bin_width) ||
            !(inverse_bin_width > 0.0)) {
            throw std::invalid_argument(
                "HBT histogram family '" +
                std::string(family_name) +
                "' has a non-finite or non-positive bin width"
            );
        }

        return {nbins, minimum, maximum, inverse_bin_width};
    }

    /**
     * @brief Parse all explicit raw-histogram binning families.
     * @param node YAML mapping containing osl, radial, and delta_t families.
     * @return Fully validated raw-histogram configuration.
     * @throws std::runtime_error If the YAML structure is invalid.
     * @throws std::invalid_argument If any scientific binning value is invalid.
     */
    HBTHistogramConfig parse_histogram_config(
        const YAML::Node& node
    ) {
        const std::map<std::string, YAML::Node> entries =
            extract_exact_mapping(
                node,
                "HBT histogram configuration",
                {"osl", "radial", "delta_t"}
            );

        return {
            parse_histogram_binning(
                entries.at("osl"),
                "osl",
                "min_fm",
                "max_fm",
                true,
                false
            ),
            parse_histogram_binning(
                entries.at("radial"),
                "radial",
                "min_fm",
                "max_fm",
                true,
                false
            ),
            parse_histogram_binning(
                entries.at("delta_t"),
                "delta_t",
                "min_fm_c",
                "max_fm_c",
                false,
                true
            )
        };
    }

    /**
     * @brief Parse the longitudinal particle-acceptance variable.
     *
     * Accepted tokens are exactly rapidity and pseudorapidity.
     *
     * @param node YAML scalar containing the configured token.
     *
     * @return Resolved longitudinal-variable identifier.
     *
     * @throws std::runtime_error if @p node is not a YAML scalar.
     * @throws std::invalid_argument if the token is not canonical.
     */
    LongitudinalVariable parse_longitudinal_variable(
        const YAML::Node& node
    ) {
        const std::string value = extract_scalar(
            node,
            "hbt_particle_acceptance.longitudinal_variable"
        );

        if (value == "rapidity") {
            return LongitudinalVariable::Rapidity;
        }

        if (value == "pseudorapidity") {
            return LongitudinalVariable::Pseudorapidity;
        }

        throw std::invalid_argument(
            "invalid hbt_particle_acceptance.longitudinal_variable '" +
            value +
            "': expected rapidity or pseudorapidity"
        );
    }

    /**
     * @brief Parse and validate one particle-group acceptance configuration.
     *
     * The group mapping contains exactly longitudinal_abs_max, pt_min_gev,
     * and pt_max_gev.
     *
     * @param node YAML mapping containing the group cuts.
     * @param group_name Canonical HBT particle-group name.
     *
     * @return Validated particle-acceptance cuts.
     *
     * @throws std::runtime_error if the YAML structure is invalid.
     * @throws std::invalid_argument if a numeric value violates the scientific
     *         acceptance contract.
     */
    ParticleAcceptanceCuts parse_particle_acceptance_cuts(
        const YAML::Node& node,
        std::string_view group_name
    ) {
        const std::string context =
            "HBT particle-acceptance group '" +
            std::string(group_name) +
            "'";

        const std::map<std::string, YAML::Node> entries =
            extract_exact_mapping(
                node,
                context,
                {
                    "longitudinal_abs_max",
                    "pt_min_gev",
                    "pt_max_gev"
                }
            );

        const std::string prefix =
            "hbt_particle_acceptance.groups." +
            std::string(group_name) +
            ".";

        const double longitudinal_abs_max =
            parse_finite_double(
                entries.at("longitudinal_abs_max"),
                prefix + "longitudinal_abs_max"
            );

        const double pt_min_gev =
            parse_finite_double(
                entries.at("pt_min_gev"),
                prefix + "pt_min_gev"
            );

        const double pt_max_gev =
            parse_finite_double(
                entries.at("pt_max_gev"),
                prefix + "pt_max_gev"
            );

        if (!(longitudinal_abs_max > 0.0)) {
            throw std::invalid_argument(
                "HBT configuration value '" +
                prefix +
                "longitudinal_abs_max' must be > 0"
            );
        }

        if (!(pt_min_gev >= 0.0)) {
            throw std::invalid_argument(
                "HBT configuration value '" +
                prefix +
                "pt_min_gev' must be >= 0"
            );
        }

        if (!(pt_max_gev > pt_min_gev)) {
            throw std::invalid_argument(
                "HBT configuration value '" +
                prefix +
                "pt_max_gev' must be > pt_min_gev"
            );
        }

        return ParticleAcceptanceCuts{
            longitudinal_abs_max,
            pt_min_gev,
            pt_max_gev
        };
    }

    /**
     * @brief Parse the complete HBT particle-acceptance configuration.
     *
     * The particle-acceptance mapping contains exactly longitudinal_variable
     * and groups. The groups mapping contains exactly pions, kaons, nucleons,
     * sigmas, and lambdas.
     *
     * @param node YAML node containing hbt_particle_acceptance.
     *
     * @return Fully resolved particle-acceptance configuration.
     *
     * @throws std::runtime_error if the YAML structure is invalid.
     * @throws std::invalid_argument if a scientific value is invalid.
     */
    ParticleAcceptanceConfig parse_particle_acceptance(
        const YAML::Node& node
    ) {
        const std::map<std::string, YAML::Node> acceptance_entries =
            extract_exact_mapping(
                node,
                "HBT particle acceptance",
                {
                    "longitudinal_variable",
                    "groups"
                }
            );

        const std::map<std::string, YAML::Node> group_entries =
            extract_exact_mapping(
                acceptance_entries.at("groups"),
                "HBT particle-acceptance groups",
                {
                    "pions",
                    "kaons",
                    "nucleons",
                    "sigmas",
                    "lambdas"
                }
            );

        return ParticleAcceptanceConfig{
            parse_longitudinal_variable(
                acceptance_entries.at("longitudinal_variable")
            ),
            parse_particle_acceptance_cuts(
                group_entries.at("pions"),
                "pions"
            ),
            parse_particle_acceptance_cuts(
                group_entries.at("kaons"),
                "kaons"
            ),
            parse_particle_acceptance_cuts(
                group_entries.at("nucleons"),
                "nucleons"
            ),
            parse_particle_acceptance_cuts(
                group_entries.at("sigmas"),
                "sigmas"
            ),
            parse_particle_acceptance_cuts(
                group_entries.at("lambdas"),
                "lambdas"
            )
        };
    }

    /**
     * @brief Parse the requested nested HBT origin-selection mode.
     *
     * Accepted tokens are exactly primordial, primordial_rescattering,
     * primordial_rescattering_decay, and all.
     *
     * @param node YAML scalar containing hbt_origin_mode.
     *
     * @return Resolved origin-selection mode.
     *
     * @throws std::runtime_error if @p node is not a YAML scalar.
     * @throws std::invalid_argument if the token is not canonical.
     */
    OriginMode parse_origin_mode(
        const YAML::Node& node
    ) {
        const std::string value =
            extract_scalar(node, "hbt_origin_mode");

        if (value == "primordial") {
            return OriginMode::Primordial;
        }

        if (value == "primordial_rescattering") {
            return OriginMode::PrimordialRescattering;
        }

        if (value == "primordial_rescattering_decay") {
            return OriginMode::PrimordialRescatteringDecay;
        }

        if (value == "all") {
            return OriginMode::All;
        }

        throw std::invalid_argument(
            "invalid hbt_origin_mode '" +
            value +
            "': expected primordial, primordial_rescattering, "
            "primordial_rescattering_decay, or all"
        );
    }

}  // namespace

/*
 * @brief Load and resolve the scientific HBT configuration from YAML.
 *
 * @param path Path to the HBT module YAML configuration file.
 *
 * @return Fully resolved scientific HBT configuration.
 */
HBTConfig load_hbt_config(
    const std::filesystem::path& path
) {
    const YAML::Node root = load_yaml_document(path);

    const std::map<std::string, YAML::Node> root_entries =
        extract_exact_mapping(
            root,
            "HBT configuration root",
            {
                "hbt_enabled_channels",
                "hbt_particle_acceptance",
                "hbt_pair_slicing",
                "hbt_histograms",
                "hbt_origin_mode"
            }
        );

    const std::string enabled_channels =
        extract_scalar(
            root_entries.at("hbt_enabled_channels"),
            "hbt_enabled_channels"
        );

    return HBTConfig{
        parse_hbt_enabled_channels(enabled_channels),
        parse_particle_acceptance(
            root_entries.at("hbt_particle_acceptance")
        ),
        parse_pair_slicing(
            root_entries.at("hbt_pair_slicing")
        ),
        parse_histogram_config(
            root_entries.at("hbt_histograms")
        ),
        parse_origin_mode(
            root_entries.at("hbt_origin_mode")
        )
    };
}

}  // namespace hbt
