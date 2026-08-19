#!/usr/bin/env python3
"""Plot SAAM production distributions without recalculating results.

Run:
    python scripts/plot_results.py /ruta/al/directorio_de_resultados
"""

from __future__ import annotations

import argparse
import csv
import sys
from pathlib import Path
from typing import Dict, List, Optional, Sequence, Tuple

try:
    import matplotlib

    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    from matplotlib.lines import Line2D
    from matplotlib.patches import Patch
except ImportError as exc:
    raise SystemExit(
        "plot_results.py requires matplotlib to be installed"
    ) from exc


HISTOGRAM_COLOR = "mediumseagreen"
BAND_COLOR = "lightgreen"
GAUSSIAN_COLOR = "royalblue"
MIXED_COLOR = "crimson"
HISTOGRAM_ALPHA = 0.45
BAND_ALPHA = 0.30
HISTOGRAM_LINEWIDTH = 1.3
FIT_LINEWIDTH = 2.2
EMPTY_BIN_MARGIN = 2

SPECIES_MATH = {
    "pi_plus": r"\pi^{+}",
    "pi_minus": r"\pi^{-}",
    "p": r"p",
    "p_bar": r"\bar{p}",
    "k_plus": r"K^{+}",
    "k_minus": r"K^{-}",
}

X_LABELS = {
    "r_out": r"$|r_{\mathrm{out}}|\;[\mathrm{fm}]$",
    "r_side": r"$|r_{\mathrm{side}}|\;[\mathrm{fm}]$",
    "r_long": r"$|r_{\mathrm{long}}|\;[\mathrm{fm}]$",
    "delta_t": r"$\Delta t\;[\mathrm{fm}/c]$",
}

RADIAL_X_LABELS = {
    "LCMS": (
        r"$r=\sqrt{r_{\mathrm{out}}^{2}+r_{\mathrm{side}}^{2}"
        r"+r_{\mathrm{long}}^{2}}\;[\mathrm{fm}]$"
    ),
    "PRF": (
        r"$r^{*}=\sqrt{(r_{\mathrm{out}}^{*})^{2}"
        r"+(r_{\mathrm{side}}^{*})^{2}"
        r"+(r_{\mathrm{long}}^{*})^{2}}\;[\mathrm{fm}]$"
    ),
}

FRAMES = {"LAB", "LCMS", "PRF"}
KINDS = {"osl", "radial", "dt"}


class PlotInputError(RuntimeError):
    """Raised when the supplied production-output tree is incomplete."""


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description=(
            "Plot SAAM distributions using only values already exported "
            "by the analysis program."
        )
    )
    parser.add_argument(
        "results_root",
        type=Path,
        help="path to the production-results root",
    )
    parser.add_argument(
        "--plots-dir",
        type=Path,
        default=None,
        help=(
            "optional figure destination; default: a sibling directory "
            "named <results_root>_plots"
        ),
    )
    return parser.parse_args()


def read_csv_rows(path: Path) -> List[Dict[str, str]]:
    try:
        with path.open("r", encoding="utf-8", newline="") as stream:
            return list(csv.DictReader(stream))
    except OSError as exc:
        raise PlotInputError(f"cannot read {path}: {exc}") from exc


def species_math(token: str) -> str:
    if token in SPECIES_MATH:
        return SPECIES_MATH[token]
    escaped = token.replace("_", r"\_")
    return rf"\mathrm{{{escaped}}}"


def load_product_labels(results_root: Path) -> Dict[int, str]:
    catalog_path = results_root / "product_catalog.csv"
    if not catalog_path.is_file():
        raise PlotInputError(
            f"missing product catalog: {catalog_path}"
        )

    rows = read_csv_rows(catalog_path)
    grouped: Dict[int, List[Tuple[int, str]]] = {}
    for row in rows:
        try:
            product_index = int(row["product_index"])
            channel_index = int(row["channel_index"])
            species_a = row["species_a"]
            species_b = row["species_b"]
        except (KeyError, ValueError) as exc:
            raise PlotInputError(
                f"invalid product catalog row in {catalog_path}"
            ) from exc

        channel = species_math(species_a) + species_math(species_b)
        grouped.setdefault(product_index, []).append(
            (channel_index, channel)
        )

    labels: Dict[int, str] = {}
    for product_index, channels in grouped.items():
        channels.sort(key=lambda item: item[0])
        expression = r" \oplus ".join(
            channel for _, channel in channels
        )
        labels[product_index] = f"${expression}$"
    return labels


def parse_product_index(path: Path) -> int:
    for part in path.parts:
        if part.startswith("product_"):
            suffix = part[len("product_") :]
            try:
                return int(suffix)
            except ValueError:
                continue
    raise PlotInputError(f"cannot determine product index from {path}")


def find_path_token(path: Path, values: set[str]) -> Optional[str]:
    for part in path.parts:
        if part in values:
            return part
    return None


def read_float_column(
    rows: Sequence[Dict[str, str]],
    column: str,
) -> Optional[List[float]]:
    if not rows or column not in rows[0]:
        return None

    values: List[float] = []
    for row in rows:
        text = row.get(column, "")
        if text == "":
            return None
        try:
            values.append(float(text))
        except ValueError as exc:
            raise PlotInputError(
                f"invalid value in column {column!r}"
            ) from exc
    return values


def require_float_column(
    rows: Sequence[Dict[str, str]],
    column: str,
    path: Path,
) -> List[float]:
    values = read_float_column(rows, column)
    if values is None:
        raise PlotInputError(
            f"missing required column {column!r} in {path}"
        )
    return values


def visible_bin_range(pdf: Sequence[float]) -> Tuple[int, int]:
    occupied = [index for index, value in enumerate(pdf) if value != 0.0]
    if not occupied:
        return 0, len(pdf) - 1

    first = max(0, occupied[0] - EMPTY_BIN_MARGIN)
    last = min(len(pdf) - 1, occupied[-1] + EMPTY_BIN_MARGIN)
    return first, last


def x_axis_label(kind: str, frame: str, observable: str) -> str:
    if kind == "radial":
        try:
            return RADIAL_X_LABELS[frame]
        except KeyError as exc:
            raise PlotInputError(
                f"radial distributions are unsupported in frame {frame!r}"
            ) from exc

    try:
        return X_LABELS[observable]
    except KeyError as exc:
        raise PlotInputError(
            f"unsupported observable {observable!r}"
        ) from exc


def y_axis_label(kind: str, observable: str) -> str:
    if kind == "osl":
        component = {
            "r_out": r"r_{\mathrm{out}}",
            "r_side": r"r_{\mathrm{side}}",
            "r_long": r"r_{\mathrm{long}}",
        }
        try:
            argument = component[observable]
        except KeyError as exc:
            raise PlotInputError(
                f"unsupported OSL observable {observable!r}"
            ) from exc
        return rf"$S_{{\mathrm{{proj}}}}(|{argument}|)\;[\mathrm{{fm}}^{{-1}}]$"

    if kind == "radial":
        return r"$S(r)\;[\mathrm{fm}^{-1}]$"

    if kind == "dt":
        return r"$P(\Delta t)\;[(\mathrm{fm}/c)^{-1}]$"

    raise PlotInputError(f"unsupported distribution type {kind!r}")


def add_channel_and_fit_legend(
    ax: plt.Axes,
    channel_label: str,
    gaussian: Optional[Sequence[float]],
    mixed: Optional[Sequence[float]],
) -> None:
    handles = [
        Patch(
            facecolor=HISTOGRAM_COLOR,
            edgecolor=HISTOGRAM_COLOR,
            alpha=HISTOGRAM_ALPHA,
            label=channel_label,
        )
    ]
    if gaussian is not None:
        handles.append(
            Line2D(
                [0],
                [0],
                color=GAUSSIAN_COLOR,
                linewidth=FIT_LINEWIDTH,
                label="Gaussian fit",
            )
        )
    if mixed is not None:
        handles.append(
            Line2D(
                [0],
                [0],
                color=MIXED_COLOR,
                linewidth=FIT_LINEWIDTH,
                label="Gauss + exponential fit",
            )
        )
    ax.legend(handles=handles, loc="best", frameon=False)


def plot_distribution(
    distribution_path: Path,
    results_root: Path,
    plots_root: Path,
    product_labels: Dict[int, str],
) -> Path:
    relative_path = distribution_path.relative_to(results_root)
    product_index = parse_product_index(relative_path)
    try:
        channel_label = product_labels[product_index]
    except KeyError as exc:
        raise PlotInputError(
            f"product_{product_index} is absent from product_catalog.csv"
        ) from exc

    frame = find_path_token(relative_path, FRAMES)
    kind = find_path_token(relative_path, KINDS)
    observable = distribution_path.parent.name
    if frame is None or kind is None:
        raise PlotInputError(
            f"cannot determine frame/type from {relative_path}"
        )
    supported_observables = set(X_LABELS) | {"r_radial"}
    if observable not in supported_observables:
        raise PlotInputError(
            f"unsupported observable {observable!r} in {relative_path}"
        )

    rows = read_csv_rows(distribution_path)
    if not rows:
        raise PlotInputError(f"empty distribution: {distribution_path}")

    lower = require_float_column(rows, "lower_edge", distribution_path)
    upper = require_float_column(rows, "upper_edge", distribution_path)
    centers = require_float_column(rows, "center", distribution_path)
    pdf = require_float_column(rows, "pdf", distribution_path)
    d_pdf = require_float_column(rows, "d_pdf", distribution_path)
    gaussian = read_float_column(rows, "gaussian_fit_pdf")
    mixed = read_float_column(rows, "mixed_fit_pdf")

    edges = [lower[0], *upper]
    band_lower = [value - error for value, error in zip(pdf, d_pdf)]
    band_upper = [value + error for value, error in zip(pdf, d_pdf)]
    if kind == "osl":
        band_lower = [max(value, 1.0e-7) for value in band_lower]
        band_upper = [max(value, 1.0e-7) for value in band_upper]
    else:
        band_lower = [max(value, 0.0) for value in band_lower]

    fig, ax = plt.subplots(figsize=(6.0, 4.5))
    ax.stairs(
        pdf,
        edges,
        baseline=0.0,
        fill=True,
        facecolor=HISTOGRAM_COLOR,
        edgecolor="none",
        alpha=HISTOGRAM_ALPHA,
        zorder=1,
    )
    ax.fill_between(
        centers,
        band_lower,
        band_upper,
        step="mid",
        color=BAND_COLOR,
        alpha=BAND_ALPHA,
        linewidth=0.0,
        zorder=2,
    )
    ax.stairs(
        pdf,
        edges,
        color=HISTOGRAM_COLOR,
        alpha=1.0,
        linewidth=HISTOGRAM_LINEWIDTH,
        zorder=3,
    )
    if gaussian is not None:
        ax.plot(
            centers,
            gaussian,
            color=GAUSSIAN_COLOR,
            linewidth=FIT_LINEWIDTH,
            zorder=4,
        )
    if mixed is not None:
        ax.plot(
            centers,
            mixed,
            color=MIXED_COLOR,
            linewidth=FIT_LINEWIDTH,
            zorder=5,
        )

    if kind == "osl":
        ax.set_yscale("log")
        ax.set_ylim(bottom=1.0e-7)

    first, last = visible_bin_range(pdf)
    ax.set_xlim(lower[first], upper[last])
    ax.margins(y=0.08)
    ax.set_xlabel(x_axis_label(kind, frame, observable))
    ax.set_ylabel(y_axis_label(kind, observable))
    ax.grid(True, alpha=0.25)

    add_channel_and_fit_legend(
        ax,
        channel_label,
        gaussian,
        mixed,
    )

    ax.text(
        0.97,
        0.97,
        frame,
        transform=ax.transAxes,
        horizontalalignment="right",
        verticalalignment="top",
    )
    fig.tight_layout()

    output_path = plots_root / relative_path.parent / "distribution.png"
    output_path.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(output_path, dpi=200, bbox_inches="tight")
    plt.close(fig)
    return output_path


def main() -> int:
    args = parse_args()
    results_root = args.results_root.expanduser().resolve()
    if not results_root.is_dir():
        print(
            f"error: results root is not a directory: {results_root}",
            file=sys.stderr,
        )
        return 2

    if args.plots_dir is None:
        plots_root = results_root.with_name(
            f"{results_root.name}_plots"
        )
    else:
        plots_root = args.plots_dir.expanduser().resolve()

    try:
        product_labels = load_product_labels(results_root)
        distributions = sorted(results_root.rglob("distribution.csv"))
        if not distributions:
            raise PlotInputError(
                f"no distribution.csv files found under {results_root}"
            )

        for distribution_path in distributions:
            plot_distribution(
                distribution_path,
                results_root,
                plots_root,
                product_labels,
            )
    except PlotInputError as exc:
        print(f"error: {exc}", file=sys.stderr)
        return 2

    print(f"generated {len(distributions)} figures in {plots_root}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
