"""SVG plotting helpers for CABOMB generated figures."""

from __future__ import annotations

from pathlib import Path
from typing import Iterable

import matplotlib as mpl

mpl.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

from cabomb_model.solver import SimulationResult


def _polyline(points: Iterable[tuple[float, float]]) -> str:
    return " ".join(f"{x:.2f},{y:.2f}" for x, y in points)


def _save_deterministic_svg(figure: plt.Figure, path: str | Path) -> Path:
    """Save and close an SVG without time-varying metadata or random IDs."""

    output = Path(path)
    output.parent.mkdir(parents=True, exist_ok=True)
    with mpl.rc_context({"svg.hashsalt": "cabomb-v1.1"}):
        figure.savefig(
            output,
            format="svg",
            bbox_inches="tight",
            metadata={"Date": None, "Creator": "CABOMB V1.1"},
        )
    svg = output.read_text(encoding="utf-8")
    normalized = "\n".join(line.rstrip(" \t") for line in svg.splitlines()) + "\n"
    output.write_text(normalized, encoding="utf-8", newline="\n")
    plt.close(figure)
    return output


def write_cargo_release_vs_alpha(
    rows: list[dict[str, object]], path: str | Path
) -> Path:
    """Plot gross cumulative cargo release across the exploratory alpha sweep."""

    figure, axis = plt.subplots(figsize=(8.0, 5.2))
    for design in ("unmodified", "ocn", "ysa", "ocn-ysa"):
        selected = sorted(
            (row for row in rows if row["design"] == design),
            key=lambda row: float(row["alpha_ocn_cell_access"]),
        )
        axis.plot(
            [float(row["alpha_ocn_cell_access"]) for row in selected],
            [float(row["cumulative_cargo_release"]) for row in selected],
            marker="o",
            label=design,
        )
    axis.set_title("Cumulative cargo release across exploratory OCN accessibility")
    axis.set_xlabel("OCN-bound cellular accessibility, alpha (dimensionless)")
    axis.set_ylabel("Cumulative cargo release (normalized cargo)")
    axis.legend(title="Carrier design")
    axis.grid(alpha=0.25)
    figure.subplots_adjust(bottom=0.24)
    figure.text(
        0.01,
        0.01,
        "Alpha is the effective YSA-EphA2 access retained by OCN-bound carrier.\n"
        "Greater access may shift carrier toward internalization and release; "
        "benefit is not guaranteed. Exploratory normalized outputs; not clinical "
        "predictions.",
        fontsize=7.5,
    )
    return _save_deterministic_svg(figure, path)


def write_auc_release_tradeoff(
    rows: list[dict[str, object]], path: str | Path
) -> Path:
    """Plot exposure AUC against gross cumulative cargo release."""

    figure, axis = plt.subplots(figsize=(8.0, 5.2))
    for design in ("unmodified", "ocn", "ysa", "ocn-ysa"):
        selected = sorted(
            (row for row in rows if row["design"] == design),
            key=lambda row: float(row["alpha_ocn_cell_access"]),
        )
        x = [float(row["tumor_auc"]) for row in selected]
        y = [float(row["cumulative_cargo_release"]) for row in selected]
        axis.plot(x, y, marker="o", label=design)
        if design == "ocn-ysa":
            for row, x_value, y_value in zip(selected, x, y):
                axis.annotate(
                    f"alpha={float(row['alpha_ocn_cell_access']):g}",
                    (x_value, y_value),
                    xytext=(4, 4),
                    textcoords="offset points",
                    fontsize=7,
                )
    axis.set_title("Tumor exposure and cumulative cargo-release trade-off")
    axis.set_xlabel("Tumor-associated carrier AUC (normalized carrier h)")
    axis.set_ylabel("Cumulative cargo release (normalized cargo)")
    axis.legend(title="Carrier design")
    axis.grid(alpha=0.25)
    figure.subplots_adjust(bottom=0.24)
    figure.text(
        0.01,
        0.01,
        "Tumor AUC measures exposure/retention; cumulative release measures "
        "cell-processing throughput.\n"
        "They are distinct processes, not efficacy or clinical outcomes. "
        "Exploratory normalized outputs; not clinical predictions.",
        fontsize=7.5,
    )
    return _save_deterministic_svg(figure, path)


def write_svg_timecourse(
    result: SimulationResult,
    path: str | Path,
    *,
    states: tuple[str, ...] = ("A_blood", "A_liver", "A_spleen", "A_tumor_free", "A_internalized", "A_cargo"),
    title: str | None = None,
) -> Path:
    """Write a simple labeled SVG time-course figure."""

    output_path = Path(path)
    output_path.parent.mkdir(parents=True, exist_ok=True)
    width, height = 900, 560
    left, right, top, bottom = 80, 190, 60, 80
    plot_w = width - left - right
    plot_h = height - top - bottom
    max_y = max(float(np.max(result.state(name))) for name in states)
    max_y = max(max_y, 1.0e-12)
    max_t = float(result.t[-1]) if result.t[-1] else 1.0
    colors = ["#1f77b4", "#d62728", "#2ca02c", "#9467bd", "#ff7f0e", "#17becf"]

    lines: list[str] = []
    legends: list[str] = []
    for idx, state_name in enumerate(states):
        color = colors[idx % len(colors)]
        xs = left + plot_w * result.t / max_t
        ys = top + plot_h * (1.0 - result.state(state_name) / max_y)
        lines.append(
            f'<polyline points="{_polyline(zip(xs, ys))}" fill="none" '
            f'stroke="{color}" stroke-width="2" />'
        )
        legend_y = top + idx * 24
        legends.append(
            f'<line x1="{width - right + 20}" y1="{legend_y}" '
            f'x2="{width - right + 48}" y2="{legend_y}" stroke="{color}" '
            f'stroke-width="3" />'
            f'<text x="{width - right + 56}" y="{legend_y + 5}" '
            f'font-size="13">{state_name}</text>'
        )

    caption = (
        "Exploratory normalized amounts; values are not clinical predictions."
    )
    figure_title = title or f"{result.design.name} carrier and cargo time courses"
    svg = f'''<svg xmlns="http://www.w3.org/2000/svg" width="{width}" height="{height}" viewBox="0 0 {width} {height}">
  <rect width="100%" height="100%" fill="white" />
  <text x="{width / 2:.0f}" y="30" text-anchor="middle" font-size="20" font-family="Arial">{figure_title}</text>
  <line x1="{left}" y1="{top + plot_h}" x2="{left + plot_w}" y2="{top + plot_h}" stroke="black" />
  <line x1="{left}" y1="{top}" x2="{left}" y2="{top + plot_h}" stroke="black" />
  <text x="{left + plot_w / 2:.0f}" y="{height - 25}" text-anchor="middle" font-size="15" font-family="Arial">Time (h)</text>
  <text x="24" y="{top + plot_h / 2:.0f}" text-anchor="middle" font-size="15" font-family="Arial" transform="rotate(-90 24 {top + plot_h / 2:.0f})">Normalized amount</text>
  <text x="{left}" y="{top + plot_h + 22}" text-anchor="middle" font-size="12">0</text>
  <text x="{left + plot_w}" y="{top + plot_h + 22}" text-anchor="middle" font-size="12">{max_t:.0f}</text>
  <text x="{left - 10}" y="{top + plot_h + 4}" text-anchor="end" font-size="12">0</text>
  <text x="{left - 10}" y="{top + 4}" text-anchor="end" font-size="12">{max_y:.3g}</text>
  {''.join(lines)}
  {''.join(legends)}
  <text x="{left}" y="{height - 5}" font-size="12" font-family="Arial">{caption}</text>
</svg>
'''
    output_path.write_text(svg, encoding="utf-8")
    return output_path

