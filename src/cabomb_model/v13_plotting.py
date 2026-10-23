"""SVG plotting helpers for the V1.3 explicit-circulation mode."""

from __future__ import annotations

from pathlib import Path

import matplotlib as mpl

mpl.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

from cabomb_model.v13_circulation import V13Result


def write_v13_ablation_timecourses(
    results: dict[str, V13Result],
    fluxes: dict[str, dict[str, np.ndarray]],
    params: dict[str, float],
    path: str | Path,
) -> Path:
    """Compare stocks and cumulative fluxes without assuming a winning design."""
    from cabomb_model.v13_circulation import v13_tumor_associated

    output = Path(path)
    output.parent.mkdir(parents=True, exist_ok=True)
    colors = ("#777777", "#4477AA", "#228833", "#AA3377")
    labels = ("Unmodified", "OCN-only", "YSA-only", "OCN-YSA")
    styles = (":", "--", "-.", "-")
    dose = params["D0"]
    cargo0 = dose * params["Y_cargo"]
    with mpl.rc_context({"svg.fonttype": "none", "font.size": 9, "axes.spines.top": False, "axes.spines.right": False}):
        fig, axes = plt.subplots(2, 2, figsize=(10, 7))
        titles = ("a  Tumor-associated carrier stock", "b  Total cumulative internalization",
                  "c  Cumulative intracellular cargo release", "d  Active cargo remaining")
        for name, label, color, style in zip(results, labels, colors, styles):
            result = results[name]
            cumulative = fluxes[name]
            zero = np.zeros_like(result.t)
            series = (
                v13_tumor_associated(result) / dose if dose else zero,
                (cumulative["internalization"] + cumulative["nonspecific_internalization"]) / dose if dose else zero,
                cumulative["deg_int"] * params["Y_cargo"] / cargo0 if cargo0 else zero,
                result.state("A_cargo") / cargo0 if cargo0 else zero,
            )
            for ax, y in zip(axes.flat, series):
                ax.plot(result.t, y, label=label, color=color, linestyle=style, linewidth=1.7)
        for i, (ax, title) in enumerate(zip(axes.flat, titles)):
            ax.set_title(title, loc="left", fontsize=10)
            ax.set_xlabel("Time (h)")
            ax.set_ylabel("Fraction of initial " + ("carrier dose" if i < 2 else "cargo load"))
            ax.set_xlim(left=0)
            ax.set_ylim(bottom=0)
        axes[0, 0].legend(fontsize=8)
        first = next(iter(results.values()))
        fig.suptitle(f"CABOMB V1.3 | {first.size} | {first.modes.targeting_model}")
        fig.tight_layout(rect=(0, 0.08, 1, 0.95))
        fig.text(0.05, 0.02, "One deterministic exploratory run per design; common parameters and initial dose. No experimental error bars.\n"
                 "Stocks and cumulative release answer different questions. Cargo is not free cytosolic calcium.", fontsize=8)
        fig.savefig(output.with_suffix(".png"), dpi=180, bbox_inches="tight")
        return _save_deterministic_svg(fig, output)

_VASCULAR_COLORS = {
    "A_central_venous": "#1f77b4",
    "A_right_heart": "#d62728",
    "A_lung_vascular": "#2ca02c",
    "A_left_heart": "#9467bd",
    "A_systemic_arterial": "#ff7f0e",
    "A_tumor_vascular": "#17becf",
}
_TUMOR_COLORS = {
    "A_tumor_free": "#1f77b4",
    "A_ocn_bound": "#d62728",
    "A_ysa_bound": "#2ca02c",
    "A_internalized": "#9467bd",
    "A_cargo": "#ff7f0e",
}


def _save_deterministic_svg(figure: plt.Figure, path: str | Path) -> Path:
    """Save and close an SVG without time-varying metadata or random IDs."""

    output = Path(path)
    output.parent.mkdir(parents=True, exist_ok=True)
    with mpl.rc_context({"svg.hashsalt": "cabomb-v1.3"}):
        figure.savefig(
            output,
            format="svg",
            bbox_inches="tight",
            metadata={"Date": None, "Creator": "CABOMB V1.3"},
        )
    svg = output.read_text(encoding="utf-8")
    normalized = "\n".join(line.rstrip(" \t") for line in svg.splitlines()) + "\n"
    output.write_text(normalized, encoding="utf-8", newline="\n")
    plt.close(figure)
    return output


def write_v13_circulation_timecourse(
    result: V13Result,
    path: str | Path,
    *,
    title: str | None = None,
) -> Path:
    """Write the V1.3 two-panel circulation and tumor-targeting time course.

    The left panel resolves the first pulmonary pass and arterial peak on the
    early time grid; the right panel shows the full tumor-targeting trajectory.
    """

    figure, (ax_vascular, ax_tumor) = plt.subplots(1, 2, figsize=(13.0, 5.0))
    early_mask = result.t <= 0.05

    for state_name, color in _VASCULAR_COLORS.items():
        ax_vascular.plot(
            result.t[early_mask], result.state(state_name)[early_mask], color=color, label=state_name
        )
    ax_vascular.set_title("Circulation (first 0.05 h)")
    ax_vascular.set_xlabel("Time (h)")
    ax_vascular.set_ylabel("Normalized carrier amount")
    ax_vascular.legend(fontsize=6.5, loc="upper right")

    for state_name, color in _TUMOR_COLORS.items():
        ax_tumor.plot(result.t, result.state(state_name), color=color, label=state_name)
    ax_tumor.set_title("Tumor targeting and cargo (full horizon)")
    ax_tumor.set_xlabel("Time (h)")
    ax_tumor.set_ylabel("Normalized amount")
    ax_tumor.legend(fontsize=6.5, loc="best")

    figure.suptitle(title or "CABOMB V1.3 explicit-circulation time courses")
    figure.subplots_adjust(bottom=0.22)
    figure.text(
        0.01,
        0.01,
        "V1.3 normalized exploratory amounts; the arterial and tumor-vascular "
        "curves are carrier states, not concentrations.\n"
        "Not clinical predictions; not a validated organ-distribution model.",
        fontsize=7.5,
    )
    return _save_deterministic_svg(figure, path)


def write_v13_mechanism_diagram(path: str | Path) -> Path:
    """Write a hand-authored SVG of the V1.3 circulation and targeting loop."""

    width, height = 1160, 640
    boxes: list[tuple[str, float, float, float, float, str]] = [
        ("IV bolus", 40, 280, 120, 48, "#fde68a"),
        ("Central venous", 190, 268, 150, 72, "#dbeafe"),
        ("Right heart", 380, 268, 140, 72, "#dbeafe"),
        ("Lung vascular", 560, 150, 150, 72, "#bbf7d0"),
        ("Left heart", 560, 300, 140, 72, "#dbeafe"),
        ("Systemic arterial", 740, 268, 160, 72, "#dbeafe"),
        ("Liver vascular", 950, 90, 150, 60, "#fef9c3"),
        ("Spleen vascular", 950, 190, 150, 60, "#fef9c3"),
        ("Tumor-feeding vessel", 950, 290, 170, 60, "#fecaca"),
        ("Rest vascular", 950, 390, 150, 60, "#fef9c3"),
        ("Tumor interstitium\n(free carrier)", 950, 470, 190, 66, "#fee2e2"),
        ("OCN-bound", 760, 470, 130, 60, "#fbcfe8"),
        ("YSA-EphA2-bound", 540, 470, 160, 60, "#e9d5ff"),
        ("Internalized", 330, 470, 140, 60, "#ddd6fe"),
        ("Cargo release", 160, 470, 130, 60, "#c7d2fe"),
        ("Lung retention\n(terminal)", 560, 40, 160, 56, "#f1f5f9"),
        ("Clearance / degradation\n(terminal)", 330, 60, 190, 56, "#f1f5f9"),
    ]
    arrows = [
        (100, 304, 190, 304, "→"),
        (340, 304, 380, 304, "→"),
        (450, 280, 560, 200, "→"),
        (560, 250, 560, 300, "→"),
        (700, 336, 740, 318, "→"),
        (820, 320, 950, 300, "→"),
        (820, 300, 950, 130, "→"),
        (820, 280, 950, 220, "→"),
        (820, 340, 950, 420, "→"),
        (950, 360, 950, 470, "extravasation\n(EPR/permeability)"),
        (950, 536, 890, 500, "→"),
        (890, 500, 760, 500, "OCN on"),
        (630, 500, 540, 500, "OCN→YSA"),
        (380, 500, 330, 500, "→"),
        (260, 500, 160, 500, "→"),
        (560, 150, 560, 96, "lung retention"),
        (820, 300, 520, 90, "venous return"),
    ]
    body: list[str] = [
        f'<svg xmlns="http://www.w3.org/2000/svg" width="{width}" height="{height}" viewBox="0 0 {width} {height}">',
        '<rect width="100%" height="100%" fill="white" />',
        f'<text x="{width / 2:.0f}" y="28" text-anchor="middle" font-size="20" font-family="Arial">CABOMB V1.3 explicit circulation and targeting</text>',
    ]
    for label, x, y, w, h, fill in boxes:
        lines = label.split("\n")
        body.append(
            f'<rect x="{x:.0f}" y="{y:.0f}" width="{w:.0f}" height="{h:.0f}" '
            f'fill="{fill}" stroke="#334155" stroke-width="1.5" />'
        )
        text_y = y + h / 2 - (len(lines) - 1) * 7
        for idx, line in enumerate(lines):
            body.append(
                f'<text x="{x + w / 2:.0f}" y="{text_y + idx * 14:.0f}" '
                f'text-anchor="middle" font-size="12" font-family="Arial">{line}</text>'
            )
    for x1, y1, x2, y2, label in arrows:
        body.append(
            f'<line x1="{x1:.0f}" y1="{y1:.0f}" x2="{x2:.0f}" y2="{y2:.0f}" '
            f'stroke="#334155" stroke-width="1.5" />'
        )
        body.append(
            f'<text x="{(x1 + x2) / 2:.0f}" y="{(y1 + y2) / 2 - 4:.0f}" '
            f'text-anchor="middle" font-size="9" font-family="Arial">{label}</text>'
        )
    body.append(
        '<line x1="1045" y1="455" x2="640" y2="455" stroke="#7c3aed" stroke-width="1.5" '
        'stroke-dasharray="6,4" />'
    )
    body.append(
        '<text x="845" y="448" text-anchor="middle" font-size="9" font-family="Arial" '
        'fill="#7c3aed">parallel_dual_access_v1_3 only: direct free-to-YSA '
        '(alpha_ysa_direct_access)</text>'
    )
    body.append(
        f'<text x="{width / 2:.0f}" y="{height - 10:.0f}" text-anchor="middle" '
        f'font-size="11" font-family="Arial">Schematic only; OCN/YSA act post-extravasation '
        f'and never change cardiac output, pulmonary transit, or initial tumor supply.</text>'
    )
    body.append("</svg>")

    output = Path(path)
    output.parent.mkdir(parents=True, exist_ok=True)
    output.write_text("\n".join(body) + "\n", encoding="utf-8")
    return output


def write_v13_mass_conservation(
    result: V13Result,
    params: Mapping[str, float],
    path: str | Path,
) -> Path:
    """Plot carrier and cargo mass ledgers against their conserved totals."""

    from cabomb_model.v13_circulation import v13_cargo_total, v13_carrier_total

    dose = params["D0"]
    carrier = v13_carrier_total(result)
    cargo = v13_cargo_total(result, params)
    figure, (ax_carrier, ax_cargo) = plt.subplots(1, 2, figsize=(12.0, 4.5))

    ax_carrier.plot(result.t, carrier, label="carrier total", color="#1f77b4")
    ax_carrier.axhline(dose, color="k", linestyle="--", linewidth=1.0, label=f"D0 = {dose:g}")
    ax_carrier.set_title("Carrier mass ledger")
    ax_carrier.set_xlabel("Time (h)")
    ax_carrier.set_ylabel("Normalized carrier amount")
    ax_carrier.legend(fontsize=8)

    cargo0 = params["Y_cargo"] * dose
    ax_cargo.plot(result.t, cargo, label="cargo total", color="#2ca02c")
    ax_cargo.axhline(cargo0, color="k", linestyle="--", linewidth=1.0, label=f"Y_cargo*D0 = {cargo0:g}")
    ax_cargo.set_title("Cargo mass ledger")
    ax_cargo.set_xlabel("Time (h)")
    ax_cargo.set_ylabel("Normalized cargo amount")
    ax_cargo.legend(fontsize=8)

    figure.suptitle("CABOMB V1.3 mass conservation")
    figure.text(
        0.01,
        0.01,
        "Carrier and cargo ledgers stay flat at their conserved totals within solver "
        "tolerance; this is software verification, not a biological result.",
        fontsize=7.5,
    )
    return _save_deterministic_svg(figure, path)


def write_v13_mechanism_comparison(rows: list[dict], path: str | Path) -> Path:
    """Bar plot of sequential vs parallel targeting for the four designs."""

    designs = [row["design"] for row in rows if row.get("size") == "100nm" and row.get("targeting_model") == "sequential_v1_2"]
    metrics = ("tumor_auc", "internalized_fraction_final", "released_cargo_fraction_final")
    figure, axes = plt.subplots(1, 3, figsize=(14.0, 4.2))
    width = 0.36
    x = np.arange(len(designs))
    for ax, metric in zip(axes, metrics):
        sequential = [
            float(row[metric])
            for row in rows
            if row.get("size") == "100nm" and row.get("targeting_model") == "sequential_v1_2"
        ]
        parallel = [
            float(row[metric])
            for row in rows
            if row.get("size") == "100nm" and row.get("targeting_model") == "parallel_dual_access_v1_3"
        ]
        ax.bar(x - width / 2, sequential, width, label="sequential_v1_2", color="#94a3b8")
        ax.bar(x + width / 2, parallel, width, label="parallel_dual_access_v1_3", color="#7c3aed")
        ax.set_xticks(x)
        ax.set_xticklabels(designs, rotation=20, fontsize=7)
        ax.set_title(metric, fontsize=8)
        ax.ticklabel_format(axis="y", style="scientific", scilimits=(0, 0))
    axes[0].legend(fontsize=7)
    figure.suptitle("Sequential vs parallel dual-access targeting (100 nm, OCN-YSA-unmodified ablations)")
    figure.text(
        0.01,
        0.01,
        "Parallel mode adds a direct free-to-YSA route gated by alpha_ysa_direct_access; "
        "both modes leave cardiac output, EPR, permeability, and initial tumor supply unchanged.",
        fontsize=7.5,
    )
    return _save_deterministic_svg(figure, path)
