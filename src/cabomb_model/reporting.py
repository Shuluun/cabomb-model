"""Report generation helpers."""

from __future__ import annotations

from html import escape
from pathlib import Path

from cabomb_model.analysis import run_all_analyses
from cabomb_model.io import ensure_dir
from cabomb_model.metrics import summary_metrics
from cabomb_model.model_modes import ModelModes, V1_1_MODES
from cabomb_model.parameters import default_model_parameters


def build_markdown_report(
    output_dir: str | Path = "outputs",
    report_path: str | Path = "reports/model_report.md",
    *,
    modes: ModelModes = V1_1_MODES,
) -> Path:
    """Generate the model report from regenerated analysis outputs."""

    results = run_all_analyses(output_dir, modes=modes)
    baseline = results["baseline"]
    baseline_result = baseline["result"] if isinstance(baseline, dict) else None
    params = default_model_parameters()
    baseline_summary = baseline_result and summary_metrics(baseline_result, params)
    design_rows = results["design_comparison"]
    scenario_rows = results["scenario_comparison"]
    alpha_rows = results["ocn_accessibility_sweep"]
    sensitivity = results["local_sensitivity"]
    monte_carlo_summary = results["monte_carlo_design_advantage_summary"]
    report = Path(report_path)
    ensure_dir(report.parent)
    if modes == V1_1_MODES:
        report_version = "V1.1"
    elif (
        modes.mps_model == "saturable_v1_2"
        and modes.targeting_model == "sequential_v1_2"
    ):
        report_version = "V1.2"
    else:
        report_version = "mixed-mode"
    output_root = Path(output_dir).as_posix().rstrip("/")
    lines = [
        f"# CABOMB {report_version} Model Report",
        "",
        "Generated deterministically from the current analysis pipeline; no wall-clock "
        "timestamp is embedded so rerunning the report does not dirty the worktree.",
        "",
        "## Scope",
        "",
        "This is a ten-state PBPK-inspired ODE model for exploratory design "
        "comparison. It is not a full PBPK model and not a clinical dose-prediction model.",
        "",
        "## Active model modes",
        "",
        "| Mechanism dimension | Selected model |",
        "| --- | --- |",
        f"| `mps_model` | `{modes.mps_model}` |",
        f"| `targeting_model` | `{modes.targeting_model}` |",
        "",
        "## Frozen model",
        "",
        "The implementation follows `docs/model_specification.md` and "
        "`docs/equations.md`: carrier and cargo ledgers are separate, OCN only "
        "affects matrix binding/retention, and YSA only affects EphA2 binding "
        "and internalization after tumor-region arrival.",
        "",
        "## Equations and ledgers",
        "",
        "The executable state vector contains circulating carrier, liver-associated "
        "carrier, spleen-associated carrier, tumor free carrier, OCN-bound carrier, "
        "YSA-EphA2-bound carrier, internalized carrier, eliminated/degraded carrier, "
        "active released cargo, and lost or inactive cargo. Carrier mass balance "
        "sums only carrier states; cargo accounting separately tracks payload still "
        "inside active carrier plus active and lost released cargo.",
        "",
        "`k_entry_eff` is computed from an effective baseline entry rate multiplied "
        "by EPR, permeability, and a dimensionless IFP factor. Blood flow, pore "
        "diameter, IFP, and permeability are not inserted directly as first-order "
        "ODE rates.",
        "",
        "## State stock, exposure AUC, and cumulative flux",
        "",
        "State stock is normalized amount present at one time, such as total "
        "tumor-associated carrier or active `A_cargo(t)` remaining after cargo loss.",
        "",
        "Exposure AUC is the time integral of a state stock. `tumor_auc` therefore "
        "has units of normalized carrier amount times hours.",
        "",
        "Cumulative flux is the time integral of process throughput. Forward binding "
        "and transfer totals are gross event counts and can include repeated events. "
        "Gross `cumulative_cargo_release` is not the same quantity as active "
        "`A_cargo(t)` remaining after loss.",
        "",
        "`C_blood = A_blood / Vblood` is normalized carrier per apparent volume, "
        "not an absolute mass concentration and not mg/L.",
        "",
        "## Assumptions and limitations",
        "",
        "The executable baseline uses normalized carrier and cargo units because "
        "the actual wet-lab dose, carrier size distribution, cargo loading, and "
        "release yield are missing. Kinetic constants without supplied measurements "
        "are marked as exploratory assumptions in `data/processed/parameters_master.csv`.",
        "",
        "Unsupported claims are intentionally excluded: docking scores are not "
        "converted into experimental `Kd`, and outputs are not presented as human "
        "clinical predictions.",
        "",
        "Unmodeled organ passage is not automatically terminal elimination; the "
        "current model represents only its named compartments and explicit terminal "
        "loss fluxes.",
        "",
        "## Generated outputs",
        "",
        f"- `{output_root}/tables/baseline_timecourse.csv`",
        f"- `{output_root}/tables/baseline_summary.csv`",
        f"- `{output_root}/tables/design_comparison.csv`",
        f"- `{output_root}/tables/scenario_comparison.csv`",
        f"- `{output_root}/tables/local_sensitivity.csv`",
        f"- `{output_root}/tables/monte_carlo.csv`",
        f"- `{output_root}/tables/ocn_accessibility_sweep.csv`",
        f"- `{output_root}/tables/monte_carlo_design_advantage_summary.csv`",
        f"- `{output_root}/figures/baseline_timecourse.svg`",
        f"- `{output_root}/figures/cumulative_cargo_release_vs_alpha.svg`",
        f"- `{output_root}/figures/tumor_auc_vs_cumulative_cargo_release.svg`",
        f"- `{output_root}/figures/parameters_master_preview.png` "
        "(registry preview generated from `data/processed/parameters_master.csv`)",
        "",
    ]
    if modes.mps_model == "linear_v1":
        lines.extend(
            [
                "## Linear MPS mechanism",
                "",
                "Liver and spleen uptake use the accepted first-order V1 rates. "
                "Intact release from the organ-associated compartments is zero in "
                "this mode, while organ degradation is an irreversible transfer to "
                "the terminal carrier ledger.",
                "",
            ]
        )
    else:
        lines.extend(
            [
                "## Saturable MPS mechanism",
                "",
                "Liver and spleen uptake use exploratory Michaelis-Menten functions. "
                "Saturable uptake and intact MPS release are unvalidated CABOMB assumptions; "
                "MPS release is one candidate mechanism that could contribute to a slow terminal phase.",
                "",
                "Organ degradation remains an irreversible transfer to the terminal "
                "carrier ledger; intact release instead returns the same carrier amount "
                "to blood.",
                "",
            ]
        )
    if modes.targeting_model == "accessibility_v1_1":
        lines.extend(
            [
                "## Accessibility targeting mechanism",
                "",
                "`alpha_ocn_cell_access` is a dimensionless, phenomenological and unvalidated "
                "assumption about effective YSA-EphA2 cellular access retained by OCN-bound "
                "carrier. It is not a measured probability. The OCN-to-YSA flux is an equal "
                "subtraction from `A_ocn_bound` and addition to `A_ysa_bound`; no "
                "`A_dual_bound` state exists.",
                "",
                "Alpha is active only for the OCN-YSA design and does not change blood flow, "
                "EPR, permeability, IFP, or tumor entry. Its baseline is zero so V1.1 "
                "reproduces accepted V1 trajectories and original metrics. No parameter was "
                "tuned to guarantee that OCN-YSA outperforms YSA-only.",
                "",
            ]
        )
    else:
        lines.extend(
            [
                "## Sequential targeting mechanism",
                "",
                "The sequential OCN-to-YSA transfer is a sequential structure assumption: "
                "the OCN-YSA design follows tumor-free carrier to OCN-bound carrier to "
                "YSA-bound carrier before internalization. `k_cross` is not a measured "
                "probability or a validated receptor-binding constant. Fixed `R_EphA2_eff` "
                "is not a dynamic receptor model and does not represent receptor saturation.",
                "",
            ]
        )
        lines = [
            line
            for line in lines
            if "ocn_accessibility_sweep.csv" not in line
            and "cumulative_cargo_release_vs_alpha.svg" not in line
            and "tumor_auc_vs_cumulative_cargo_release.svg" not in line
        ]
    lines.extend(["## Baseline result", ""])
    if baseline_summary:
        lines.extend(
            [
                "| Metric | Value |",
                "| --- | --- |",
                f"| Design | `{baseline_summary['design']}` |",
                f"| Scenario | `{baseline_summary['scenario']}` |",
                f"| Solver | `{baseline_summary['solver_method']}` |",
                f"| Peak tumor-associated fraction | {float(baseline_summary['max_tumor_fraction']):.6g} |",
                f"| Tumor AUC | {float(baseline_summary['tumor_auc']):.6g} |",
                f"| Cumulative internalization | {float(baseline_summary['cumulative_internalization']):.6g} |",
                f"| Cumulative cargo release | {float(baseline_summary['cumulative_cargo_release']):.6g} |",
                f"| Cumulative OCN-to-YSA transfer | {float(baseline_summary['cumulative_ocn_to_ysa_transfer']):.6g} |",
                f"| Final internalized fraction | {float(baseline_summary['internalized_fraction_final']):.6g} |",
                f"| Final active released cargo fraction | {float(baseline_summary['released_cargo_fraction_final']):.6g} |",
                f"| Maximum carrier residual | {float(baseline_summary['max_abs_carrier_residual']):.3e} |",
                f"| Maximum cargo residual | {float(baseline_summary['max_abs_cargo_residual']):.3e} |",
                "",
            ]
        )
    else:
        lines.extend(["Baseline summary was not available.", ""])
    lines.extend(
        [
            "## Four-design comparison",
            "",
            "| Design | Peak tumor-associated fraction | Tumor AUC | Cumulative internalization | Cumulative cargo release | OCN-to-YSA transfer |",
            "| --- | ---: | ---: | ---: | ---: | ---: |",
        ]
    )
    if isinstance(design_rows, list):
        for row in design_rows:
            lines.append(
                f"| `{row['design']}` | {float(row['max_tumor_fraction']):.6g} | "
                f"{float(row['tumor_auc']):.6g} | "
                f"{float(row['cumulative_internalization']):.6g} | "
                f"{float(row['cumulative_cargo_release']):.6g} | "
                f"{float(row['cumulative_ocn_to_ysa_transfer']):.6g} |"
            )
    if modes.mps_model == "saturable_v1_2" and baseline_summary:
        lines.extend(
            [
                "",
                "## V1.2 cumulative MPS processes",
                "",
                "These gross cumulative fluxes are normalized carrier amounts. Uptake, intact "
                "release back to blood, and irreversible degradation are distinct processes.",
                "",
                "| Organ | Cumulative uptake | Cumulative intact release | Cumulative degradation |",
                "| --- | ---: | ---: | ---: |",
                f"| Liver | {float(baseline_summary['cumulative_liver_mps_uptake']):.6g} | {float(baseline_summary['cumulative_liver_mps_release']):.6g} | {float(baseline_summary['cumulative_liver_mps_degradation']):.6g} |",
                f"| Spleen | {float(baseline_summary['cumulative_spleen_mps_uptake']):.6g} | {float(baseline_summary['cumulative_spleen_mps_release']):.6g} | {float(baseline_summary['cumulative_spleen_mps_degradation']):.6g} |",
                "",
            ]
        )
    lines.extend(
        [
            "",
            "## EPR and MPS scenarios",
            "",
            "Low, medium, and high EPR and MPS settings are generated for all four "
            "carrier designs. These are uncertainty scenarios, not measured "
            "osteosarcoma delivery rates.",
            "",
        ]
    )
    if isinstance(scenario_rows, list) and scenario_rows:
        tumor_aucs = [float(row["tumor_auc"]) for row in scenario_rows]
        lines.extend(
            [
                f"Across generated scenario/design rows, tumor AUC ranges from "
                f"{min(tumor_aucs):.6g} to {max(tumor_aucs):.6g} in normalized units.",
                "",
            ]
        )
    if modes.targeting_model == "accessibility_v1_1":
        lines.extend(
            [
                "## OCN accessibility sweep",
                "",
            ]
        )
    if modes.targeting_model == "accessibility_v1_1" and isinstance(alpha_rows, list) and alpha_rows:
        combo_rows = [row for row in alpha_rows if row["design"] == "ocn-ysa"]
        lines.extend(
            [
                "The deterministic sweep contains all four designs at alpha 0, 0.25, "
                "0.5, 0.75, and 1. Only OCN-YSA is mechanistically responsive; no "
                "monotonic benefit or superiority rule is imposed.",
                "",
                "| Alpha | OCN-YSA tumor AUC | OCN-YSA cumulative cargo release |",
                "| ---: | ---: | ---: |",
            ]
        )
        for row in combo_rows:
            lines.append(
                f"| {float(row['alpha_ocn_cell_access']):.2f} | "
                f"{float(row['tumor_auc']):.6g} | "
                f"{float(row['cumulative_cargo_release']):.6g} |"
            )
    elif modes.targeting_model == "accessibility_v1_1":
        lines.append("OCN accessibility sweep results were not available.")
    lines.extend(["", "## Multi-target local sensitivity", ""])
    if isinstance(sensitivity, list) and sensitivity:
        target_names = list(dict.fromkeys(str(row["target_metric"]) for row in sensitivity))
        for target_name in target_names:
            lines.extend(
                [
                    f"### `{target_name}`",
                    "",
                    "| Rank | Parameter | Method | Abs normalized sensitivity |",
                    "| ---: | --- | --- | ---: |",
                ]
            )
            target_rows = [
                row for row in sensitivity if row["target_metric"] == target_name
            ][:5]
            for idx, row in enumerate(target_rows, start=1):
                value = row["abs_normalized_sensitivity"]
                formatted = "not defined" if value is None else f"{float(value):.4g}"
                lines.append(
                    f"| {idx} | `{row['parameter_symbol']}` | "
                    f"`{row['sensitivity_method']}` | {formatted} |"
                )
            lines.append("")
        sensitivity_boundary = (
            "Boundary one-sided alpha sensitivities and central relative "
            "sensitivities are different local exploratory estimates and are not "
            "assumed to be biologically interchangeable."
            if modes.targeting_model == "accessibility_v1_1"
            else "Local sensitivities cover only active parameters for the selected "
            f"`mps_model={modes.mps_model}` and "
            f"`targeting_model={modes.targeting_model}` mechanisms."
        )
        lines.extend([sensitivity_boundary, ""])
    else:
        lines.extend(["Sensitivity results were not available.", ""])
    lines.extend(["## Paired four-design Monte Carlo comparison", ""])
    if isinstance(monte_carlo_summary, dict):
        n = int(monte_carlo_summary["sample_count"])
        lines.extend(
            [
                "Each sampled active-parameter vector is reused across all four "
                "designs. Strict OCN-YSA wins over YSA-only are summarized below; "
                "ties do not count as wins.",
                "",
                "| Comparison | Win count | Sample count | Fraction |",
                "| --- | ---: | ---: | ---: |",
                f"| Cumulative cargo release | {int(monte_carlo_summary['cargo_release_win_count'])} | {n} | {float(monte_carlo_summary['cargo_release_win_fraction']):.6g} |",
                f"| Tumor AUC | {int(monte_carlo_summary['tumor_auc_win_count'])} | {n} | {float(monte_carlo_summary['tumor_auc_win_fraction']):.6g} |",
                f"| Both metrics | {int(monte_carlo_summary['both_win_count'])} | {n} | {float(monte_carlo_summary['both_win_fraction']):.6g} |",
                "",
                "These are paired-sample fractions under exploratory parameter distributions, not probabilities of clinical success.",
            ]
        )
    else:
        lines.append("Paired Monte Carlo summary was not available.")
    lines.extend(
        [
            "",
            "## Wet-lab implications",
            "",
            "Grouped sensitivity results should be used to choose measurements for "
            "the specific target of interest.",
            "",
        ]
    )
    if modes.targeting_model == "accessibility_v1_1":
        lines.extend(
            [
                "OCN-bound cellular accessibility and `alpha_ocn_cell_access` are "
                "measurement priorities because the accessibility term is exploratory, "
                "not because the modeled transfer has been experimentally validated.",
                "",
            ]
        )
    else:
        lines.extend(
            [
                "`k_cross` and the strict OCN-to-YSA sequence assumption are "
                "measurement priorities. Experiments should test the transfer rate "
                "and whether matrix retention precedes EphA2 engagement, rather than "
                "treating the imposed order as established biology.",
                "",
            ]
        )
    lines.extend(
        [
            "## Interpretation boundary",
            "",
            "All outputs are normalized exploratory values because key CABOMB-specific "
            "parameters remain missing. Results should guide wet-lab measurement "
            "priority and design logic, not clinical claims.",
            "",
        ]
    )
    report.write_text("\n".join(lines), encoding="utf-8")
    return report


def build_html_report(
    markdown_path: str | Path,
    html_path: str | Path,
    *,
    title: str = "CABOMB V1.1 Model Report",
) -> Path:
    """Generate a small standalone rendered HTML report from markdown text."""

    md = Path(markdown_path).read_text(encoding="utf-8")
    body = _markdown_to_html(md)
    html = f"""<!doctype html>
<html lang="en">
<head>
  <meta charset="utf-8" />
  <title>{title}</title>
  <style>
    body {{ font-family: Arial, sans-serif; max-width: 980px; margin: 40px auto; line-height: 1.5; color: #1F2937; }}
    h1, h2 {{ color: #0F3B57; }}
    table {{ border-collapse: collapse; width: 100%; margin: 1rem 0; }}
    th, td {{ border: 1px solid #D1D5DB; padding: 0.45rem 0.6rem; text-align: left; }}
    th {{ background: #EAF2F8; }}
    code {{ background: #F3F4F6; padding: 0.1rem 0.25rem; border-radius: 3px; }}
  </style>
</head>
<body>
{body}
</body>
</html>
"""
    out = Path(html_path)
    ensure_dir(out.parent)
    out.write_text(html, encoding="utf-8")
    return out


def _inline_markdown(text: str) -> str:
    """Render a tiny, dependency-free subset of inline markdown."""

    parts = text.split("`")
    rendered: list[str] = []
    for idx, part in enumerate(parts):
        if idx % 2:
            rendered.append(f"<code>{escape(part)}</code>")
        else:
            rendered.append(escape(part))
    return "".join(rendered)


def _markdown_to_html(markdown: str) -> str:
    """Render the report's limited markdown subset to HTML."""

    lines = markdown.splitlines()
    html: list[str] = []
    in_ul = False
    in_table = False
    table_align: list[str] = []

    def close_blocks() -> None:
        nonlocal in_ul, in_table, table_align
        if in_ul:
            html.append("</ul>")
            in_ul = False
        if in_table:
            html.append("</tbody></table>")
            in_table = False
            table_align = []

    idx = 0
    while idx < len(lines):
        line = lines[idx]
        stripped = line.strip()
        if not stripped:
            close_blocks()
            idx += 1
            continue
        if stripped.startswith("### "):
            close_blocks()
            html.append(f"<h3>{_inline_markdown(stripped[4:])}</h3>")
        elif stripped.startswith("## "):
            close_blocks()
            html.append(f"<h2>{_inline_markdown(stripped[3:])}</h2>")
        elif stripped.startswith("# "):
            close_blocks()
            html.append(f"<h1>{_inline_markdown(stripped[2:])}</h1>")
        elif stripped.startswith("- "):
            if in_table:
                close_blocks()
            if not in_ul:
                html.append("<ul>")
                in_ul = True
            html.append(f"<li>{_inline_markdown(stripped[2:])}</li>")
        elif stripped.startswith("|") and stripped.endswith("|"):
            cells = [cell.strip() for cell in stripped.strip("|").split("|")]
            next_line = lines[idx + 1].strip() if idx + 1 < len(lines) else ""
            is_header = next_line.startswith("|") and set(next_line.replace("|", "").replace(":", "").replace("-", "").strip()) == set()
            is_separator = set(stripped.replace("|", "").replace(":", "").replace("-", "").strip()) == set()
            if is_separator:
                idx += 1
                continue
            if not in_table:
                close_blocks()
                html.append("<table>")
                in_table = True
            if is_header:
                table_align = cells
                html.append("<thead><tr>" + "".join(f"<th>{_inline_markdown(cell)}</th>" for cell in cells) + "</tr></thead><tbody>")
                idx += 1
            else:
                tag = "td" if table_align else "th"
                html.append("<tr>" + "".join(f"<{tag}>{_inline_markdown(cell)}</{tag}>" for cell in cells) + "</tr>")
        else:
            close_blocks()
            html.append(f"<p>{_inline_markdown(stripped)}</p>")
        idx += 1
    close_blocks()
    return "\n".join(f"  {line}" for line in html)
