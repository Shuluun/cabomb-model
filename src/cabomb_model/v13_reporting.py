"""Markdown report generation for the optional V1.3 explicit-circulation mode."""

from __future__ import annotations

from pathlib import Path

from cabomb_model.io import ensure_dir
from cabomb_model.v13_analysis import run_v13_all_analyses

# Metric keys shown in the V1.3 comparison tables.
_DESIGN_KEYS = (
    "design",
    "size",
    "max_tumor_fraction",
    "tumor_auc",
    "internalized_fraction_final",
    "released_cargo_fraction_final",
    "net_extravasated_fraction",
)


def _fmt(value: object, digits: int = 6) -> str:
    try:
        return f"{float(value):.{digits}g}"
    except (TypeError, ValueError):
        return str(value)


def write_v13_ablation_summary(rows: list[dict[str, object]], path: str | Path) -> Path:
    """Write a scoped Chinese report using already computed four-design rows."""
    labels = {"unmodified": "无修饰", "ocn": "OCN-only", "ysa": "YSA-only", "ocn-ysa": "OCN–YSA"}
    by_design = {str(row["design"]): row for row in rows}
    first = rows[0]
    lines = [
        "# CABOMB 四组计算消融：条件情景汇报", "",
        f"共同设置：{first['size']} 设计情景（不是实测制备尺寸），{first['targeting_model']}，"
        f"EPR/MPS medium，{float(first['time_final_h']):g} h，D0={first['initial_dose_normalized']} 归一化载体量。",
        "四组共用参数与时间网格，仅改变 OCN/YSA 开关。每组一个确定性计算，无实验重复或置信区间。", "",
        "| 设计 | 肿瘤 AUC | 累计特异内化/D0 | 累计非特异内化/D0 | 累计总内化/D0 | 累计释放/初始cargo | 末端活性cargo/初始cargo |",
        "|---|---:|---:|---:|---:|---:|---:|",
    ]
    keys = ("tumor_auc", "cumulative_specific_internalization_fraction", "cumulative_nonspecific_internalization_fraction",
            "cumulative_total_internalization_fraction", "cargo_release_efficiency_gross", "active_cargo_remaining_fraction_final")
    for row in rows:
        lines.append("| " + labels[str(row["design"])] + " | " + " | ".join(_fmt(row[k]) for k in keys) + " |")
    lines += ["", "## 如何理解", "",
              "AUC 是肿瘤相关载体存量的时间积分，单位为归一化载体量·h；其余表列是比例。"
              "累计内化计入特异与非特异两条路径；末端内化存量还要扣除降解。"
              "累计释放只积分内化区室的降解源，末端活性 cargo 还受 cargo 损失影响。", ""]
    for key, label in (("tumor_auc", "肿瘤 AUC"), ("cumulative_total_internalization_fraction", "累计总内化"),
                       ("cargo_release_efficiency_gross", "累计 cargo 释放"), ("active_cargo_remaining_fraction_final", "末端活性 cargo")):
        ref, dual = float(by_design["ysa"][key]), float(by_design["ocn-ysa"][key])
        change = f"{(dual / ref - 1) * 100:+.3f}%（相对变化）" if ref else "不可计算（参照为零）"
        lines.append(f"- OCN–YSA 相对 YSA-only：{label} {change}。")
    lines += ["", "不同指标可能给出不同排序，不能据此宣称双修饰全面更优。参数主要为探索性假设，"
              "结果支持机制比较，尚未经实际样品验证。100 nm 保留为设计基线；200 nm 改良目标和历史微米级样品不代表本次样品测量。",
              "", "## 证据与边界", "",
              "曲线见 figures/four_design_timecourses.svg；原始轨迹、累计通量和参数见 tables/；"
              "run_manifest.json 记录提交、工作区状态、源文件与参数哈希、模式、求解器及完整时间网格。",
              "旧滞留效率从核心汇报排除；末端肿瘤存量和累计外渗分别保留在结果表。"
              "供钙当量是独立条件换算，不从归一化 A_cargo 推导每细胞钙质量。"
              "本节仅汇报四设计基线，未在此重跑可行域扫描，也未构建真实胞质游离钙模型。", ""]
    output = Path(path)
    output.parent.mkdir(parents=True, exist_ok=True)
    output.write_text("\n".join(lines), encoding="utf-8")
    return output


def build_v13_markdown_report(
    output_dir: str | Path = "outputs/v1_3",
    report_path: str | Path = "reports/model_report_v1_3.md",
) -> Path:
    """Generate the V1.3 report from regenerated analysis outputs."""

    results = run_v13_all_analyses(output_dir)
    baseline = results["baseline"]
    baseline_result = baseline["result"] if isinstance(baseline, dict) else None
    params = baseline["params"] if isinstance(baseline, dict) else None

    output_root = Path(output_dir).as_posix().rstrip("/")
    lines: list[str] = [
        "# CABOMB V1.3 Explicit-Circulation Model Report",
        "",
        "Generated deterministically from the current V1.3 analysis pipeline; no "
        "wall-clock timestamp is embedded so rerunning the report does not dirty "
        "the worktree.",
        "",
        "## Scope",
        "",
        "V1.3 is an **optional** mode that replaces the lumped `A_blood` compartment "
        "with an explicit closed circulation (IV bolus -> central venous -> right heart "
        "-> pulmonary -> left heart -> systemic artery -> organ branches -> tumor-feeding "
        "vessel -> interstitium -> venous return). It keeps V1.2's saturable MPS and "
        "sequential OCN-to-YSA targeting unchanged and does not modify or overwrite the "
        "frozen V1.1/V1.2 path, its 128 tests, or its outputs.",
        "",
        "## Active model modes",
        "",
        "| Mechanism dimension | Selected model |",
        "| --- | --- |",
        "| `mps_model` | `saturable_v1_2` |",
        "| `targeting_model` (baseline) | `sequential_v1_2` |",
        "| `targeting_model` (mechanism comparison) | `sequential_v1_2` vs `parallel_dual_access_v1_3` |",
        "| `circulation_model` | `explicit_v1_3` |",
        "",
        "## Interpretation boundary",
        "",
        "V1.3 remains a mechanistic exploration and design-comparison model. It is not a "
        "clinical digital twin, efficacy predictor, human dose calculator, or validated "
        "organ-distribution model. Every V1.3 circulation volume and size factor is a "
        "transparent human-physiology reference or an exploratory assumption, recorded in "
        "`docs/parameter_evidence_v1_3.md`; none is a measured CABOMB value.",
        "",
        "## Baseline result",
        "",
    ]
    if baseline_result is not None and params is not None:
        from cabomb_model.v13_analysis import v13_mass_distribution
        from cabomb_model.v13_circulation import (
            v13_max_abs_cargo_residual,
            v13_max_abs_carrier_residual,
            v13_summary_metrics,
        )

        summary = v13_summary_metrics(baseline_result, params)
        distribution = v13_mass_distribution(baseline_result, params)
        lines.extend(
            [
                "| Metric | Value |",
                "| --- | --- |",
                f"| Design | `{baseline_result.design.name}` |",
                f"| Size scenario | `{baseline_result.size}` |",
                f"| Scenario | `{baseline_result.scenario.name}` |",
                f"| Solver | `{baseline_result.solver_method}` (rtol={baseline_result.rtol:g}, atol={baseline_result.atol:g}) |",
                f"| Peak tumor-associated fraction | {_fmt(summary['max_tumor_fraction'])} |",
                f"| Tumor AUC | {_fmt(summary['tumor_auc'])} |",
                f"| Net extravasated fraction | {_fmt(summary['net_extravasated_fraction'])} |",
                f"| Final internalized fraction | {_fmt(summary['internalized_fraction_final'])} |",
                f"| Final active released cargo fraction | {_fmt(summary['released_cargo_fraction_final'])} |",
                f"| Maximum carrier residual | {_fmt(v13_max_abs_carrier_residual(baseline_result, params))} |",
                f"| Maximum cargo residual | {_fmt(v13_max_abs_cargo_residual(baseline_result, params))} |",
                "",
            ]
        )
        # Final mass distribution table.
        lines.extend(
            [
                "### Final carrier mass distribution",
                "",
                "| Compartment / loss route | Final fraction |",
                "| --- | ---: |",
            ]
        )
        labels = {
            "vascular_fraction_final": "Vascular (all nine vascular states)",
            "liver_fraction_final": "Liver-MPS-associated",
            "spleen_fraction_final": "Spleen-MPS-associated",
            "tumor_associated_fraction_final": "Tumor-associated (free+OCN+YSA+internalized)",
            "internalized_fraction_final": "Internalized",
            "elim_deg_fraction_final": "Terminal eliminated/degraded (total)",
            "cumulative_lung_ret_fraction": "  └ cumulative lung retention over all passes (terminal)",
            "cumulative_blood_clear_fraction": "  └ systemic blood clearance (terminal)",
            "cumulative_liver_elim_fraction": "  └ liver MPS degradation (terminal)",
            "cumulative_spleen_elim_fraction": "  └ spleen MPS degradation (terminal)",
            "cumulative_tumor_deg_fraction": "  └ tumor interstitial degradation (terminal)",
            "cumulative_ocn_loss_fraction": "  └ OCN-bound loss (terminal)",
            "cumulative_deg_int_fraction": "  └ intracellular degradation (terminal)",
        }
        for key, label in labels.items():
            lines.append(f"| {label} | {_fmt(distribution[key])} |")
        lines.append("")
    else:
        lines.extend(["Baseline summary was not available.", ""])

    # Four-design comparison.
    design_rows = results["design_comparison"]
    lines.extend(
        [
            "## Four-design ablation",
            "",
            "| Design | Peak tumor fraction | Tumor AUC | Internalized (final) | Released cargo (final) | Net extravasated |",
            "| --- | ---: | ---: | ---: | ---: | ---: |",
        ]
    )
    if isinstance(design_rows, list):
        for row in design_rows:
            lines.append(
                f"| `{row['design']}` | {_fmt(row['max_tumor_fraction'])} | "
                f"{_fmt(row['tumor_auc'])} | {_fmt(row['internalized_fraction_final'])} | "
                f"{_fmt(row['released_cargo_fraction_final'])} | {_fmt(row['net_extravasated_fraction'])} |"
            )
    lines.append("")

    # Six-size comparison.
    size_rows = results["size_comparison"]
    lines.extend(
        [
            "## Size-scenario comparison (OCN-YSA)",
            "",
            "| Size | Peak tumor fraction | Tumor AUC | Internalized (final) | Released cargo (final) | Net extravasated | Lung retention (cumulative) |",
            "| --- | ---: | ---: | ---: | ---: | ---: | ---: |",
        ]
    )
    if isinstance(size_rows, list):
        for row in size_rows:
            lines.append(
                f"| `{row['size']}` | {_fmt(row['max_tumor_fraction'])} | "
                f"{_fmt(row['tumor_auc'])} | {_fmt(row['internalized_fraction_final'])} | "
                f"{_fmt(row['released_cargo_fraction_final'])} | {_fmt(row['net_extravasated_fraction'])} | "
                f"{_fmt(row['cumulative_lung_retention_fraction'])} |"
            )
    lines.extend(
        [
            "",
            "100 nm is the original design scenario and 200 nm the optimization target; "
            "500/800 nm, 2 um, and 3 um are exploratory extensions from the oral/EM "
            "observation, not validated distribution data. The 2 um and 3 um scenarios "
            "predict heavy lung/MPS loss, so a 100 nm result must not be used to represent "
            "the current experimental microparticles.",
            "",
        ]
    )

    # First-pass and delivery metrics (baseline detailed).
    detailed_rows = results.get("baseline_detailed", [])
    if detailed_rows:
        detailed = detailed_rows[0]
        lines.extend(
            [
                "## First-pass, delivery, and temporal metrics (baseline)",
                "",
                "First-pass metrics come from a single-pass run (venous return diverted to a "
                "sink, no recirculation) and are therefore unique-particle first-event "
                "fractions, not gross throughput. Gross binding metrics count repeated "
                "binding and may exceed 100 %. Total internalization sums specific and "
                "nonspecific paths and is distinct from remaining internalized stock.",
                "",
                "| Metric | Value |",
                "| --- | ---: |",
                f"| First-pass lung arrival fraction | {_fmt(detailed['first_pass_lung_arrival_fraction'])} |",
                f"| First-pass lung retention fraction | {_fmt(detailed['first_pass_lung_retention_fraction'])} |",
                f"| First-pass systemic availability | {_fmt(detailed['first_pass_systemic_availability'])} |",
                f"| Tumor-vessel first-pass arrival fraction | {_fmt(detailed['tumor_vascular_first_pass_arrival_fraction'])} |",
                f"| Tumor first extravasation fraction | {_fmt(detailed['tumor_first_extravasation_fraction'])} |",
                f"| First-pass extravasation efficiency | {_fmt(detailed['first_pass_extravasation_efficiency'])} |",
                f"| Final tumor-associated stock / initial dose | {_fmt(detailed['tumor_fraction_final'])} |",
                f"| Apparent blood half-life (h) | {_fmt(detailed['apparent_blood_half_life_h'])} |",
                f"| OCN binding fraction (gross) | {_fmt(detailed['ocn_binding_fraction_gross'])} |",
                f"| YSA binding fraction (gross) | {_fmt(detailed['ysa_binding_fraction_gross'])} |",
                f"| OCN-to-YSA transfer fraction (gross) | {_fmt(detailed['ocn_to_ysa_transfer_fraction_gross'])} |",
                f"| Total internalization / gross entry (specific + nonspecific) | {_fmt(detailed['internalization_efficiency_gross'])} |",
                f"| Specific internalization / gross entry | {_fmt(detailed['ysa_specific_internalization_fraction_gross'])} |",
                f"| Nonspecific internalization / gross entry | {_fmt(detailed['nonspecific_internalization_fraction_gross'])} |",
                f"| Cargo release efficiency (gross) | {_fmt(detailed['cargo_release_efficiency_gross'])} |",
                f"| Active cargo remaining (final) | {_fmt(detailed['active_cargo_remaining_fraction_final'])} |",
                "",
            ]
        )

    # Sequential vs parallel mechanism comparison.
    mechanism_rows = results.get("mechanism_comparison", [])
    if isinstance(mechanism_rows, list) and mechanism_rows:
        lines.extend(
            [
                "## Sequential vs parallel dual-access targeting (100 nm)",
                "",
                "`parallel_dual_access_v1_3` adds a direct free-to-YSA route (gated by "
                "`alpha_ysa_direct_access`) alongside the sequential free-to-OCN-to-YSA route. "
                "Neither mechanism changes cardiac output, pulmonary transit, EPR, "
                "permeability, or initial tumor-vessel supply.",
                "",
                "| Targeting | Design | Tumor AUC | Internalized (final) | Released cargo (final) | Net extravasated |",
                "| --- | --- | ---: | ---: | ---: | ---: |",
            ]
        )
        for row in mechanism_rows:
            if row.get("size") != "100nm":
                continue
            lines.append(
                f"| `{row['targeting_model']}` | `{row['design']}` | {_fmt(row['tumor_auc'])} | "
                f"{_fmt(row['internalized_fraction_final'])} | "
                f"{_fmt(row['released_cargo_fraction_final'])} | "
                f"{_fmt(row['net_extravasated_fraction'])} |"
            )
        lines.append("")

    # Feasibility domain summary.
    feasibility_rows = results.get("feasibility_domain", [])
    if isinstance(feasibility_rows, list) and feasibility_rows:
        lines.extend(["## Feasibility domain (dual vs YSA-only)", ""])
        from pathlib import Path as _Path

        summary_path = _Path(output_dir) / "tables" / "feasibility_domain_summary.csv"
        if summary_path.exists():
            import csv

            with summary_path.open(newline="", encoding="utf-8") as handle:
                summary_row = next(csv.DictReader(handle))
            lines.extend(
                [
                    f"Dual targeting (OCN-YSA) vs YSA-only, over a {summary_row['grid_dimension']}-parameter "
                    f"x {summary_row['levels_per_parameter']}-level grid ({summary_row['grid_points']} points), "
                    f"using `{summary_row['targeting_model']}`. Parameters: "
                    f"{summary_row['grid_symbols'].replace(';', ', ')}. Winning fractions "
                    "describe the sampled parameter settings.",
                    "",
                    "| Criterion (dual wins) | Sample fraction |",
                    "| --- | ---: |",
                    f"| Tumor AUC | {_fmt(summary_row['dual_wins_tumor_auc_fraction'], 4)} |",
                    f"| Final internalization | {_fmt(summary_row['dual_wins_internalized_fraction'], 4)} |",
                    f"| Final released cargo | {_fmt(summary_row['dual_wins_released_cargo_fraction'], 4)} |",
                    "",
                    "Across this grid, the dual-targeting win fractions are "
                    f"{float(summary_row['dual_wins_tumor_auc_fraction']):.1%} for tumor AUC, "
                    f"{float(summary_row['dual_wins_internalized_fraction']):.1%} for final internalization "
                    f"and {float(summary_row['dual_wins_released_cargo_fraction']):.1%} for final cargo. "
                    "OCN-mediated retention and the time required for sequential internalization "
                    "are mechanisms to examine when interpreting differences between these metrics.",
                    "",
                ]
            )
        else:
            lines.extend(["Feasibility summary was not available.", ""])

    # EPR/MPS scenarios.
    scenario_rows = results["scenario_comparison"]
    lines.extend(
        [
            "## EPR and MPS scenarios (OCN-YSA)",
            "",
            "Low, medium, and high EPR/MPS combinations are uncertainty scenarios, not "
            "measured osteosarcoma delivery rates.",
            "",
            "| EPR | MPS | Peak tumor fraction | Net extravasated | Internalized (final) |",
            "| --- | --- | ---: | ---: | ---: |",
        ]
    )
    if isinstance(scenario_rows, list):
        for row in scenario_rows:
            lines.append(
                f"| `{row['epr_label']}` | `{row['mps_label']}` | "
                f"{_fmt(row['max_tumor_fraction'])} | {_fmt(row['net_extravasated_fraction'])} | "
                f"{_fmt(row['internalized_fraction_final'])} |"
            )
    lines.append("")

    # Local sensitivity.
    sensitivity = results["local_sensitivity"]
    lines.extend(["## Multi-target local sensitivity", ""])
    if isinstance(sensitivity, list) and sensitivity:
        target_names = list(dict.fromkeys(str(row["target_metric"]) for row in sensitivity))
        for target_name in target_names:
            lines.extend(
                [
                    f"### `{target_name}`",
                    "",
                    "| Rank | Parameter | Abs normalized sensitivity |",
                    "| ---: | --- | ---: |",
                ]
            )
            target_rows = [
                row for row in sensitivity if row["target_metric"] == target_name
            ][:8]
            for idx, row in enumerate(target_rows, start=1):
                value = row["abs_normalized_sensitivity"]
                formatted = "not defined" if value is None else f"{float(value):.4g}"
                lines.append(
                    f"| {idx} | `{row['parameter_symbol']}` | {formatted} |"
                )
            lines.append("")
        lines.extend(
            [
                "Local sensitivities are central 20 % relative perturbations over the "
                "registry bounds and cover only V1.3-active uncertain parameters; "
                "physiology-fixed flows (Qco/Qliv/Qspl) are not perturbed. "
                "The default sequential dual analysis excludes alpha_ysa_direct_access "
                "and k_ysa_on because their routes are inactive for that design. "
                "The dual-vs-YSA-only grid retains k_ysa_on because it acts in YSA-only.",
                "",
            ]
        )
    else:
        lines.extend(["Sensitivity results were not available.", ""])

    # Monte Carlo.
    monte_carlo = results["monte_carlo"]
    lines.extend(["## Monte Carlo uncertainty scan", ""])
    if isinstance(monte_carlo, list) and monte_carlo:
        tumor_aucs = [float(row["tumor_auc"]) for row in monte_carlo]
        internalized = [float(row["internalized_fraction_final"]) for row in monte_carlo]
        lines.extend(
            [
                f"`{len(monte_carlo)}` samples (fixed seed) over the V1.3 registry bounds. "
                f"Tumor AUC ranges {min(tumor_aucs):.6g}-{max(tumor_aucs):.6g}; "
                f"internalized fraction ranges {min(internalized):.6g}-{max(internalized):.6g}.",
                "",
                "These are exploratory uncertainty ranges, not probabilities of clinical success.",
                "",
            ]
        )
    else:
        lines.extend(["Monte Carlo results were not available.", ""])

    # Mass-balance audit.
    audit = results["mass_balance_audit"]
    lines.extend(["## Validation audit", ""])
    if isinstance(audit, list) and audit:
        lines.extend(
            [
                "| Check | Passed | Detail |",
                "| --- | --- | --- |",
            ]
        )
        for row in audit:
            lines.append(f"| `{row['check']}` | `{row['passed']}` | {row['detail']} |")
        lines.append("")

    lines.extend(
        [
            "## Generated outputs",
            "",
            f"- `{output_root}/tables/baseline_timecourse.csv`",
            f"- `{output_root}/tables/baseline_summary.csv`",
            f"- `{output_root}/tables/baseline_detailed_metrics.csv`",
            f"- `{output_root}/tables/design_comparison.csv`",
            f"- `{output_root}/tables/size_comparison.csv`",
            f"- `{output_root}/tables/scenario_comparison.csv`",
            f"- `{output_root}/tables/local_sensitivity.csv`",
            f"- `{output_root}/tables/monte_carlo.csv`",
            f"- `{output_root}/tables/mass_balance_audit.csv`",
            f"- `{output_root}/tables/baseline_mass_distribution.csv`",
            f"- `{output_root}/tables/metric_dictionary.csv`",
            f"- `{output_root}/tables/mechanism_comparison.csv`",
            f"- `{output_root}/tables/feasibility_domain.csv`",
            f"- `{output_root}/tables/feasibility_domain_summary.csv`",
            f"- `{output_root}/figures/baseline_circulation_timecourse.svg`",
            f"- `{output_root}/figures/circulation_targeting_diagram.svg`",
            f"- `{output_root}/figures/mass_conservation.svg`",
            f"- `{output_root}/figures/mechanism_comparison.svg`",
            "",
            "## Wet-lab implications",
            "",
            "`k_entry` (extravasation), the tumor-vessel volume, and the size-dependent "
            "extravasation factor are measurement priorities because the explicit "
            "circulation makes the tumor-vessel-to-interstitium step the limiting gate. "
            "No targeting switch changes initial tumor supply, so improved internalization "
            "can only follow extravasation.",
            "",
        ]
    )

    report = Path(report_path)
    ensure_dir(report.parent)
    report.write_text("\n".join(lines), encoding="utf-8")
    return report
