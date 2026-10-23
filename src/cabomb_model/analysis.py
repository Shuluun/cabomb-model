"""Analysis workflows for CABOMB V1."""

from __future__ import annotations

from pathlib import Path
from typing import Mapping

import numpy as np

from cabomb_model.designs import DESIGNS
from cabomb_model.io import ensure_dir, write_csv
from cabomb_model.metrics import (
    max_abs_cargo_residual,
    max_abs_carrier_residual,
    minimum_state_value,
    summary_metrics,
)
from cabomb_model.model_modes import ModelModes, V1_1_MODES, V1_2_MODES
from cabomb_model.parameters import default_model_parameters, parameter_bounds
from cabomb_model.plotting import (
    write_auc_release_tradeoff,
    write_cargo_release_vs_alpha,
    write_svg_timecourse,
)
from cabomb_model.registry_artifacts import write_parameter_registry_preview
from cabomb_model.scenarios import make_scenario
from cabomb_model.solver import SimulationResult, simulate

SENSITIVITY_SYMBOLS: tuple[str, ...] = (
    "k_entry",
    "k_liver_uptake",
    "k_spleen_uptake",
    "k_blood_clear",
    "k_ocn_on",
    "k_ocn_off",
    "k_ysa_on",
    "k_ysa_off",
    "R_EphA2_eff",
    "k_int",
    "k_deg_int_acid",
    "Y_cargo",
    "alpha_ocn_cell_access",
)

MONTE_CARLO_SYMBOLS: tuple[str, ...] = (
    "k_entry",
    "k_liver_uptake",
    "k_spleen_uptake",
    "k_blood_clear",
    "k_ocn_on",
    "k_ocn_off",
    "k_ysa_on",
    "k_ysa_off",
    "R_EphA2_eff",
    "k_int",
    "k_deg_int_acid",
    "Y_cargo",
    "alpha_ocn_cell_access",
)

V12_SENSITIVITY_SYMBOLS: tuple[str, ...] = (
    "Vmax_liver",
    "Km_liver",
    "k_liver_release",
    "Vmax_spleen",
    "Km_spleen",
    "k_spleen_release",
    "k_liver_elim",
    "k_spleen_elim",
    "k_cross",
    "k_entry",
    "k_return",
    "k_ocn_on",
    "k_ocn_off",
    "R_EphA2_eff",
    "k_int",
    "k_deg_int_acid",
    "Y_cargo",
)
V12_MONTE_CARLO_SYMBOLS = V12_SENSITIVITY_SYMBOLS

_LINEAR_MPS_SYMBOLS: tuple[str, ...] = (
    "k_liver_uptake",
    "k_spleen_uptake",
    "k_blood_clear",
)
_SATURABLE_MPS_SYMBOLS: tuple[str, ...] = V12_SENSITIVITY_SYMBOLS[:8]
_ACCESSIBILITY_TARGETING_SYMBOLS: tuple[str, ...] = (
    "k_ocn_on",
    "k_ocn_off",
    "k_ysa_on",
    "k_ysa_off",
    "R_EphA2_eff",
    "k_int",
    "k_deg_int_acid",
    "Y_cargo",
    "alpha_ocn_cell_access",
)
_SEQUENTIAL_TARGETING_SYMBOLS: tuple[str, ...] = (
    "k_ocn_on",
    "k_ocn_off",
    "k_cross",
    "R_EphA2_eff",
    "k_int",
    "k_deg_int_acid",
    "Y_cargo",
)
_HYBRID_SHARED_SYMBOLS: tuple[str, ...] = ("k_entry", "k_return")

TARGET_METRICS: tuple[str, ...] = (
    "tumor_auc",
    "cumulative_internalization",
    "cumulative_cargo_release",
)

NONNEGATIVE_TOLERANCE = -1.0e-8
MASS_RESIDUAL_TOLERANCE = 1.0e-7
ALPHA_SWEEP_VALUES: tuple[float, ...] = (0.0, 0.25, 0.5, 0.75, 1.0)


def validate_simulation_result(
    result: SimulationResult,
    params: Mapping[str, float],
    *,
    context: str,
) -> None:
    """Raise if a simulation is not safe to use in saved analyses."""

    if not result.success:
        raise RuntimeError(f"{context} solver failed: {result.message}")
    minimum_value = minimum_state_value(result)
    if minimum_value < NONNEGATIVE_TOLERANCE:
        raise RuntimeError(
            f"{context} produced minimum state {minimum_value}, below "
            f"tolerance {NONNEGATIVE_TOLERANCE}."
        )
    carrier_residual = max_abs_carrier_residual(result, params)
    if carrier_residual > MASS_RESIDUAL_TOLERANCE:
        raise RuntimeError(
            f"{context} carrier residual {carrier_residual} exceeds "
            f"tolerance {MASS_RESIDUAL_TOLERANCE}."
        )
    cargo_residual = max_abs_cargo_residual(result, params)
    if cargo_residual > MASS_RESIDUAL_TOLERANCE:
        raise RuntimeError(
            f"{context} cargo residual {cargo_residual} exceeds "
            f"tolerance {MASS_RESIDUAL_TOLERANCE}."
        )


def _validated_summary(
    result: SimulationResult,
    params: Mapping[str, float],
    *,
    context: str,
) -> dict[str, float | str]:
    validate_simulation_result(result, params, context=context)
    return summary_metrics(result, params)


def _symbols_for_modes(
    modes: ModelModes,
    *,
    v11_symbols: tuple[str, ...],
    v12_symbols: tuple[str, ...],
) -> tuple[str, ...]:
    """Compose the active uncertain-parameter set for the selected modes."""

    if modes == V1_1_MODES:
        return v11_symbols
    if modes == V1_2_MODES:
        return v12_symbols

    mps_symbols = (
        _SATURABLE_MPS_SYMBOLS
        if modes.mps_model == "saturable_v1_2"
        else _LINEAR_MPS_SYMBOLS
    )
    targeting_symbols = (
        _SEQUENTIAL_TARGETING_SYMBOLS
        if modes.targeting_model == "sequential_v1_2"
        else _ACCESSIBILITY_TARGETING_SYMBOLS
    )
    return tuple(
        dict.fromkeys((*mps_symbols, *_HYBRID_SHARED_SYMBOLS, *targeting_symbols))
    )


def run_baseline(
    output_dir: str | Path = "outputs",
    *,
    modes: ModelModes = V1_1_MODES,
) -> dict[str, Path | SimulationResult]:
    """Run and save the baseline OCN-YSA medium-scenario simulation."""

    params = default_model_parameters()
    out = ensure_dir(output_dir)
    tables = ensure_dir(out / "tables")
    figures = ensure_dir(out / "figures")
    result = simulate(
        params, design="ocn-ysa", epr="medium", mps="medium", modes=modes
    )
    baseline_summary = _validated_summary(result, params, context="baseline")
    trajectory_rows = result.as_rows()
    for row in trajectory_rows:
        row["mps_model"] = result.modes.mps_model
        row["targeting_model"] = result.modes.targeting_model
    trajectory_path = write_csv(tables / "baseline_timecourse.csv", trajectory_rows)
    summary_path = write_csv(tables / "baseline_summary.csv", [baseline_summary])
    figure_path = write_svg_timecourse(
        result,
        figures / "baseline_timecourse.svg",
        title="Baseline OCN-YSA carrier and cargo time courses",
    )
    return {
        "result": result,
        "trajectory_path": trajectory_path,
        "summary_path": summary_path,
        "figure_path": figure_path,
    }


def run_design_comparison(
    output_dir: str | Path = "outputs",
    *,
    modes: ModelModes = V1_1_MODES,
) -> list[dict[str, object]]:
    """Run all four carrier designs under the medium exploratory scenario."""

    params = default_model_parameters()
    rows: list[dict[str, object]] = []
    out = ensure_dir(output_dir)
    for design_name in DESIGNS:
        result = simulate(
            params, design=design_name, epr="medium", mps="medium", modes=modes
        )
        rows.append(
            _validated_summary(result, params, context=f"design comparison {design_name}")
        )
    write_csv(out / "tables" / "design_comparison.csv", rows)
    return rows


def run_epr_mps_scenarios(
    output_dir: str | Path = "outputs",
    *,
    modes: ModelModes = V1_1_MODES,
) -> list[dict[str, object]]:
    """Run low, medium, and high EPR/MPS scenario combinations."""

    params = default_model_parameters()
    rows: list[dict[str, object]] = []
    out = ensure_dir(output_dir)
    for epr in ("low", "medium", "high"):
        for mps in ("low", "medium", "high"):
            for design_name in DESIGNS:
                result = simulate(
                    params, design=design_name, epr=epr, mps=mps, modes=modes
                )
                row = _validated_summary(
                    result,
                    params,
                    context=f"scenario epr={epr} mps={mps} design={design_name}",
                )
                row["epr_label"] = epr
                row["mps_label"] = mps
                rows.append(row)
    write_csv(out / "tables" / "scenario_comparison.csv", rows)
    return rows


def run_ocn_accessibility_sweep(
    output_dir: str | Path = "outputs",
    *,
    modes: ModelModes = V1_1_MODES,
) -> list[dict[str, object]]:
    """Run the four designs across five exploratory OCN-accessibility levels."""

    if modes.targeting_model != "accessibility_v1_1":
        raise ValueError(
            "OCN accessibility sweep requires targeting_model='accessibility_v1_1'."
        )

    params = default_model_parameters()
    out = ensure_dir(output_dir)
    rows: list[dict[str, object]] = []
    for alpha in ALPHA_SWEEP_VALUES:
        run_params = dict(params, alpha_ocn_cell_access=alpha)
        for design_name in DESIGNS:
            result = simulate(
                run_params,
                design=design_name,
                epr="medium",
                mps="medium",
                modes=modes,
            )
            row = _validated_summary(
                result,
                run_params,
                context=(
                    f"OCN accessibility alpha={alpha} design={design_name}"
                ),
            )
            row["alpha_ocn_cell_access"] = alpha
            rows.append(row)
    write_csv(out / "tables" / "ocn_accessibility_sweep.csv", rows)
    write_cargo_release_vs_alpha(
        rows, out / "figures" / "cumulative_cargo_release_vs_alpha.svg"
    )
    write_auc_release_tradeoff(
        rows, out / "figures" / "tumor_auc_vs_cumulative_cargo_release.svg"
    )
    return rows


def _target_metric(
    params: Mapping[str, float], target_metric: str, *, modes: ModelModes
) -> float:
    if target_metric not in TARGET_METRICS:
        allowed = ", ".join(TARGET_METRICS)
        raise ValueError(f"Unknown target_metric {target_metric!r}; choose {allowed}.")
    result = simulate(
        params, design="ocn-ysa", epr="medium", mps="medium", modes=modes
    )
    validate_simulation_result(result, params, context=f"sensitivity {target_metric}")
    return float(summary_metrics(result, params)[target_metric])


def run_local_sensitivity(
    output_dir: str | Path = "outputs",
    target_metrics: tuple[str, ...] = TARGET_METRICS,
    *,
    modes: ModelModes = V1_1_MODES,
) -> list[dict[str, object]]:
    """Run one-at-a-time local sensitivity for the requested target metrics."""

    params = default_model_parameters()
    bounds = parameter_bounds()
    rows: list[dict[str, object]] = []
    for target_metric in target_metrics:
        base_metric = _target_metric(params, target_metric, modes=modes)
        target_rows: list[dict[str, object]] = []
        sensitivity_symbols = _symbols_for_modes(
            modes,
            v11_symbols=SENSITIVITY_SYMBOLS,
            v12_symbols=V12_SENSITIVITY_SYMBOLS,
        )
        for symbol in sensitivity_symbols:
            lower, upper = bounds[symbol]
            baseline = params[symbol]
            if symbol == "alpha_ocn_cell_access":
                low_value = 0.0
                high_value = 0.1
                method = "one_sided_boundary"
            else:
                low_value = max(lower, baseline * 0.8)
                high_value = min(upper, baseline * 1.2)
                method = "central_relative_20pct"

            low_params = dict(params)
            high_params = dict(params)
            low_params[symbol] = low_value
            high_params[symbol] = high_value
            low_metric = _target_metric(low_params, target_metric, modes=modes)
            high_metric = _target_metric(high_params, target_metric, modes=modes)
            denominator = high_value - low_value
            slope = (high_metric - low_metric) / denominator if denominator else 0.0
            if method == "one_sided_boundary":
                range_fraction = denominator / (upper - lower) if upper != lower else 0.0
                normalized = (
                    ((high_metric - low_metric) / base_metric) / range_fraction
                    if base_metric and range_fraction
                    else None
                )
            else:
                normalized = slope * baseline / base_metric if base_metric else None
            abs_normalized = abs(normalized) if normalized is not None else None
            target_rows.append(
                {
                    "mps_model": modes.mps_model,
                    "targeting_model": modes.targeting_model,
                    "target_metric": target_metric,
                    "parameter_symbol": symbol,
                    "sensitivity_method": method,
                    "baseline": baseline,
                    "low_value": low_value,
                    "high_value": high_value,
                    "baseline_metric": base_metric,
                    "low_metric": low_metric,
                    "high_metric": high_metric,
                    "finite_difference_slope": slope,
                    "normalized_sensitivity": normalized,
                    "abs_normalized_sensitivity": abs_normalized,
                }
            )
        target_rows.sort(
            key=lambda row: (
                row["abs_normalized_sensitivity"] is None,
                -float(row["abs_normalized_sensitivity"])
                if row["abs_normalized_sensitivity"] is not None
                else 0.0,
            )
        )
        rows.extend(target_rows)
    write_csv(Path(output_dir) / "tables" / "local_sensitivity.csv", rows)
    return rows


def run_monte_carlo(
    output_dir: str | Path = "outputs",
    *,
    n_samples: int = 128,
    seed: int = 20260709,
    modes: ModelModes = V1_1_MODES,
) -> list[dict[str, object]]:
    """Run a reproducible uncertainty scan using registry bounds."""

    params = default_model_parameters()
    bounds = parameter_bounds()
    rng = np.random.default_rng(seed)
    rows: list[dict[str, object]] = []
    monte_carlo_symbols = _symbols_for_modes(
        modes,
        v11_symbols=MONTE_CARLO_SYMBOLS,
        v12_symbols=V12_MONTE_CARLO_SYMBOLS,
    )
    for sample_idx in range(n_samples):
        sample_params = dict(params)
        for symbol in monte_carlo_symbols:
            lower, upper = bounds[symbol]
            sample_params[symbol] = float(rng.uniform(lower, upper))
        sampled_columns = {
            symbol: sample_params[symbol] for symbol in monte_carlo_symbols
        }
        for design_name in DESIGNS:
            result = simulate(
                sample_params,
                design=design_name,
                epr="medium",
                mps="medium",
                modes=modes,
            )
            row = _validated_summary(
                result,
                sample_params,
                context=f"monte carlo sample {sample_idx} design {design_name}",
            )
            row.update(
                {
                    "sample_id": sample_idx,
                    "seed": seed,
                    **sampled_columns,
                }
            )
            rows.append(row)
    summary = summarize_paired_design_advantage(rows)
    write_csv(Path(output_dir) / "tables" / "monte_carlo.csv", rows)
    write_csv(
        Path(output_dir) / "tables" / "monte_carlo_design_advantage_summary.csv",
        [summary],
    )
    return rows


def summarize_paired_design_advantage(
    rows: list[dict[str, object]],
) -> dict[str, object]:
    """Summarize strict OCN-YSA wins over YSA within paired samples."""

    if not rows:
        raise ValueError("Paired Monte Carlo summary requires at least one sample.")

    grouped: dict[int, dict[str, dict[str, object]]] = {}
    mode_pairs: set[tuple[str, str]] = set()
    for row in rows:
        mps_model = str(row.get("mps_model", ""))
        targeting_model = str(row.get("targeting_model", ""))
        if not mps_model or not targeting_model:
            raise ValueError(
                "Paired Monte Carlo rows must record one consistent mode pair."
            )
        mode_pairs.add((mps_model, targeting_model))
        sample_id = int(row["sample_id"])
        design_name = str(row["design"])
        sample = grouped.setdefault(sample_id, {})
        if design_name in sample:
            raise ValueError(f"Duplicate design {design_name} for sample {sample_id}.")
        sample[design_name] = row

    if len(mode_pairs) != 1:
        raise ValueError("Paired Monte Carlo rows must record one consistent mode pair.")
    mps_model, targeting_model = next(iter(mode_pairs))

    cargo_wins = 0
    auc_wins = 0
    both_wins = 0
    for sample_id, designs in grouped.items():
        if set(designs) != set(DESIGNS):
            raise ValueError(f"Sample {sample_id} does not contain all four designs.")
        combo = designs["ocn-ysa"]
        ysa = designs["ysa"]
        combo_cargo = float(combo["cumulative_cargo_release"])
        ysa_cargo = float(ysa["cumulative_cargo_release"])
        combo_auc = float(combo["tumor_auc"])
        ysa_auc = float(ysa["tumor_auc"])
        if not all(np.isfinite((combo_cargo, ysa_cargo, combo_auc, ysa_auc))):
            raise ValueError(
                f"Sample {sample_id} must contain finite comparison metrics."
            )
        cargo_win = combo_cargo > ysa_cargo
        auc_win = combo_auc > ysa_auc
        cargo_wins += int(cargo_win)
        auc_wins += int(auc_win)
        both_wins += int(cargo_win and auc_win)

    sample_count = len(grouped)
    return {
        "mps_model": mps_model,
        "targeting_model": targeting_model,
        "sample_count": sample_count,
        "cargo_release_win_count": cargo_wins,
        "cargo_release_win_fraction": cargo_wins / sample_count,
        "tumor_auc_win_count": auc_wins,
        "tumor_auc_win_fraction": auc_wins / sample_count,
        "both_win_count": both_wins,
        "both_win_fraction": both_wins / sample_count,
        "interpretation": (
            "These are paired-sample fractions under exploratory parameter distributions, "
            "not probabilities of clinical success."
        ),
    }


def run_all_analyses(
    output_dir: str | Path = "outputs",
    *,
    modes: ModelModes = V1_1_MODES,
) -> dict[str, object]:
    """Run all currently implemented analyses and write outputs."""

    baseline = run_baseline(output_dir, modes=modes)
    design_rows = run_design_comparison(output_dir, modes=modes)
    scenario_rows = run_epr_mps_scenarios(output_dir, modes=modes)
    sweep_rows = (
        run_ocn_accessibility_sweep(output_dir, modes=modes)
        if modes.targeting_model == "accessibility_v1_1"
        else []
    )
    sensitivity_rows = run_local_sensitivity(output_dir, modes=modes)
    monte_carlo_rows = run_monte_carlo(output_dir, modes=modes)
    monte_carlo_summary = summarize_paired_design_advantage(monte_carlo_rows)
    preview_path = write_parameter_registry_preview(output_dir)
    return {
        "baseline": baseline,
        "design_comparison": design_rows,
        "scenario_comparison": scenario_rows,
        "ocn_accessibility_sweep": sweep_rows,
        "local_sensitivity": sensitivity_rows,
        "monte_carlo": monte_carlo_rows,
        "monte_carlo_design_advantage_summary": monte_carlo_summary,
        "parameter_registry_preview_path": preview_path,
    }
