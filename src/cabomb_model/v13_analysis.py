"""Analysis workflows for the optional V1.3 explicit-circulation mode.

This module mirrors ``analysis.py`` for the V1.1/V1.2 path but drives the V1.3
eighteen-state solver and its own metrics. It is separate so the frozen V1.1/V1.2
pipeline and its outputs are never modified.
"""

from __future__ import annotations

from pathlib import Path
from typing import Mapping
import hashlib
import json
import platform
import subprocess
from importlib.metadata import version

import numpy as np

from cabomb_model.designs import DESIGNS
from cabomb_model.io import ensure_dir, write_csv
from cabomb_model.model_modes import ModelModes, V1_3_MODES
from cabomb_model.parameters import load_registry
from cabomb_model.v13_circulation import (
    SIZE_KEYS,
    V13Result,
    v13_carrier_total,
    v13_cargo_total,
    v13_circulation_metrics,
    v13_cumulative_fluxes,
    v13_binding_internalization_metrics,
    v13_detailed_metrics,
    v13_first_pass_metrics,
    v13_max_abs_cargo_residual,
    v13_max_abs_carrier_residual,
    v13_minimum_state_value,
    v13_simulate,
    v13_summary_metrics,
    v13_tumor_associated,
)
from cabomb_model.v13_plotting import (
    write_v13_circulation_timecourse,
    write_v13_mass_conservation,
    write_v13_mechanism_comparison,
    write_v13_mechanism_diagram,
)
from cabomb_model.v13_states import V13_CARRIER_STATES, V13_VASCULAR_STATES

# Uncertainty candidates for V1.3. The default local/Monte Carlo workflows use
# the sequential dual design and filter structurally inactive candidates below.
# physiology-fixed flows (Qco/Qliv/Qspl) are excluded because their registry
# bounds are single-point constants. alpha_ysa_direct_access is exclusive to
# the parallel targeting mode and is evaluated in the mechanism comparison.
V13_BOUND_SYMBOLS: tuple[str, ...] = (
    "k_entry",
    "k_return",
    "k_int",
    "k_int_ns",
    "k_ocn_on",
    "k_ocn_off",
    "k_cross",
    "k_ysa_on",
    "k_ysa_off",
    "alpha_ysa_direct_access",
    "R_EphA2_eff",
    "k_deg_int_acid",
    "k_blood_clear",
    "k_tumor_deg",
    "Vmax_liver",
    "Km_liver",
    "k_liver_release",
    "Vmax_spleen",
    "Km_spleen",
    "k_spleen_release",
    "Qtum",
    "V_tumor_vascular",
    "Y_cargo",
)
V13_SENSITIVITY_SYMBOLS: tuple[str, ...] = tuple(
    symbol for symbol in V13_BOUND_SYMBOLS
    if symbol not in {"alpha_ysa_direct_access", "k_ysa_on"}
)
V13_MONTE_CARLO_SYMBOLS: tuple[str, ...] = V13_SENSITIVITY_SYMBOLS

V13_TARGET_METRICS: tuple[str, ...] = (
    "tumor_auc",
    "internalized_fraction_final",
    "released_cargo_fraction_final",
)

NONNEGATIVE_TOLERANCE = -1.0e-8
MASS_RESIDUAL_TOLERANCE = 1.0e-7

# States whose time courses are required by the V1.3 specification.
V13_VASCULAR_TIMECOURSE_STATES: tuple[str, ...] = (
    "A_central_venous",
    "A_right_heart",
    "A_lung_vascular",
    "A_left_heart",
    "A_systemic_arterial",
    "A_tumor_vascular",
)
V13_TUMOR_TIMECOURSE_STATES: tuple[str, ...] = (
    "A_tumor_free",
    "A_ocn_bound",
    "A_ysa_bound",
    "A_internalized",
    "A_cargo",
)


def v13_parameter_bounds(path: str | Path | None = None) -> dict[str, tuple[float, float]]:
    """Return registry lower/upper bounds for the V1.3 uncertain symbols."""

    records = load_registry(path)
    by_symbol = {row["parameter_symbol"]: row for row in records}
    bounds: dict[str, tuple[float, float]] = {}
    for symbol in V13_BOUND_SYMBOLS:
        if symbol not in by_symbol:
            raise ValueError(f"V1.3 uncertain parameter {symbol} is absent from registry.")
        bounds[symbol] = (
            float(by_symbol[symbol]["lower_bound"]),
            float(by_symbol[symbol]["upper_bound"]),
        )
    return bounds


def validate_v13_result(result: V13Result, params: Mapping[str, float], *, context: str) -> None:
    """Raise if a V1.3 simulation is not safe to use in saved analyses."""

    if not result.success:
        raise RuntimeError(f"{context} solver failed: {result.message}")
    minimum_value = v13_minimum_state_value(result)
    if minimum_value < NONNEGATIVE_TOLERANCE:
        raise RuntimeError(
            f"{context} produced minimum state {minimum_value}, below "
            f"tolerance {NONNEGATIVE_TOLERANCE}."
        )
    carrier_residual = v13_max_abs_carrier_residual(result, params)
    if carrier_residual > MASS_RESIDUAL_TOLERANCE:
        raise RuntimeError(
            f"{context} carrier residual {carrier_residual} exceeds "
            f"tolerance {MASS_RESIDUAL_TOLERANCE}."
        )
    cargo_residual = v13_max_abs_cargo_residual(result, params)
    if cargo_residual > MASS_RESIDUAL_TOLERANCE:
        raise RuntimeError(
            f"{context} cargo residual {cargo_residual} exceeds "
            f"tolerance {MASS_RESIDUAL_TOLERANCE}."
        )


def _validated_v13_summary(
    result: V13Result, params: Mapping[str, float], *, context: str
) -> dict[str, float | str]:
    validate_v13_result(result, params, context=context)
    return v13_summary_metrics(result, params)


def run_v13_baseline(
    output_dir: str | Path = "outputs/v1_3",
    *,
    size: str = "100nm",
    design: str = "ocn-ysa",
    epr: str = "medium",
    mps: str = "medium",
    modes: ModelModes = V1_3_MODES,
) -> dict[str, object]:
    """Run and save the V1.3 baseline explicit-circulation simulation."""

    from cabomb_model.v13_circulation import v13_default_parameters

    params = v13_default_parameters()
    out = ensure_dir(output_dir)
    tables = ensure_dir(out / "tables")
    figures = ensure_dir(out / "figures")
    result = v13_simulate(params, design=design, epr=epr, mps=mps, size=size, modes=modes)
    baseline_summary = _validated_v13_summary(result, params, context="V1.3 baseline")
    trajectory_rows = result.as_rows()
    trajectory_path = write_csv(tables / "baseline_timecourse.csv", trajectory_rows)
    summary_path = write_csv(tables / "baseline_summary.csv", [baseline_summary])
    figure_path = write_v13_circulation_timecourse(
        result, figures / "baseline_circulation_timecourse.svg"
    )
    diagram_path = write_v13_mechanism_diagram(figures / "circulation_targeting_diagram.svg")
    mass_path = write_v13_mass_conservation(
        result, params, figures / "mass_conservation.svg"
    )
    return {
        "result": result,
        "params": params,
        "trajectory_path": trajectory_path,
        "summary_path": summary_path,
        "figure_path": figure_path,
        "diagram_path": diagram_path,
        "mass_path": mass_path,
    }


def run_v13_design_comparison(
    output_dir: str | Path = "outputs/v1_3",
    *,
    size: str = "100nm",
    modes: ModelModes = V1_3_MODES,
) -> list[dict[str, object]]:
    """Save a scoped four-design comparison, trajectories and traceable report.

    All designs reuse one parameter dictionary, size, mode and saved time grid.
    This entry point does not run the feasibility scan or Monte Carlo workflow.
    """

    from cabomb_model.v13_circulation import v13_default_parameters

    params = v13_default_parameters()
    out = ensure_dir(output_dir)
    results: dict[str, V13Result] = {}
    fluxes: dict[str, dict[str, np.ndarray]] = {}
    rows: list[dict[str, object]] = []
    for design_name in DESIGNS:
        result = v13_simulate(params, design=design_name, size=size, modes=modes)
        row = _validated_v13_summary(
            result, params, context=f"V1.3 design comparison {design_name}"
        )
        row.update(v13_binding_internalization_metrics(result, params))
        row["time_final_h"] = float(result.t[-1])
        row["initial_dose_normalized"] = params["D0"]
        row["final_tumor_stock_fraction"] = float(v13_tumor_associated(result)[-1] / params["D0"]) if params["D0"] else 0.0
        rows.append(row)
        results[design_name] = result
        fluxes[design_name] = v13_cumulative_fluxes(result, params)
        write_csv(out / "tables" / f"{design_name}_timecourse.csv", result.as_rows())
        write_csv(out / "tables" / f"{design_name}_cumulative_fluxes.csv", [
            {"time_h": float(t), **{k: float(v[i]) for k, v in fluxes[design_name].items()}}
            for i, t in enumerate(result.t)
        ])
    write_csv(out / "tables" / "design_comparison.csv", rows)
    write_csv(out / "tables" / "effective_parameters.csv", [
        {"parameter_symbol": key, "value": value} for key, value in sorted(params.items())
    ])
    write_csv(out / "tables" / "parameter_registry_snapshot.csv", load_registry())
    from cabomb_model.v13_plotting import write_v13_ablation_timecourses
    from cabomb_model.v13_reporting import write_v13_ablation_summary

    write_v13_ablation_timecourses(results, fluxes, params, out / "figures" / "four_design_timecourses.svg")
    write_v13_ablation_summary(rows, out / "summary_zh.md")
    from cabomb_model.resource_paths import source_root
    from cabomb_model.parameters import default_registry_path
    checkout = source_root()
    root = Path(__file__).resolve().parent
    sources = sorted(root.glob("*.py"))
    registry = default_registry_path()

    def git_read(*args: str) -> str:
        if checkout is None or not (checkout / ".git").exists():
            return "unavailable"
        try:
            return subprocess.check_output(["git", *args], cwd=checkout, text=True, stderr=subprocess.DEVNULL).strip()
        except (OSError, subprocess.CalledProcessError):
            return "unavailable"

    first = next(iter(results.values()))
    manifest = {
        "schema_version": 1,
        "entry_point": "cabomb_model.v13_analysis.run_v13_design_comparison",
        "git_head_at_run": git_read("rev-parse", "HEAD"),
        "git_status_at_run": git_read("status", "--porcelain=v1", "-uall"),
        "source_sha256": {
            **{str(p.relative_to(root)).replace("\\", "/"): hashlib.sha256(p.read_bytes()).hexdigest() for p in sources},
            "parameter_registry.csv": hashlib.sha256(registry.read_bytes()).hexdigest(),
        },
        "parameters": params,
        "design_switches": {name: {"i_ocn": d.i_ocn, "i_ysa": d.i_ysa} for name, d in DESIGNS.items()},
        "size": size, "design_baseline_not_measured_size": True,
        "modes": {"mps": modes.mps_model, "targeting": modes.targeting_model, "circulation": modes.circulation_model},
        "scenario": vars(first.scenario), "time_grid_h": first.t.tolist(),
        "solver": {"method": first.solver_method, "rtol": first.rtol, "atol": first.atol},
        "versions": {"python": platform.python_version(), **{name: version(name) for name in ("numpy", "scipy", "matplotlib")}},
        "interpretation": "Exploratory normalized IV carrier model; no per-cell absolute uptake or free cytosolic calcium prediction.",
        "output_sha256": {str(p.relative_to(out)).replace("\\", "/"): hashlib.sha256(p.read_bytes()).hexdigest() for p in sorted(out.rglob("*")) if p.is_file() and p.name != "run_manifest.json"},
    }
    (out / "run_manifest.json").write_text(json.dumps(manifest, ensure_ascii=False, indent=2, allow_nan=False) + "\n", encoding="utf-8")
    return rows


def run_v13_size_comparison(
    output_dir: str | Path = "outputs/v1_3",
) -> list[dict[str, object]]:
    """Run the OCN-YSA design across all six size scenarios."""

    from cabomb_model.v13_circulation import v13_default_parameters

    params = v13_default_parameters()
    rows: list[dict[str, object]] = []
    for size in SIZE_KEYS:
        result = v13_simulate(params, design="ocn-ysa", size=size)
        rows.append(
            _validated_v13_summary(result, params, context=f"V1.3 size comparison {size}")
        )
    write_csv(Path(output_dir) / "tables" / "size_comparison.csv", rows)
    return rows


def run_v13_epr_mps_scenarios(
    output_dir: str | Path = "outputs/v1_3",
    *,
    size: str = "100nm",
) -> list[dict[str, object]]:
    """Run low, medium, and high EPR/MPS combinations for the OCN-YSA design."""

    from cabomb_model.v13_circulation import v13_default_parameters

    params = v13_default_parameters()
    rows: list[dict[str, object]] = []
    for epr in ("low", "medium", "high"):
        for mps in ("low", "medium", "high"):
            result = v13_simulate(params, design="ocn-ysa", epr=epr, mps=mps, size=size)
            row = _validated_v13_summary(
                result, params, context=f"V1.3 scenario epr={epr} mps={mps}"
            )
            row["epr_label"] = epr
            row["mps_label"] = mps
            rows.append(row)
    write_csv(Path(output_dir) / "tables" / "scenario_comparison.csv", rows)
    return rows


def _v13_target_metric(
    params: Mapping[str, float], target_metric: str, *, size: str = "100nm"
) -> float:
    if target_metric not in V13_TARGET_METRICS:
        allowed = ", ".join(V13_TARGET_METRICS)
        raise ValueError(f"Unknown V1.3 target_metric {target_metric!r}; choose {allowed}.")
    result = v13_simulate(params, design="ocn-ysa", size=size, modes=V1_3_MODES)
    validate_v13_result(result, params, context=f"V1.3 sensitivity {target_metric}")
    return float(v13_summary_metrics(result, params)[target_metric])


def run_v13_local_sensitivity(
    output_dir: str | Path = "outputs/v1_3",
    target_metrics: tuple[str, ...] = V13_TARGET_METRICS,
    *,
    size: str = "100nm",
) -> list[dict[str, object]]:
    """Run one-at-a-time local sensitivity for the requested target metrics."""

    from cabomb_model.v13_circulation import v13_default_parameters

    params = v13_default_parameters()
    bounds = v13_parameter_bounds()
    rows: list[dict[str, object]] = []
    for target_metric in target_metrics:
        base_metric = _v13_target_metric(params, target_metric, size=size)
        target_rows: list[dict[str, object]] = []
        for symbol in V13_SENSITIVITY_SYMBOLS:
            lower, upper = bounds[symbol]
            baseline = params[symbol]
            low_value = max(lower, baseline * 0.8)
            high_value = min(upper, baseline * 1.2)
            low_params = dict(params)
            high_params = dict(params)
            low_params[symbol] = low_value
            high_params[symbol] = high_value
            low_metric = _v13_target_metric(low_params, target_metric, size=size)
            high_metric = _v13_target_metric(high_params, target_metric, size=size)
            denominator = high_value - low_value
            slope = (high_metric - low_metric) / denominator if denominator else 0.0
            normalized = slope * baseline / base_metric if base_metric else None
            abs_normalized = abs(normalized) if normalized is not None else None
            target_rows.append(
                {
                    "circulation_model": "explicit_v1_3",
                    "mps_model": V1_3_MODES.mps_model,
                    "targeting_model": V1_3_MODES.targeting_model,
                    "target_metric": target_metric,
                    "parameter_symbol": symbol,
                    "sensitivity_method": "central_relative_20pct",
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


def run_v13_monte_carlo(
    output_dir: str | Path = "outputs/v1_3",
    *,
    n_samples: int = 96,
    seed: int = 20260709,
    size: str = "100nm",
) -> list[dict[str, object]]:
    """Run a reproducible uncertainty scan over V1.3 registry bounds."""

    from cabomb_model.v13_circulation import v13_default_parameters

    params = v13_default_parameters()
    bounds = v13_parameter_bounds()
    rng = np.random.default_rng(seed)
    rows: list[dict[str, object]] = []
    for sample_idx in range(n_samples):
        sample_params = dict(params)
        for symbol in V13_MONTE_CARLO_SYMBOLS:
            lower, upper = bounds[symbol]
            sample_params[symbol] = float(rng.uniform(lower, upper))
        sampled_columns = {
            symbol: sample_params[symbol] for symbol in V13_MONTE_CARLO_SYMBOLS
        }
        result = v13_simulate(sample_params, design="ocn-ysa", size=size, modes=V1_3_MODES)
        row = _validated_v13_summary(
            result, sample_params, context=f"V1.3 monte carlo sample {sample_idx}"
        )
        row.update({"sample_id": sample_idx, "seed": seed, **sampled_columns})
        rows.append(row)
    write_csv(Path(output_dir) / "tables" / "monte_carlo.csv", rows)
    return rows


def v13_mass_distribution(
    result: V13Result, params: Mapping[str, float]
) -> dict[str, float]:
    """Return the final carrier mass distribution by compartment and loss route."""

    dose = params["D0"]
    cumulative = v13_cumulative_fluxes(result, params)
    final = result.y[:, -1]
    vascular = sum(final[result.state_names.index(name)] for name in V13_VASCULAR_STATES)
    tumor = v13_tumor_associated(result)[-1]
    internalized = float(result.state("A_internalized")[-1])

    def frac(value: float) -> float:
        return float(value / dose) if dose else 0.0

    return {
        "vascular_fraction_final": frac(vascular),
        "liver_fraction_final": frac(result.state("A_liver")[-1]),
        "spleen_fraction_final": frac(result.state("A_spleen")[-1]),
        "tumor_associated_fraction_final": frac(tumor),
        "internalized_fraction_final": frac(internalized),
        "elim_deg_fraction_final": frac(result.state("A_elim_deg")[-1]),
        "cumulative_lung_ret_fraction": frac(float(cumulative["lung_ret"][-1])),
        "cumulative_blood_clear_fraction": frac(float(cumulative["blood_clear"][-1])),
        "cumulative_liver_elim_fraction": frac(float(cumulative["liver_elim"][-1])),
        "cumulative_spleen_elim_fraction": frac(float(cumulative["spleen_elim"][-1])),
        "cumulative_tumor_deg_fraction": frac(float(cumulative["tumor_deg"][-1])),
        "cumulative_ocn_loss_fraction": frac(float(cumulative["ocn_loss"][-1])),
        "cumulative_deg_int_fraction": frac(float(cumulative["deg_int"][-1])),
    }


def run_v13_mass_balance_audit(
    output_dir: str | Path = "outputs/v1_3",
) -> list[dict[str, object]]:
    """Run the V1.3 validation checks and write a single audit table."""

    from cabomb_model.v13_circulation import v13_default_parameters

    params = v13_default_parameters()
    rows: list[dict[str, object]] = []
    out = ensure_dir(output_dir)

    # 1. Zero-dose simulation: every state stays at zero.
    zero_params = dict(params, D0=0.0)
    zero = v13_simulate(zero_params, design="ocn-ysa")
    rows.append(
        {
            "check": "zero_dose_zeros_all_states",
            "passed": str(bool(zero.success) and v13_minimum_state_value(zero) == 0.0
                         and np.max(np.abs(zero.y)) == 0.0),
            "carrier_residual": v13_max_abs_carrier_residual(zero, zero_params),
            "cargo_residual": v13_max_abs_cargo_residual(zero, zero_params),
            "minimum_state_value": v13_minimum_state_value(zero),
            "detail": "IV bolus of zero must leave all 18 states at zero.",
        }
    )

    # 2. Baseline mass conservation and non-negativity.
    baseline = v13_simulate(params, design="ocn-ysa")
    rows.append(
        {
            "check": "carrier_mass_conservation",
            "passed": str(
                v13_max_abs_carrier_residual(baseline, params) < MASS_RESIDUAL_TOLERANCE
            ),
            "carrier_residual": v13_max_abs_carrier_residual(baseline, params),
            "cargo_residual": v13_max_abs_cargo_residual(baseline, params),
            "minimum_state_value": v13_minimum_state_value(baseline),
            "detail": "Carrier ledger must close to D0 within tolerance.",
        }
    )
    rows.append(
        {
            "check": "cargo_mass_conservation",
            "passed": str(
                v13_max_abs_cargo_residual(baseline, params) < MASS_RESIDUAL_TOLERANCE
            ),
            "carrier_residual": v13_max_abs_carrier_residual(baseline, params),
            "cargo_residual": v13_max_abs_cargo_residual(baseline, params),
            "minimum_state_value": v13_minimum_state_value(baseline),
            "detail": "Cargo ledger must close to Y_cargo * D0 within tolerance.",
        }
    )
    rows.append(
        {
            "check": "nonnegativity",
            "passed": str(v13_minimum_state_value(baseline) >= NONNEGATIVE_TOLERANCE),
            "carrier_residual": v13_max_abs_carrier_residual(baseline, params),
            "cargo_residual": v13_max_abs_cargo_residual(baseline, params),
            "minimum_state_value": v13_minimum_state_value(baseline),
            "detail": "No state may be negative beyond numerical tolerance.",
        }
    )
    rows.append(
        {
            "check": "solver_success",
            "passed": str(baseline.success),
            "carrier_residual": v13_max_abs_carrier_residual(baseline, params),
            "cargo_residual": v13_max_abs_cargo_residual(baseline, params),
            "minimum_state_value": v13_minimum_state_value(baseline),
            "detail": baseline.message,
        }
    )

    # 3. Flow branch sum equals cardiac output (structural identity).
    qco = params["Qco"] * 60.0
    qrest = qco - params["Qliv"] * 60.0 - params["Qspl"] * 60.0 - params["Qtum"] * 60.0
    rows.append(
        {
            "check": "flow_branch_sum_equals_cardiac_output",
            "passed": str(qrest >= 0.0),
            "carrier_residual": 0.0,
            "cargo_residual": 0.0,
            "minimum_state_value": 0.0,
            "detail": f"Qrest = Qco - Qliv - Qspl - Qtum = {qrest:.6g} L/h (>= 0).",
        }
    )

    # 4. Deterministic reproducibility (same inputs, same trajectory).
    repeat = v13_simulate(params, design="ocn-ysa")
    identical = bool(np.allclose(baseline.y, repeat.y, rtol=0.0, atol=0.0))
    rows.append(
        {
            "check": "deterministic_reproducibility",
            "passed": str(identical),
            "carrier_residual": 0.0,
            "cargo_residual": 0.0,
            "minimum_state_value": 0.0,
            "detail": "Identical parameters and seed must reproduce identical trajectories.",
        }
    )

    # 5. OCN/YSA targeting cannot change arterial supply or initial tumor-vessel
    # arrival. Two layers: (a) at the flux-function level the circulation and
    # tumor-entry fluxes are exactly identical for a fixed state vector across
    # all four designs; (b) on the trajectory, the early (first-pass) arterial
    # and tumor-vessel states are identical before extravasation can feed back.
    from cabomb_model.v13_circulation import v13_fluxes

    probe_state = baseline.y[:, 0].copy()
    probe_state[baseline.state_names.index("A_tumor_free")] = 0.1
    probe_state[baseline.state_names.index("A_ocn_bound")] = 0.05
    probe_state[baseline.state_names.index("A_ysa_bound")] = 0.05
    circulation_keys = (
        "cv_to_rh",
        "rh_to_lung",
        "lung_pass",
        "lung_ret",
        "lh_to_sa",
        "sa_to_liver",
        "sa_to_spleen",
        "sa_to_tumor",
        "sa_to_rest",
        "liver_venous",
        "spleen_venous",
        "tumor_venous",
        "rest_venous",
        "entry",
    )
    base_fluxes = v13_fluxes(probe_state, params, "unmodified", baseline.scenario, baseline.size)
    flux_invariant = all(
        all(
            v13_fluxes(probe_state, params, design_name, baseline.scenario, baseline.size)[key]
            == base_fluxes[key]
            for key in circulation_keys
        )
        for design_name in ("ocn", "ysa", "ocn-ysa")
    )

    unmodified = v13_simulate(params, design="unmodified")
    early_mask = unmodified.t <= 0.05
    arterial_identical = bool(
        np.allclose(
            unmodified.state("A_systemic_arterial")[early_mask],
            baseline.state("A_systemic_arterial")[early_mask],
            rtol=1e-6,
            atol=1e-9,
        )
    )
    tumor_vessel_identical = bool(
        np.allclose(
            unmodified.state("A_tumor_vascular")[early_mask],
            baseline.state("A_tumor_vascular")[early_mask],
            rtol=1e-6,
            atol=1e-9,
        )
    )
    rows.append(
        {
            "check": "targeting_does_not_change_supply",
            "passed": str(flux_invariant and arterial_identical and tumor_vessel_identical),
            "carrier_residual": 0.0,
            "cargo_residual": 0.0,
            "minimum_state_value": 0.0,
            "detail": (
                "Circulation and entry fluxes are design-invariant for a fixed state "
                "vector, and the first-pass arterial/tumor-vessel states are identical "
                "across designs before extravasation feedback."
            ),
        }
    )

    # 6. Extravasation off blocks all tumor interstitial and targeting states.
    no_entry_params = dict(params, k_entry=0.0)
    no_entry = v13_simulate(no_entry_params, design="ocn-ysa")
    tumor_blocked = bool(
        np.max(np.abs(no_entry.state("A_tumor_free"))) == 0.0
        and np.max(np.abs(no_entry.state("A_ocn_bound"))) == 0.0
        and np.max(np.abs(no_entry.state("A_ysa_bound"))) == 0.0
        and np.max(np.abs(no_entry.state("A_internalized"))) == 0.0
        and np.max(np.abs(no_entry.state("A_cargo"))) == 0.0
    )
    rows.append(
        {
            "check": "extravasation_off_blocks_tumor",
            "passed": str(tumor_blocked),
            "carrier_residual": v13_max_abs_carrier_residual(no_entry, no_entry_params),
            "cargo_residual": v13_max_abs_cargo_residual(no_entry, no_entry_params),
            "minimum_state_value": v13_minimum_state_value(no_entry),
            "detail": "With k_entry = 0, no carrier may reach the tumor interstitium or any downstream targeting state.",
        }
    )

    # 7. Organ-vascular un-uptaken carrier returns to the venous pool (no vanishing).
    return_params = dict(
        params,
        k_entry=0.0,
        k_blood_clear=0.0,
        Vmax_liver=0.0,
        Vmax_spleen=0.0,
        k_liver_elim=0.0,
        k_spleen_elim=0.0,
    )
    closed = v13_simulate(return_params, design="ocn-ysa")
    cumulative = v13_cumulative_fluxes(closed, return_params)
    arterial_out = float(cumulative["lh_to_sa"][-1])
    venous_in = float(
        cumulative["liver_venous"][-1]
        + cumulative["spleen_venous"][-1]
        + cumulative["tumor_venous"][-1]
        + cumulative["rest_venous"][-1]
    )
    rows.append(
        {
            "check": "organ_uptake_returns_to_venous",
            "passed": str(np.isclose(venous_in, arterial_out, rtol=1e-6, atol=1e-9)),
            "carrier_residual": v13_max_abs_carrier_residual(closed, return_params),
            "cargo_residual": v13_max_abs_cargo_residual(closed, return_params),
            "minimum_state_value": v13_minimum_state_value(closed),
            "detail": f"With all uptake/clearance/extravasation off, cumulative venous return ({venous_in:.6g}) must equal cumulative arterial delivery ({arterial_out:.6g}).",
        }
    )

    # Append the baseline final mass distribution as a final row for the report.
    distribution = v13_mass_distribution(baseline, params)
    rows.append(
        {
            "check": "baseline_final_mass_distribution",
            "passed": "n/a",
            "carrier_residual": v13_max_abs_carrier_residual(baseline, params),
            "cargo_residual": v13_max_abs_cargo_residual(baseline, params),
            "minimum_state_value": v13_minimum_state_value(baseline),
            "detail": (
                f"vascular {distribution['vascular_fraction_final']:.4g}, "
                f"liver {distribution['liver_fraction_final']:.4g}, "
                f"spleen {distribution['spleen_fraction_final']:.4g}, "
                f"tumor-associated {distribution['tumor_associated_fraction_final']:.4g}, "
                f"terminal {distribution['elim_deg_fraction_final']:.4g} "
                f"(lung {distribution['cumulative_lung_ret_fraction']:.4g}, "
                f"liver-elim {distribution['cumulative_liver_elim_fraction']:.4g}, "
                f"spleen-elim {distribution['cumulative_spleen_elim_fraction']:.4g})"
            ),
        }
    )

    write_csv(out / "tables" / "mass_balance_audit.csv", rows)
    write_csv(out / "tables" / "baseline_mass_distribution.csv", [distribution])
    return rows


def v13_metric_dictionary() -> list[dict[str, str]]:
    """Return the metric dictionary: stock/gross/first-event classification per metric."""

    def row(
        name: str,
        name_zh: str,
        numerator: str,
        denominator: str,
        unit: str,
        klass: str,
        can_exceed: str,
        interpretation: str,
        source: str,
        confidence: str,
        comparable: str,
    ) -> dict[str, str]:
        return {
            "metric_name": name,
            "name_zh": name_zh,
            "numerator": numerator,
            "denominator": denominator,
            "unit": unit,
            "metric_class": klass,
            "can_exceed_100pct": can_exceed,
            "biological_interpretation": interpretation,
            "parameter_source": source,
            "confidence_level": confidence,
            "experimental_comparable": comparable,
        }

    return [
        row("first_pass_lung_arrival_fraction", "首次肺到达比例", "cumulative rh_to_lung", "D0", "dimensionless",
            "first-event", "no", "首轮循环进入肺血管的剂量比例（结构恒等于约 1）", "structure", "high", "partial"),
        row("first_pass_lung_retention_fraction", "首次肺截留比例", "cumulative lung_ret", "D0", "dimensionless",
            "first-event", "no", "首轮肺循环后滞留肺部的剂量比例（≈ f_lung_ret）", "f_lung_ret (assumed_exploratory)", "low", "partial"),
        row("first_pass_systemic_availability", "首过肺后系统可用比例", "cumulative lh_to_sa", "D0", "dimensionless",
            "first-event", "no", "首次通过肺后进入系统动脉的剂量比例（≈ 1 - f_lung_ret）", "f_lung_ret (assumed_exploratory)", "low", "partial"),
        row("tumor_vascular_first_pass_arrival_fraction", "肿瘤血管首次到达比例", "cumulative sa_to_tumor", "D0", "dimensionless",
            "first-event", "no", "首次循环进入肿瘤血管的剂量比例", "Qtum/Qco flow fraction", "medium", "partial"),
        row("tumor_first_extravasation_fraction", "肿瘤首次外渗比例", "cumulative entry", "D0", "dimensionless",
            "first-event", "no", "首次循环进入肿瘤间质的剂量比例", "k_entry (assumed_exploratory)", "low", "partial"),
        row("first_pass_extravasation_efficiency", "首次外渗效率", "cumulative entry", "cumulative sa_to_tumor", "dimensionless",
            "first-event", "no", "首次进入肿瘤间质/首次到达肿瘤血管", "k_entry + Qtum", "low", "partial"),
        row("arterial_auc", "动脉AUC", "integral A_systemic_arterial dt", "n/a", "fraction*h",
            "stock (exposure)", "yes", "系统动脉内载体暴露量（stock 时间积分）", "circulation volumes", "medium", "partial"),
        row("tumor_vascular_auc", "肿瘤血管AUC", "integral A_tumor_vascular dt", "n/a", "fraction*h",
            "stock (exposure)", "yes", "肿瘤血管内载体暴露量", "V_tumor_vascular + Qtum", "low", "partial"),
        row("tumor_interstitial_auc", "肿瘤间质AUC", "integral A_tumor_free dt", "n/a", "fraction*h",
            "stock (exposure)", "yes", "肿瘤间质游离载体暴露量", "k_entry/k_return", "low", "partial"),
        row("tumor_auc", "肿瘤相关总AUC", "integral (free+OCN+YSA+int) dt", "n/a", "fraction*h",
            "stock (exposure)", "yes", "肿瘤相关状态总暴露量", "targeting params", "low", "partial"),
        row("max_tumor_fraction", "肿瘤相关峰值比例", "peak (free+OCN+YSA+int)", "D0", "dimensionless",
            "stock", "no", "某一时刻肿瘤相关状态实际存量峰值", "targeting params", "low", "partial"),
        row("tumor_retention_efficiency", "末端存量/累计外渗（旧字段）", "tumor-associated(final)", "cumulative entry (repeat entry included)", "dimensionless",
            "stock / gross throughput", "no", "非首次进入滞留概率；从核心汇报排除", "k_return/k_ocn/k_int", "low", "partial"),
        row("net_extravasated_fraction", "净外渗比例", "cumulative entry - return", "D0", "dimensionless",
            "gross", "no", "累计外渗减去累计回流（净进入间质）", "k_entry/k_return", "low", "partial"),
        row("ocn_binding_fraction_gross", "OCN结合比例(gross)", "cumulative ocn_on", "cumulative entry", "dimensionless",
            "gross", "yes", "累计OCN结合通量/累计外渗（可重复计数）", "k_ocn_on (assumed_exploratory)", "low", "no"),
        row("ysa_binding_fraction_gross", "YSA结合比例(gross)", "cumulative (ysa_from_free+ocn_to_ysa)", "cumulative entry", "dimensionless",
            "gross", "yes", "累计YSA结合通量/累计外渗（可重复计数）", "k_ysa_on/k_cross (assumed_exploratory)", "low", "no"),
        row("ocn_to_ysa_transfer_fraction_gross", "OCN→YSA转移比例(gross)", "cumulative ocn_to_ysa", "cumulative ocn_on", "dimensionless",
            "gross", "yes", "OCN结合后转移至YSA的累计通量比例", "k_cross (assumed_exploratory)", "low", "no"),
        row("internalization_efficiency_gross", "总内化/累计外渗(gross)", "cumulative (internalization + nonspecific_internalization)", "cumulative entry", "dimensionless",
            "gross", "no", "特异与非特异累计内化之和/累计外渗", "k_int + k_int_ns (assumed_exploratory)", "low", "no"),
        row("ysa_specific_internalization_fraction_gross", "YSA特异内化比例(gross)", "cumulative internalization", "cumulative entry", "dimensionless",
            "gross", "yes", "YSA介导特异内化通量/累计外渗", "k_int", "low", "no"),
        row("nonspecific_internalization_fraction_gross", "非特异内化比例(gross)", "cumulative nonspecific_internalization", "cumulative entry", "dimensionless",
            "gross", "yes", "非特异内化通量/累计外渗", "k_int_ns (assumed_exploratory)", "low", "no"),
        row("intracellular_degradation_fraction_gross", "胞内降解比例(gross)", "cumulative deg_int", "cumulative entry", "dimensionless",
            "gross", "yes", "胞内酸性降解通量/累计外渗", "k_deg_int_acid", "low", "no"),
        row("cargo_release_efficiency_gross", "cargo释放效率(gross)", "Y_cargo * cumulative deg_int", "Y_cargo * D0", "dimensionless",
            "gross", "yes", "胞内释放cargo总量/初始装载cargo", "Y_cargo + k_deg_int_acid", "low", "no"),
        row("active_cargo_remaining_fraction_final", "活性cargo剩余比例", "A_cargo(final)", "Y_cargo * D0", "dimensionless",
            "stock", "no", "末端仍存活的活性cargo比例", "k_cargo_loss", "low", "partial"),
        row("lung_fraction_final", "肺分布比例(末)", "A_lung_vascular(final)", "D0", "dimensionless",
            "stock", "no", "末端肺血管内载体比例", "f_lung_ret", "low", "partial"),
        row("liver_fraction_final", "肝分布比例(末)", "A_liver(final)", "D0", "dimensionless",
            "stock", "no", "末端肝MPS相关载体比例", "Vmax_liver/Km_liver", "low", "partial"),
        row("spleen_fraction_final", "脾分布比例(末)", "A_spleen(final)", "D0", "dimensionless",
            "stock", "no", "末端脾MPS相关载体比例", "Vmax_spleen/Km_spleen", "low", "partial"),
        row("rest_of_body_fraction_final", "其余组织比例(末)", "A_rest_vascular(final)", "D0", "dimensionless",
            "stock", "no", "末端其余组织血管内载体比例", "circulation volumes", "medium", "partial"),
        row("tumor_fraction_final", "肿瘤比例(末)", "tumor-associated(final)", "D0", "dimensionless",
            "stock", "no", "末端肿瘤相关状态载体比例", "targeting params", "low", "partial"),
        row("eliminated_fraction_final", "清除/降解比例(末)", "A_elim_deg(final)", "D0", "dimensionless",
            "stock", "no", "末端不可逆清除/降解载体比例", "clearance params", "low", "partial"),
        row("apparent_blood_half_life_h", "表观血液半衰期", "t(血管总量降至峰值一半) - t(peak)", "n/a", "h",
            "temporal", "no", "总血管载体量从峰值衰减一半所需时间", "circulation + clearance", "low", "partial"),
        row("internalized_fraction_final", "内化比例(末)", "A_internalized(final)", "D0", "dimensionless",
            "stock", "no", "末端胞内内化载体比例", "k_int", "low", "partial"),
        row("released_cargo_fraction_final", "释放cargo比例(末)", "A_cargo(final)", "Y_cargo * D0", "dimensionless",
            "stock", "no", "末端仍存活活性cargo比例", "k_deg_int_acid + k_cargo_loss", "low", "partial"),
    ]


def run_v13_metric_dictionary(output_dir: str | Path = "outputs/v1_3") -> list[dict[str, str]]:
    """Write the metric dictionary CSV."""

    rows = v13_metric_dictionary()
    write_csv(Path(output_dir) / "tables" / "metric_dictionary.csv", rows)
    return rows


def run_v13_baseline_detailed(
    output_dir: str | Path = "outputs/v1_3",
    *,
    size: str = "100nm",
    modes: ModelModes = V1_3_MODES,
) -> list[dict[str, object]]:
    """Run the baseline with the full metric set (first-event + gross + temporal)."""

    from cabomb_model.v13_circulation import v13_default_parameters

    params = v13_default_parameters()
    result = v13_simulate(params, design="ocn-ysa", size=size, modes=modes)
    first_pass = v13_simulate(params, design="ocn-ysa", size=size, modes=modes, recirculate=False)
    validate_v13_result(result, params, context="V1.3 baseline detailed")
    detailed = v13_detailed_metrics(result, params, first_pass=first_pass)
    write_csv(Path(output_dir) / "tables" / "baseline_detailed_metrics.csv", [detailed])
    return [detailed]


def run_v13_mechanism_comparison(
    output_dir: str | Path = "outputs/v1_3",
    *,
    sizes: tuple[str, ...] = SIZE_KEYS,
) -> list[dict[str, object]]:
    """Compare sequential and parallel targeting across designs and sizes."""

    from cabomb_model.v13_circulation import v13_default_parameters

    params = v13_default_parameters()
    sequential = ModelModes("saturable_v1_2", "sequential_v1_2", "explicit_v1_3")
    parallel = ModelModes("saturable_v1_2", "parallel_dual_access_v1_3", "explicit_v1_3")
    rows: list[dict[str, object]] = []
    for modes in (sequential, parallel):
        for size in sizes:
            for design_name in DESIGNS:
                result = v13_simulate(params, design=design_name, size=size, modes=modes)
                first_pass = v13_simulate(
                    params, design=design_name, size=size, modes=modes, recirculate=False
                )
                validate_v13_result(
                    result, params, context=f"mechanism {modes.targeting_model} {size} {design_name}"
                )
                detailed = v13_detailed_metrics(result, params, first_pass=first_pass)
                rows.append(detailed)
    write_csv(Path(output_dir) / "tables" / "mechanism_comparison.csv", rows)
    write_v13_mechanism_comparison(
        rows, Path(output_dir) / "figures" / "mechanism_comparison.svg"
    )
    return rows


def run_v13_feasibility_domain(
    output_dir: str | Path = "outputs/v1_3",
    *,
    size: str = "100nm",
    levels: tuple[float, float, float] = (0.5, 1.0, 1.5),
) -> list[dict[str, object]]:
    """Grid-sweep the targeting parameters to locate the dual-targeting win region.

    For each grid point we compare the OCN-YSA design against the YSA-only design
    and record whether dual targeting improves tumor AUC, final internalization,
    and final cargo release. The winning sample fraction is reported as an
    exploratory parameter-sample fraction, never as a clinical success probability.
    """

    from cabomb_model.v13_circulation import v13_default_parameters

    params = v13_default_parameters()
    grid_symbols = (
        "k_ocn_on",
        "k_cross",
        "k_ysa_on",
        "k_int",
        "R_EphA2_eff",
    )
    bounds = v13_parameter_bounds()
    rows: list[dict[str, object]] = []
    win_counts = {metric: 0 for metric in ("tumor_auc", "internalized_fraction_final", "released_cargo_fraction_final")}
    total = 0

    def sweep(index: int, current: dict[str, float]) -> None:
        nonlocal total
        if index == len(grid_symbols):
            sample = dict(params, **current)
            dual = v13_simulate(sample, design="ocn-ysa", size=size, modes=V1_3_MODES)
            ysa_only = v13_simulate(sample, design="ysa", size=size, modes=V1_3_MODES)
            dual_summary = _validated_v13_summary(dual, sample, context="V1.3 feasibility dual")
            ysa_summary = _validated_v13_summary(ysa_only, sample, context="V1.3 feasibility YSA-only")
            record = {
                "size": size,
                "targeting_model": V1_3_MODES.targeting_model,
                "mps_model": V1_3_MODES.mps_model,
                "circulation_model": V1_3_MODES.circulation_model,
                **{symbol: current[symbol] for symbol in grid_symbols},
                "dual_tumor_auc": dual_summary["tumor_auc"],
                "ysa_only_tumor_auc": ysa_summary["tumor_auc"],
                "dual_internalized_fraction_final": dual_summary["internalized_fraction_final"],
                "ysa_only_internalized_fraction_final": ysa_summary["internalized_fraction_final"],
                "dual_released_cargo_fraction_final": dual_summary["released_cargo_fraction_final"],
                "ysa_only_released_cargo_fraction_final": ysa_summary["released_cargo_fraction_final"],
            }
            for metric in win_counts:
                record[f"dual_wins_{metric}"] = str(
                    float(record[f"dual_{metric}"]) > float(record[f"ysa_only_{metric}"])
                )
            rows.append(record)
            for metric in win_counts:
                if record[f"dual_wins_{metric}"] == "True":
                    win_counts[metric] += 1
            total += 1
            return
        symbol = grid_symbols[index]
        lower, upper = bounds[symbol]
        baseline = params[symbol]
        for level in levels:
            value = min(upper, max(lower, baseline * level))
            current[symbol] = value
            sweep(index + 1, current)

    sweep(0, {})
    summary_row: dict[str, object] = {
        "size": size,
        "grid_points": total,
        "grid_dimension": len(grid_symbols),
        "levels_per_parameter": len(levels),
        "grid_symbols": ";".join(grid_symbols),
        "targeting_model": V1_3_MODES.targeting_model,
        "mps_model": V1_3_MODES.mps_model,
        "circulation_model": V1_3_MODES.circulation_model,
        "dual_wins_tumor_auc_fraction": float(win_counts["tumor_auc"] / total) if total else 0.0,
        "dual_wins_internalized_fraction": float(win_counts["internalized_fraction_final"] / total) if total else 0.0,
        "dual_wins_released_cargo_fraction": float(win_counts["released_cargo_fraction_final"] / total) if total else 0.0,
    }
    write_csv(Path(output_dir) / "tables" / "feasibility_domain.csv", rows)
    write_csv(Path(output_dir) / "tables" / "feasibility_domain_summary.csv", [summary_row])
    return rows


def run_v13_all_analyses(
    output_dir: str | Path = "outputs/v1_3",
) -> dict[str, object]:
    """Run every V1.3 analysis and write all tables and figures."""

    baseline = run_v13_baseline(output_dir)
    design_rows = run_v13_design_comparison(output_dir)
    size_rows = run_v13_size_comparison(output_dir)
    scenario_rows = run_v13_epr_mps_scenarios(output_dir)
    sensitivity_rows = run_v13_local_sensitivity(output_dir)
    monte_carlo_rows = run_v13_monte_carlo(output_dir)
    audit_rows = run_v13_mass_balance_audit(output_dir)
    detailed_rows = run_v13_baseline_detailed(output_dir)
    dictionary_rows = run_v13_metric_dictionary(output_dir)
    mechanism_rows = run_v13_mechanism_comparison(output_dir)
    feasibility_rows = run_v13_feasibility_domain(output_dir)
    return {
        "baseline": baseline,
        "design_comparison": design_rows,
        "size_comparison": size_rows,
        "scenario_comparison": scenario_rows,
        "local_sensitivity": sensitivity_rows,
        "monte_carlo": monte_carlo_rows,
        "mass_balance_audit": audit_rows,
        "baseline_detailed": detailed_rows,
        "metric_dictionary": dictionary_rows,
        "mechanism_comparison": mechanism_rows,
        "feasibility_domain": feasibility_rows,
    }
