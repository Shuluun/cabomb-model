"""V1.3 explicit-circulation ODE, solver, and mass-ledger metrics.

This module implements the optional ``explicit_v1_3`` circulation mode. It is
separate from ``model.py``/``solver.py`` so the frozen V1.1/V1.2 ten-state path
is never modified. The tumor-targeting and saturable-MPS terms are the V1.2
forms re-expressed over the nine vascular states.
"""

from __future__ import annotations

import math
from dataclasses import dataclass
from typing import Mapping, Sequence

import numpy as np

from cabomb_model.designs import CarrierDesign, get_design
from cabomb_model.fluxes import compute_k_entry, degradation_rate
from cabomb_model.model import assert_nonnegative_parameters
from cabomb_model.model_modes import ModelModes, V1_3_MODES
from cabomb_model.parameters import (
    V13_EXECUTABLE_SYMBOLS,
    V13_FLOW_SYMBOLS,
    V13_LUNG_RET_SYMBOLS,
    V13_SIZE_EXTRA_SYMBOLS,
    V13_SIZE_INT_SYMBOLS,
    V13_SIZE_MPS_SYMBOLS,
    V13_VOLUME_SYMBOLS,
    default_model_parameters,
    load_registry,
)
from cabomb_model.scenarios import Scenario, make_scenario
from cabomb_model.v13_states import (
    V13_ACTIVE_CARRIER_STATES,
    V13_CARRIER_STATES,
    V13_STATE_INDEX,
    V13_STATES,
    V13_TUMOR_ASSOCIATED_STATES,
    V13_VASCULAR_STATES,
    v13_initial_state,
)

SIZE_KEYS: tuple[str, ...] = ("100nm", "200nm", "500nm", "800nm", "2um", "3um")

_LUNG_RET_KEYS = ("f_lung_ret",)
_MPS_KEYS = ("liver_mps", "spleen_mps", "liver_release", "spleen_release")
_CIRCULATION_KEYS = (
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
)
_TUMOR_KEYS = (
    "entry",
    "return",
    "ocn_on",
    "ocn_off",
    "ocn_loss",
    "ysa_from_free",
    "ocn_to_ysa",
    "ysa_off",
    "internalization",
    "nonspecific_internalization",
    "deg_int",
    "tumor_deg",
)
_TERMINAL_KEYS = ("blood_clear", "liver_elim", "spleen_elim")


def v13_default_parameters(path: str | None = None) -> dict[str, float]:
    """Return the shared registry parameters plus the V1.3-only parameters."""

    params = default_model_parameters(path)
    records = load_registry(path)
    by_symbol = {row["parameter_symbol"]: row for row in records}
    for symbol in V13_EXECUTABLE_SYMBOLS:
        if symbol not in by_symbol:
            raise ValueError(f"Executable parameter {symbol} is absent from registry.")
        params[symbol] = float(by_symbol[symbol]["baseline_value"])
    return params


def validate_v13_parameters(params: Mapping[str, float]) -> None:
    """Validate V1.3-only parameter domains not expressible as non-negativity."""

    for symbol in V13_VOLUME_SYMBOLS:
        value = float(params[symbol])
        if not math.isfinite(value) or value <= 0.0:
            raise ValueError(f"{symbol} must be finite and positive.")
    for symbol in V13_FLOW_SYMBOLS:
        value = float(params[symbol])
        if not math.isfinite(value) or value <= 0.0:
            raise ValueError(f"{symbol} must be finite and positive.")
    for symbol in V13_LUNG_RET_SYMBOLS:
        value = float(params[symbol])
        if not math.isfinite(value) or not 0.0 <= value < 1.0:
            raise ValueError(f"{symbol} must be finite and within [0, 1).")
    for symbol in V13_SIZE_MPS_SYMBOLS + V13_SIZE_EXTRA_SYMBOLS + V13_SIZE_INT_SYMBOLS:
        value = float(params[symbol])
        if not math.isfinite(value) or value < 0.0:
            raise ValueError(f"{symbol} must be finite and nonnegative.")
    alpha = float(params["alpha_ysa_direct_access"])
    if not math.isfinite(alpha) or not 0.0 <= alpha <= 1.0:
        raise ValueError("alpha_ysa_direct_access must be finite and within [0, 1].")


def v13_fluxes(
    state: Sequence[float],
    params: Mapping[str, float],
    design: str | CarrierDesign,
    scenario: Scenario,
    size: str,
    targeting_model: str = "sequential_v1_2",
) -> dict[str, float]:
    """Evaluate all V1.3 process fluxes for one state vector (hours).

    ``targeting_model`` selects between ``sequential_v1_2`` (tumor_free -> OCN ->
    YSA only) and ``parallel_dual_access_v1_3`` (tumor_free -> YSA directly in
    parallel with tumor_free -> OCN -> YSA). Circulation and entry fluxes are
    identical for both targeting models.
    """

    y = np.asarray(state, dtype=float)
    idx = V13_STATE_INDEX
    carrier_design = get_design(design)
    i_ocn = carrier_design.i_ocn
    i_ysa = carrier_design.i_ysa

    qco = params["Qco"] * 60.0
    qliv = params["Qliv"] * 60.0
    qspl = params["Qspl"] * 60.0
    qtum = params["Qtum"] * 60.0
    qrest = qco - qliv - qspl - qtum

    c_cv = y[idx["A_central_venous"]] / params["V_central_venous"]
    c_rh = y[idx["A_right_heart"]] / params["V_right_heart"]
    c_lung = y[idx["A_lung_vascular"]] / params["V_lung_vascular"]
    c_lh = y[idx["A_left_heart"]] / params["V_left_heart"]
    c_sa = y[idx["A_systemic_arterial"]] / params["V_systemic_arterial"]
    c_liv_v = y[idx["A_liver_vascular"]] / params["V_liver_vascular"]
    c_spl_v = y[idx["A_spleen_vascular"]] / params["V_spleen_vascular"]
    c_tum_v = y[idx["A_tumor_vascular"]] / params["V_tumor_vascular"]
    c_rest_v = y[idx["A_rest_vascular"]] / params["V_rest_vascular"]

    f_lung_ret = params[f"f_lung_ret_{size}"]
    size_mps = params[f"size_mps_{size}"]
    size_extra = params[f"size_extra_{size}"]
    size_int = params[f"size_int_{size}"]

    a_liver = y[idx["A_liver"]]
    a_spleen = y[idx["A_spleen"]]
    a_tumor_free = y[idx["A_tumor_free"]]
    a_ocn = y[idx["A_ocn_bound"]]
    a_ysa = y[idx["A_ysa_bound"]]
    a_int = y[idx["A_internalized"]]
    a_rest_v = y[idx["A_rest_vascular"]]

    cv_to_rh = qco * c_cv
    rh_to_lung = qco * c_rh
    lung_pass = qco * c_lung * (1.0 - f_lung_ret)
    lung_ret = qco * c_lung * f_lung_ret
    lh_to_sa = qco * c_lh
    sa_to_liver = qliv * c_sa
    sa_to_spleen = qspl * c_sa
    sa_to_tumor = qtum * c_sa
    sa_to_rest = qrest * c_sa
    liver_venous = qliv * c_liv_v
    spleen_venous = qspl * c_spl_v
    tumor_venous = qtum * c_tum_v
    rest_venous = qrest * c_rest_v

    liver_mps = (
        scenario.mps_factor
        * size_mps
        * params["Vmax_liver"]
        * c_liv_v
        / (params["Km_liver"] + c_liv_v)
    )
    spleen_mps = (
        scenario.mps_factor
        * size_mps
        * params["Vmax_spleen"]
        * c_spl_v
        / (params["Km_spleen"] + c_spl_v)
    )
    liver_release = params["k_liver_release"] * a_liver
    spleen_release = params["k_spleen_release"] * a_spleen

    # Permeability-surface clearance PS_eff = k_entry_eff * size_extra * Vblood,
    # so extravasation is flow/permeability-limited on the tumor-vascular
    # concentration (E = PS_eff / (Qtum + PS_eff)) and reproduces the V1.2
    # lumped tumor-entry scale in the well-mixed limit.
    entry = compute_k_entry(params, scenario) * size_extra * params["Vblood"] * c_tum_v
    return_flux = params["k_return"] * a_tumor_free

    ocn_on = i_ocn * params["k_ocn_on"] * a_tumor_free
    ocn_off = i_ocn * params["k_ocn_off"] * a_ocn
    ocn_loss = i_ocn * params["k_ocn_loss"] * a_ocn
    ysa_rate = params["k_ysa_on"] * params["R_EphA2_eff"]
    if targeting_model == "parallel_dual_access_v1_3":
        # Direct free -> YSA access runs in parallel with the OCN-first route.
        # The direct-access coefficient is full for the YSA-only design and
        # governed by alpha_ysa_direct_access for the dual (OCN-YSA) design so
        # the parallel assumption is an independent, interpretable parameter
        # rather than a hidden re-tuning of k_ysa_on.
        direct_access = (1 - i_ocn) + i_ocn * params["alpha_ysa_direct_access"]
        ysa_from_free = i_ysa * direct_access * ysa_rate * a_tumor_free
    else:  # sequential_v1_2
        ysa_from_free = i_ysa * (1 - i_ocn) * ysa_rate * a_tumor_free
    ocn_to_ysa = i_ocn * i_ysa * params["k_cross"] * params["R_EphA2_eff"] * a_ocn
    ysa_off = i_ysa * params["k_ysa_off"] * a_ysa
    internalization = i_ysa * params["k_int"] * size_int * a_ysa
    nonspecific_internalization = params["k_int_ns"] * size_int * a_tumor_free
    deg_int = degradation_rate(scenario.intracellular_pH, params) * a_int
    tumor_deg = params["k_tumor_deg"] * a_tumor_free

    blood_clear = params["k_blood_clear"] * a_rest_v
    liver_elim = params["k_liver_elim"] * a_liver
    spleen_elim = params["k_spleen_elim"] * a_spleen

    return {
        "cv_to_rh": cv_to_rh,
        "rh_to_lung": rh_to_lung,
        "lung_pass": lung_pass,
        "lung_ret": lung_ret,
        "lh_to_sa": lh_to_sa,
        "sa_to_liver": sa_to_liver,
        "sa_to_spleen": sa_to_spleen,
        "sa_to_tumor": sa_to_tumor,
        "sa_to_rest": sa_to_rest,
        "liver_venous": liver_venous,
        "spleen_venous": spleen_venous,
        "tumor_venous": tumor_venous,
        "rest_venous": rest_venous,
        "liver_mps": liver_mps,
        "spleen_mps": spleen_mps,
        "liver_release": liver_release,
        "spleen_release": spleen_release,
        "entry": entry,
        "return": return_flux,
        "ocn_on": ocn_on,
        "ocn_off": ocn_off,
        "ocn_loss": ocn_loss,
        "ysa_from_free": ysa_from_free,
        "ocn_to_ysa": ocn_to_ysa,
        "ysa_off": ysa_off,
        "internalization": internalization,
        "nonspecific_internalization": nonspecific_internalization,
        "deg_int": deg_int,
        "tumor_deg": tumor_deg,
        "blood_clear": blood_clear,
        "liver_elim": liver_elim,
        "spleen_elim": spleen_elim,
    }


def v13_rhs(
    t: float,
    y: Sequence[float],
    params: Mapping[str, float],
    design: str | CarrierDesign = "ocn-ysa",
    scenario: Scenario | None = None,
    size: str = "100nm",
    modes: ModelModes = V1_3_MODES,
    recirculate: bool = True,
) -> np.ndarray:
    """Evaluate the V1.3 eighteen-state ODE right-hand side.

    ``recirculate=False`` diverts venous return to the terminal sink so the
    injected bolus transits the circulation exactly once; this single-pass mode
    is used only to compute first-event fractions (first lung pass, first
    systemic availability, first tumor-vessel arrival, first extravasation).
    """

    del t
    run_scenario = scenario or make_scenario(params)
    f = v13_fluxes(y, params, design, run_scenario, size, modes.targeting_model)
    idx = V13_STATE_INDEX
    dy = np.zeros(len(V13_STATES), dtype=float)

    venous_return = (
        f["liver_venous"]
        + f["spleen_venous"]
        + f["tumor_venous"]
        + f["rest_venous"]
        + f["return"]
    )
    if recirculate:
        dy[idx["A_central_venous"]] = -f["cv_to_rh"] + venous_return
    else:
        dy[idx["A_central_venous"]] = -f["cv_to_rh"]
    dy[idx["A_right_heart"]] = f["cv_to_rh"] - f["rh_to_lung"]
    dy[idx["A_lung_vascular"]] = f["rh_to_lung"] - f["lung_pass"] - f["lung_ret"]
    dy[idx["A_left_heart"]] = f["lung_pass"] - f["lh_to_sa"]
    dy[idx["A_systemic_arterial"]] = (
        f["lh_to_sa"] - f["sa_to_liver"] - f["sa_to_spleen"] - f["sa_to_tumor"] - f["sa_to_rest"]
    )
    dy[idx["A_liver_vascular"]] = (
        f["sa_to_liver"] + f["liver_release"] - f["liver_venous"] - f["liver_mps"]
    )
    dy[idx["A_spleen_vascular"]] = (
        f["sa_to_spleen"] + f["spleen_release"] - f["spleen_venous"] - f["spleen_mps"]
    )
    dy[idx["A_tumor_vascular"]] = f["sa_to_tumor"] - f["tumor_venous"] - f["entry"]
    dy[idx["A_rest_vascular"]] = f["sa_to_rest"] - f["rest_venous"] - f["blood_clear"]

    dy[idx["A_liver"]] = f["liver_mps"] - f["liver_release"] - f["liver_elim"]
    dy[idx["A_spleen"]] = f["spleen_mps"] - f["spleen_release"] - f["spleen_elim"]

    dy[idx["A_tumor_free"]] = (
        f["entry"]
        - f["return"]
        - f["ocn_on"]
        + f["ocn_off"]
        - f["ysa_from_free"]
        + f["ysa_off"]
        - f["nonspecific_internalization"]
        - f["tumor_deg"]
    )
    dy[idx["A_ocn_bound"]] = f["ocn_on"] - f["ocn_off"] - f["ocn_loss"] - f["ocn_to_ysa"]
    dy[idx["A_ysa_bound"]] = f["ysa_from_free"] + f["ocn_to_ysa"] - f["ysa_off"] - f["internalization"]
    dy[idx["A_internalized"]] = f["internalization"] + f["nonspecific_internalization"] - f["deg_int"]

    dy[idx["A_elim_deg"]] = (
        f["lung_ret"]
        + f["blood_clear"]
        + f["liver_elim"]
        + f["spleen_elim"]
        + f["tumor_deg"]
        + f["ocn_loss"]
        + f["deg_int"]
        + (0.0 if recirculate else venous_return)
    )
    dy[idx["A_cargo"]] = params["Y_cargo"] * f["deg_int"] - params["k_cargo_loss"] * y[idx["A_cargo"]]
    dy[idx["A_cargo_lost"]] = (
        params["Y_cargo"]
        * (
            f["lung_ret"]
            + f["blood_clear"]
            + f["liver_elim"]
            + f["spleen_elim"]
            + f["tumor_deg"]
            + f["ocn_loss"]
        )
        + params["k_cargo_loss"] * y[idx["A_cargo"]]
    )
    return dy


def v13_time_grid(t_end_h: float = 72.0) -> np.ndarray:
    """Dense early grid for the first pulmonary pass and arterial peak."""

    early = np.linspace(0.0, 0.05, 400)
    late = np.linspace(0.05, float(t_end_h), 240)[1:]
    return np.concatenate((early, late))


@dataclass(frozen=True)
class V13Result:
    """Container for one V1.3 explicit-circulation simulation."""

    t: np.ndarray
    y: np.ndarray
    state_names: tuple[str, ...]
    design: CarrierDesign
    scenario: Scenario
    modes: ModelModes
    size: str
    solver_method: str
    success: bool
    message: str
    rtol: float
    atol: float

    def state(self, name: str) -> np.ndarray:
        return self.y[V13_STATE_INDEX[name], :]

    def as_rows(self) -> list[dict[str, float | str]]:
        rows: list[dict[str, float | str]] = []
        for col, time_h in enumerate(self.t):
            row: dict[str, float | str] = {
                "time_h": float(time_h),
                "design": self.design.name,
                "scenario": self.scenario.name,
                "size": self.size,
                "solver_method": self.solver_method,
            }
            for state_name in self.state_names:
                row[state_name] = float(self.y[V13_STATE_INDEX[state_name], col])
            rows.append(row)
        return rows


def v13_simulate(
    params: Mapping[str, float] | None = None,
    *,
    design: str | CarrierDesign = "ocn-ysa",
    scenario: Scenario | None = None,
    modes: ModelModes | None = None,
    size: str = "100nm",
    epr: str = "medium",
    mps: str = "medium",
    t_end_h: float = 72.0,
    method: str = "LSODA",
    rtol: float = 1.0e-8,
    atol: float = 1.0e-10,
    recirculate: bool = True,
) -> V13Result:
    """Run one V1.3 explicit-circulation simulation with SciPy solve_ivp."""

    run_params = v13_default_parameters() if params is None else dict(params)
    run_modes = modes or V1_3_MODES
    validate_v13_parameters(run_params)
    assert_nonnegative_parameters(run_params)
    if size not in SIZE_KEYS:
        raise ValueError(f"Unknown size scenario {size!r}. Valid sizes: {', '.join(SIZE_KEYS)}")
    run_design = get_design(design)
    run_scenario = scenario or make_scenario(run_params, epr=epr, mps=mps)
    t_eval = v13_time_grid(t_end_h)
    y0 = np.asarray(v13_initial_state(run_params["D0"]), dtype=float)

    from scipy.integrate import solve_ivp  # type: ignore

    solution = solve_ivp(
        fun=lambda t, y: v13_rhs(
            t, y, run_params, run_design, run_scenario, size, run_modes, recirculate
        ),
        t_span=(0.0, float(t_end_h)),
        y0=y0,
        method=method,
        t_eval=t_eval,
        rtol=rtol,
        atol=atol,
    )
    return V13Result(
        t=t_eval,
        y=np.asarray(solution.y, dtype=float),
        state_names=V13_STATES,
        design=run_design,
        scenario=run_scenario,
        modes=run_modes,
        size=size,
        solver_method=method,
        success=bool(solution.success),
        message=str(solution.message),
        rtol=rtol,
        atol=atol,
    )


def _sum_states(result: V13Result, names: tuple[str, ...]) -> np.ndarray:
    indexes = [V13_STATE_INDEX[name] for name in names]
    return np.sum(result.y[indexes, :], axis=0)


def v13_carrier_total(result: V13Result) -> np.ndarray:
    return _sum_states(result, V13_CARRIER_STATES)


def v13_active_carrier_total(result: V13Result) -> np.ndarray:
    return _sum_states(result, V13_ACTIVE_CARRIER_STATES)


def v13_tumor_associated(result: V13Result) -> np.ndarray:
    return _sum_states(result, V13_TUMOR_ASSOCIATED_STATES)


def v13_cargo_total(result: V13Result, params: Mapping[str, float]) -> np.ndarray:
    return (
        params["Y_cargo"] * v13_active_carrier_total(result)
        + result.state("A_cargo")
        + result.state("A_cargo_lost")
    )


def v13_max_abs_carrier_residual(result: V13Result, params: Mapping[str, float]) -> float:
    return float(np.max(np.abs(v13_carrier_total(result) - params["D0"])))


def v13_max_abs_cargo_residual(result: V13Result, params: Mapping[str, float]) -> float:
    return float(
        np.max(np.abs(v13_cargo_total(result, params) - params["Y_cargo"] * params["D0"]))
    )


def v13_minimum_state_value(result: V13Result) -> float:
    return float(np.min(result.y))


def v13_auc(result: V13Result, state_name: str) -> float:
    return float(np.trapezoid(result.state(state_name), result.t))


def v13_cumulative_fluxes(
    result: V13Result, params: Mapping[str, float]
) -> dict[str, np.ndarray]:
    """Postprocess gross cumulative process fluxes on the saved time grid."""

    keys = _CIRCULATION_KEYS + _MPS_KEYS + _TUMOR_KEYS + _TERMINAL_KEYS
    values: dict[str, np.ndarray] = {key: np.zeros_like(result.t) for key in keys}
    for index in range(result.t.size):
        state = np.maximum(result.y[:, index], 0.0)
        flux = v13_fluxes(
            state, params, result.design, result.scenario, result.size,
            result.modes.targeting_model,
        )
        for key in keys:
            values[key][index] = max(0.0, float(flux[key]))
    out: dict[str, np.ndarray] = {}
    for key, series in values.items():
        increments = 0.5 * (series[1:] + series[:-1]) * np.diff(result.t)
        out[key] = np.concatenate((np.array([0.0]), np.cumsum(increments)))
    return out


def v13_circulation_metrics(
    result: V13Result, params: Mapping[str, float]
) -> dict[str, float]:
    """Return the V1.3 first-pass, delivery, extravasation, and AUC metrics."""

    dose = params["D0"]
    cumulative = v13_cumulative_fluxes(result, params)
    net_extravasated = float(cumulative["entry"][-1] - cumulative["return"][-1])
    return {
        "fraction_surviving_first_lung_pass": 1.0 - params[f"f_lung_ret_{result.size}"],
        "cumulative_lung_retention_fraction": float(cumulative["lung_ret"][-1] / dose) if dose else 0.0,
        "gross_arterial_delivery_fraction": float(cumulative["sa_to_tumor"][-1] / dose) if dose else 0.0,
        "gross_extravasated_fraction": float(cumulative["entry"][-1] / dose) if dose else 0.0,
        "net_extravasated_fraction": float(net_extravasated / dose) if dose else 0.0,
        "arterial_auc": v13_auc(result, "A_systemic_arterial"),
        "tumor_vascular_auc": v13_auc(result, "A_tumor_vascular"),
        "tumor_interstitial_auc": v13_auc(result, "A_tumor_free"),
    }


def v13_summary_metrics(
    result: V13Result, params: Mapping[str, float]
) -> dict[str, float | str]:
    """Summarize one V1.3 simulation for comparison tables."""

    dose = params["D0"]
    tumor = v13_tumor_associated(result)
    base: dict[str, float | str] = {
        "design": result.design.name,
        "scenario": result.scenario.name,
        "size": result.size,
        "mps_model": result.modes.mps_model,
        "targeting_model": result.modes.targeting_model,
        "circulation_model": result.modes.circulation_model,
        "solver_method": result.solver_method,
        "success": str(result.success),
        "max_tumor_fraction": float(np.max(tumor) / dose) if dose else 0.0,
        "tumor_auc": float(np.trapezoid(tumor, result.t)),
        "internalized_fraction_final": float(result.state("A_internalized")[-1] / dose) if dose else 0.0,
        "released_cargo_fraction_final": float(result.state("A_cargo")[-1] / params["Y_cargo"] / dose) if dose and params["Y_cargo"] else 0.0,
        "elim_deg_fraction_final": float(result.state("A_elim_deg")[-1] / dose) if dose else 0.0,
        "max_abs_carrier_residual": v13_max_abs_carrier_residual(result, params),
        "max_abs_cargo_residual": v13_max_abs_cargo_residual(result, params),
        "minimum_state_value": v13_minimum_state_value(result),
    }
    base.update(v13_circulation_metrics(result, params))
    return base


def v13_vascular_total(result: V13Result) -> np.ndarray:
    """Return the total carrier across all nine vascular states over time."""

    return _sum_states(result, V13_VASCULAR_STATES)


def _peak_and_time(series: np.ndarray, time: np.ndarray) -> tuple[float, float]:
    index = int(np.argmax(series))
    return float(series[index]), float(time[index])


def _apparent_half_life(series: np.ndarray, time: np.ndarray) -> float:
    """Return time after peak for ``series`` to first fall to half its peak."""

    peak_value = float(np.max(series))
    if peak_value <= 0.0:
        return float("nan")
    peak_index = int(np.argmax(series))
    after = np.where(series[peak_index:] <= 0.5 * peak_value)[0]
    if after.size == 0:
        return float("nan")
    return float(time[peak_index + int(after[0])] - time[peak_index])


def v13_first_pass_metrics(result: V13Result, params: Mapping[str, float]) -> dict[str, float]:
    """First-event fractions from a single-pass (``recirculate=False``) run.

    Because the single-pass run never recirculates venous return, its cumulative
    fluxes are unique-particle first-event fractions rather than gross throughput.
    """

    dose = params["D0"]
    cumulative = v13_cumulative_fluxes(result, params)
    arrival = float(cumulative["sa_to_tumor"][-1])
    extravasated = float(cumulative["entry"][-1])
    return {
        "first_pass_lung_arrival_fraction": float(cumulative["rh_to_lung"][-1] / dose) if dose else 0.0,
        "first_pass_lung_retention_fraction": float(cumulative["lung_ret"][-1] / dose) if dose else 0.0,
        "first_pass_systemic_availability": float(cumulative["lh_to_sa"][-1] / dose) if dose else 0.0,
        "tumor_vascular_first_pass_arrival_fraction": float(arrival / dose) if dose else 0.0,
        "tumor_first_extravasation_fraction": float(extravasated / dose) if dose else 0.0,
        "first_pass_extravasation_efficiency": float(extravasated / arrival) if arrival else 0.0,
    }


def v13_organ_distribution_metrics(
    result: V13Result, params: Mapping[str, float]
) -> dict[str, float]:
    """Final and peak organ/compartment fractions plus apparent vascular half-life."""

    dose = params["D0"]

    def frac(value: float) -> float:
        return float(value / dose) if dose else 0.0

    vascular = v13_vascular_total(result)
    tumor = v13_tumor_associated(result)
    lung = result.state("A_lung_vascular")
    liver = result.state("A_liver")
    spleen = result.state("A_spleen")
    rest = result.state("A_rest_vascular")
    elim = result.state("A_elim_deg")

    lung_peak, lung_tpeak = _peak_and_time(lung, result.t)
    liver_peak, liver_tpeak = _peak_and_time(liver, result.t)
    spleen_peak, spleen_tpeak = _peak_and_time(spleen, result.t)
    rest_peak, _ = _peak_and_time(rest, result.t)
    tumor_peak, tumor_tpeak = _peak_and_time(tumor, result.t)
    vascular_peak, vascular_tpeak = _peak_and_time(vascular, result.t)

    return {
        "vascular_peak_fraction": frac(vascular_peak),
        "vascular_peak_time_h": vascular_tpeak,
        "apparent_blood_half_life_h": _apparent_half_life(vascular, result.t),
        "lung_fraction_final": frac(lung[-1]),
        "lung_fraction_peak": frac(lung_peak),
        "lung_peak_time_h": lung_tpeak,
        "liver_fraction_final": frac(liver[-1]),
        "liver_fraction_peak": frac(liver_peak),
        "liver_peak_time_h": liver_tpeak,
        "spleen_fraction_final": frac(spleen[-1]),
        "spleen_fraction_peak": frac(spleen_peak),
        "spleen_peak_time_h": spleen_tpeak,
        "rest_of_body_fraction_final": frac(rest[-1]),
        "rest_of_body_fraction_peak": frac(rest_peak),
        "tumor_fraction_final": frac(tumor[-1]),
        "tumor_fraction_peak": frac(tumor_peak),
        "tumor_peak_time_h": tumor_tpeak,
        "eliminated_fraction_final": frac(elim[-1]),
    }


def v13_binding_internalization_metrics(
    result: V13Result, params: Mapping[str, float]
) -> dict[str, float]:
    """Gross (throughput) binding, internalization, degradation, and cargo fractions.

    Binding ratios are cumulative process fluxes normalized to gross entry and
    can count repeated binding/unbinding; they are not unique-particle binding
    probabilities. Internalization is irreversible in this model and its total
    includes both the YSA-specific and nonspecific paths. Cargo release is
    normalized to initial cargo load and is distinct from final active stock.
    """

    dose = params["D0"]
    cargo0 = params["Y_cargo"] * dose
    cumulative = v13_cumulative_fluxes(result, params)
    gross_entry = float(cumulative["entry"][-1])

    def per_entry(flux_key: str) -> float:
        return float(cumulative[flux_key][-1] / gross_entry) if gross_entry else 0.0

    ocn_on = float(cumulative["ocn_on"][-1])
    ysa_total = float(cumulative["ysa_from_free"][-1] + cumulative["ocn_to_ysa"][-1])
    return {
        "ocn_binding_fraction_gross": float(ocn_on / gross_entry) if gross_entry else 0.0,
        "ysa_binding_fraction_gross": float(ysa_total / gross_entry) if gross_entry else 0.0,
        "ocn_to_ysa_transfer_fraction_gross": (
            float(cumulative["ocn_to_ysa"][-1] / ocn_on) if ocn_on else 0.0
        ),
        "internalization_efficiency_gross": (
            per_entry("internalization") + per_entry("nonspecific_internalization")
        ),
        "cumulative_specific_internalization_fraction": float(cumulative["internalization"][-1] / dose) if dose else 0.0,
        "cumulative_nonspecific_internalization_fraction": float(cumulative["nonspecific_internalization"][-1] / dose) if dose else 0.0,
        "cumulative_total_internalization_fraction": float(
            (cumulative["internalization"][-1] + cumulative["nonspecific_internalization"][-1]) / dose
        ) if dose else 0.0,
        "ysa_specific_internalization_fraction_gross": per_entry("internalization"),
        "nonspecific_internalization_fraction_gross": per_entry("nonspecific_internalization"),
        "intracellular_degradation_fraction_gross": per_entry("deg_int"),
        "cargo_release_efficiency_gross": float(cumulative["deg_int"][-1] * params["Y_cargo"] / cargo0) if cargo0 else 0.0,
        "active_cargo_remaining_fraction_final": float(result.state("A_cargo")[-1] / cargo0) if cargo0 else 0.0,
    }


def v13_retention_efficiency(result: V13Result, params: Mapping[str, float]) -> float:
    """Legacy ratio: final tumor-associated amount / gross cumulative entry.

    The denominator includes repeat entry. This is not first-entry retention or
    a unique-particle probability; retained for compatibility, not core reporting.
    """

    dose = params["D0"]
    cumulative = v13_cumulative_fluxes(result, params)
    extravasated = float(cumulative["entry"][-1])
    tumor_final = float(v13_tumor_associated(result)[-1])
    return float(tumor_final / extravasated) if extravasated else 0.0


def v13_temporal_state_metrics(result: V13Result, params: Mapping[str, float]) -> dict[str, float]:
    """Peak, time-to-peak, and 24/48/72 h values for the key reporting states."""

    dose = params["D0"]

    def frac(series: np.ndarray, index: int) -> float:
        return float(series[index] / dose) if dose else 0.0

    out: dict[str, float] = {}
    for state_name in ("A_central_venous", "A_systemic_arterial", "A_tumor_vascular"):
        series = result.state(state_name)
        peak, tpeak = _peak_and_time(series, result.t)
        out[f"{state_name}_peak_fraction"] = frac(series, int(np.argmax(series)))
        out[f"{state_name}_time_to_peak_h"] = tpeak
    for state_name in ("A_tumor_free", "A_internalized", "A_cargo"):
        series = result.state(state_name)
        out[f"{state_name}_peak"] = float(np.max(series))
        out[f"{state_name}_time_to_peak_h"] = float(result.t[int(np.argmax(series))])
    for horizon in (24.0, 48.0, 72.0):
        index = int(np.argmin(np.abs(result.t - horizon)))
        out[f"tumor_associated_fraction_{int(horizon)}h"] = frac(v13_tumor_associated(result), index)
        out[f"internalized_fraction_{int(horizon)}h"] = frac(result.state("A_internalized"), index)
        out[f"active_cargo_fraction_{int(horizon)}h"] = float(
            result.state("A_cargo")[index] / (params["Y_cargo"] * dose)
        ) if dose and params["Y_cargo"] else 0.0
    return out


def v13_detailed_metrics(
    result: V13Result, params: Mapping[str, float], *, first_pass: V13Result | None = None
) -> dict[str, float | str]:
    """Return the full V1.3 metric set (first-event, gross, stock, temporal)."""

    base: dict[str, float | str] = {
        "design": result.design.name,
        "scenario": result.scenario.name,
        "size": result.size,
        "targeting_model": result.modes.targeting_model,
        "circulation_model": result.modes.circulation_model,
    }
    summary = v13_summary_metrics(result, params)
    for key in (
        "max_tumor_fraction",
        "tumor_auc",
        "internalized_fraction_final",
        "released_cargo_fraction_final",
        "elim_deg_fraction_final",
        "net_extravasated_fraction",
    ):
        base[key] = summary[key]
    base.update(v13_circulation_metrics(result, params))
    base.update(v13_organ_distribution_metrics(result, params))
    base.update(v13_binding_internalization_metrics(result, params))
    base.update(v13_temporal_state_metrics(result, params))
    base["tumor_retention_efficiency"] = v13_retention_efficiency(result, params)
    if first_pass is not None:
        base.update(v13_first_pass_metrics(first_pass, params))
    return base
