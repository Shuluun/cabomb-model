"""Mass ledgers and output metrics for CABOMB simulations."""

from __future__ import annotations

from typing import Mapping

import numpy as np

from cabomb_model.fluxes import cumulative_flux_metrics
from cabomb_model.solver import SimulationResult
from cabomb_model.states import ACTIVE_CARRIER_STATES, CARRIER_STATES, STATE_INDEX


def _sum_states(result: SimulationResult, names: tuple[str, ...]) -> np.ndarray:
    indexes = [STATE_INDEX[name] for name in names]
    return np.sum(result.y[indexes, :], axis=0)


def carrier_total(result: SimulationResult) -> np.ndarray:
    """Return carrier ledger total over time."""

    return _sum_states(result, CARRIER_STATES)


def active_carrier_total(result: SimulationResult) -> np.ndarray:
    """Return intact active carrier over time, excluding eliminated/degraded."""

    return _sum_states(result, ACTIVE_CARRIER_STATES)


def cargo_total(result: SimulationResult, params: Mapping[str, float]) -> np.ndarray:
    """Return separate payload ledger total over time."""

    return (
        params["Y_cargo"] * active_carrier_total(result)
        + result.state("A_cargo")
        + result.state("A_cargo_lost")
    )


def carrier_mass_residual(result: SimulationResult, params: Mapping[str, float]) -> np.ndarray:
    """Return carrier ledger residual relative to the initial dose."""

    return carrier_total(result) - params["D0"]


def cargo_mass_residual(result: SimulationResult, params: Mapping[str, float]) -> np.ndarray:
    """Return cargo ledger residual relative to initial payload."""

    return cargo_total(result, params) - params["Y_cargo"] * params["D0"]


def max_abs_carrier_residual(result: SimulationResult, params: Mapping[str, float]) -> float:
    """Return maximum absolute carrier mass residual."""

    return float(np.max(np.abs(carrier_mass_residual(result, params))))


def max_abs_cargo_residual(result: SimulationResult, params: Mapping[str, float]) -> float:
    """Return maximum absolute cargo mass residual."""

    return float(np.max(np.abs(cargo_mass_residual(result, params))))


def minimum_state_value(result: SimulationResult) -> float:
    """Return the minimum state value over all states and times."""

    return float(np.min(result.y))


def auc(result: SimulationResult, state_name: str) -> float:
    """Return trapezoidal area under a state-time curve."""

    return float(np.trapezoid(result.state(state_name), result.t))


def tumor_associated(result: SimulationResult) -> np.ndarray:
    """Return total tumor-associated carrier over time."""

    names = ("A_tumor_free", "A_ocn_bound", "A_ysa_bound", "A_internalized")
    return _sum_states(result, names)


def summary_metrics(result: SimulationResult, params: Mapping[str, float]) -> dict[str, float | str]:
    """Summarize one simulation for comparison tables.

    Cumulative process metrics appended below are gross, process-specific
    normalized amounts obtained from the shared postprocessed trapezoidal flux
    integral. MPS metrics use normalized carrier amount; cargo release uses
    normalized payload amount.
    """

    dose = params["D0"]
    tumor = tumor_associated(result)
    summary: dict[str, float | str] = {
        "design": result.design.name,
        "scenario": result.scenario.name,
        "mps_model": result.modes.mps_model,
        "targeting_model": result.modes.targeting_model,
        "solver_method": result.solver_method,
        "success": str(result.success),
        "max_tumor_fraction": float(np.max(tumor) / dose) if dose else 0.0,
        "tumor_auc": float(np.trapezoid(tumor, result.t)),
        "internalized_fraction_final": float(result.state("A_internalized")[-1] / dose)
        if dose
        else 0.0,
        "released_cargo_fraction_final": float(result.state("A_cargo")[-1] / params["Y_cargo"] / dose)
        if dose and params["Y_cargo"]
        else 0.0,
        "elim_deg_fraction_final": float(result.state("A_elim_deg")[-1] / dose) if dose else 0.0,
        "max_abs_carrier_residual": max_abs_carrier_residual(result, params),
        "max_abs_cargo_residual": max_abs_cargo_residual(result, params),
        "minimum_state_value": minimum_state_value(result),
    }
    summary.update(cumulative_flux_metrics(result, params))
    return summary
