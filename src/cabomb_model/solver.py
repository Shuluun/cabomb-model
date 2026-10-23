"""Numerical integration for the CABOMB V1 ODE model."""

from __future__ import annotations

from dataclasses import dataclass
from typing import Mapping

import numpy as np

from cabomb_model.designs import CarrierDesign, get_design
from cabomb_model.model import assert_nonnegative_parameters, rhs
from cabomb_model.model_modes import ModelModes, V1_1_MODES
from cabomb_model.parameters import (
    default_model_parameters,
    inactive_mode_parameter_symbols,
    validate_executable_parameters,
)
from cabomb_model.scenarios import Scenario, make_scenario
from cabomb_model.states import STATE_INDEX, STATES, initial_state


@dataclass(frozen=True)
class SimulationResult:
    """Container for one ODE simulation."""

    t: np.ndarray
    y: np.ndarray
    state_names: tuple[str, ...]
    design: CarrierDesign
    scenario: Scenario
    modes: ModelModes
    solver_method: str
    success: bool
    message: str
    rtol: float
    atol: float

    def state(self, name: str) -> np.ndarray:
        """Return one state trajectory by name."""

        return self.y[STATE_INDEX[name], :]

    def as_rows(self) -> list[dict[str, float | str]]:
        """Return tidy rows for CSV export."""

        rows: list[dict[str, float | str]] = []
        for col, time_h in enumerate(self.t):
            row: dict[str, float | str] = {
                "time_h": float(time_h),
                "design": self.design.name,
                "scenario": self.scenario.name,
                "solver_method": self.solver_method,
            }
            for state_name in self.state_names:
                row[state_name] = float(self.y[STATE_INDEX[state_name], col])
            rows.append(row)
        return rows


def _rk4_integrate(
    params: Mapping[str, float],
    design: CarrierDesign,
    scenario: Scenario,
    modes: ModelModes,
    t_eval: np.ndarray,
) -> np.ndarray:
    """Fixed-step RK4 fallback used only when SciPy is unavailable."""

    y = np.zeros((len(STATES), len(t_eval)), dtype=float)
    y[:, 0] = np.asarray(initial_state(params["D0"]), dtype=float)
    for idx in range(1, len(t_eval)):
        t_prev = float(t_eval[idx - 1])
        h = float(t_eval[idx] - t_eval[idx - 1])
        y_prev = y[:, idx - 1]
        k1 = rhs(t_prev, y_prev, params, design, scenario, modes)
        k2 = rhs(t_prev + h / 2.0, y_prev + h * k1 / 2.0, params, design, scenario, modes)
        k3 = rhs(t_prev + h / 2.0, y_prev + h * k2 / 2.0, params, design, scenario, modes)
        k4 = rhs(t_prev + h, y_prev + h * k3, params, design, scenario, modes)
        y[:, idx] = y_prev + h * (k1 + 2.0 * k2 + 2.0 * k3 + k4) / 6.0
    return y


def simulate(
    params: Mapping[str, float] | None = None,
    *,
    design: str | CarrierDesign = "ocn-ysa",
    scenario: Scenario | None = None,
    modes: ModelModes | None = None,
    epr: str = "medium",
    mps: str = "medium",
    t_end_h: float = 72.0,
    n_timepoints: int = 241,
    method: str = "LSODA",
    rtol: float = 1.0e-8,
    atol: float = 1.0e-10,
    allow_rk4_fallback: bool = True,
) -> SimulationResult:
    """Run one simulation with SciPy `solve_ivp` when available."""

    run_params = dict(default_model_parameters() if params is None else params)
    run_params.setdefault("alpha_ocn_cell_access", 0.0)
    run_modes = modes or V1_1_MODES
    validate_executable_parameters(run_params, modes=run_modes)
    inactive_symbols = inactive_mode_parameter_symbols(run_modes)
    assert_nonnegative_parameters(
        {
            symbol: value
            for symbol, value in run_params.items()
            if symbol not in inactive_symbols
        }
    )
    run_design = get_design(design)
    run_scenario = scenario or make_scenario(run_params, epr=epr, mps=mps)
    t_eval = np.linspace(0.0, float(t_end_h), int(n_timepoints))
    y0 = np.asarray(initial_state(run_params["D0"]), dtype=float)

    try:
        from scipy.integrate import solve_ivp  # type: ignore
    except ImportError:
        if not allow_rk4_fallback:
            raise
        y = _rk4_integrate(run_params, run_design, run_scenario, run_modes, t_eval)
        return SimulationResult(
            t=t_eval,
            y=y,
            state_names=STATES,
            design=run_design,
            scenario=run_scenario,
            modes=run_modes,
            solver_method="RK4-fallback",
            success=True,
            message="SciPy unavailable; used documented fixed-step RK4 fallback.",
            rtol=rtol,
            atol=atol,
        )

    solution = solve_ivp(
        fun=lambda t, y: rhs(t, y, run_params, run_design, run_scenario, run_modes),
        t_span=(0.0, float(t_end_h)),
        y0=y0,
        method=method,
        t_eval=t_eval,
        rtol=rtol,
        atol=atol,
    )
    return SimulationResult(
        t=t_eval,
        y=np.asarray(solution.y, dtype=float),
        state_names=STATES,
        design=run_design,
        scenario=run_scenario,
        modes=run_modes,
        solver_method=method,
        success=bool(solution.success),
        message=str(solution.message),
        rtol=rtol,
        atol=atol,
    )

