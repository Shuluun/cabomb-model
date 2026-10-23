"""Frozen CABOMB V1 ODE right-hand side and rate functions."""

from __future__ import annotations

from collections.abc import Mapping, Sequence

import numpy as np

from cabomb_model.designs import CarrierDesign, get_design
from cabomb_model.fluxes import (
    compute_k_entry,
    degradation_rate,
    instantaneous_fluxes,
)
from cabomb_model.model_modes import ModelModes, V1_1_MODES
from cabomb_model.scenarios import Scenario, make_scenario
from cabomb_model.states import STATE_INDEX


def rhs(
    t: float,
    y: Sequence[float],
    params: Mapping[str, float],
    design: str | CarrierDesign = "ocn-ysa",
    scenario: Scenario | None = None,
    modes: ModelModes = V1_1_MODES,
) -> np.ndarray:
    """Evaluate the frozen ODE right-hand side."""

    del t
    carrier_design = get_design(design)
    run_scenario = scenario or make_scenario(params)
    dy = np.zeros(len(STATE_INDEX), dtype=float)
    state = np.asarray(y, dtype=float)

    a_blood = state[STATE_INDEX["A_blood"]]
    a_tumor_free = state[STATE_INDEX["A_tumor_free"]]
    a_ocn_bound = state[STATE_INDEX["A_ocn_bound"]]
    a_ysa_bound = state[STATE_INDEX["A_ysa_bound"]]
    a_cargo = state[STATE_INDEX["A_cargo"]]

    i_ocn = carrier_design.i_ocn
    i_ysa = carrier_design.i_ysa
    flux = instantaneous_fluxes(state, params, carrier_design, run_scenario, modes)

    dy[STATE_INDEX["A_blood"]] = (
        -flux.liver_mps_uptake
        - flux.spleen_mps_uptake
        + flux.liver_mps_release
        + flux.spleen_mps_release
        - flux.tumor_entry
        - params["k_blood_clear"] * a_blood
        + flux.tumor_return
    )

    dy[STATE_INDEX["A_liver"]] = (
        flux.liver_mps_uptake
        - flux.liver_mps_release
        - flux.liver_mps_degradation
    )

    dy[STATE_INDEX["A_spleen"]] = (
        flux.spleen_mps_uptake
        - flux.spleen_mps_release
        - flux.spleen_mps_degradation
    )

    dy[STATE_INDEX["A_tumor_free"]] = (
        flux.tumor_entry
        - flux.tumor_return
        - flux.ocn_binding
        + i_ocn * params["k_ocn_off"] * a_ocn_bound
        - flux.ysa_binding_from_free
        + i_ysa * params["k_ysa_off"] * a_ysa_bound
        - flux.nonspecific_internalization
        - params["k_tumor_deg"] * a_tumor_free
    )

    dy[STATE_INDEX["A_ocn_bound"]] = (
        flux.ocn_binding
        - i_ocn * params["k_ocn_off"] * a_ocn_bound
        - i_ocn * params["k_ocn_loss"] * a_ocn_bound
        - flux.ocn_to_ysa_transfer
    )

    dy[STATE_INDEX["A_ysa_bound"]] = (
        flux.ysa_binding_from_free
        + flux.ocn_to_ysa_transfer
        - i_ysa * params["k_ysa_off"] * a_ysa_bound
        - flux.ysa_internalization
    )

    dy[STATE_INDEX["A_internalized"]] = (
        flux.internalization - flux.intracellular_degradation
    )

    dy[STATE_INDEX["A_elim_deg"]] = (
        params["k_blood_clear"] * a_blood
        + flux.liver_mps_degradation
        + flux.spleen_mps_degradation
        + params["k_tumor_deg"] * a_tumor_free
        + i_ocn * params["k_ocn_loss"] * a_ocn_bound
        + flux.intracellular_degradation
    )

    dy[STATE_INDEX["A_cargo"]] = (
        flux.cargo_release - params["k_cargo_loss"] * a_cargo
    )

    dy[STATE_INDEX["A_cargo_lost"]] = (
        params["Y_cargo"]
        * (
            params["k_blood_clear"] * a_blood
            + flux.liver_mps_degradation
            + flux.spleen_mps_degradation
            + params["k_tumor_deg"] * a_tumor_free
            + i_ocn * params["k_ocn_loss"] * a_ocn_bound
        )
        + params["k_cargo_loss"] * a_cargo
    )

    return dy


def assert_nonnegative_parameters(params: Mapping[str, float]) -> None:
    """Reject negative executable parameters before integration."""

    negative = sorted(name for name, value in params.items() if float(value) < 0)
    if negative:
        raise ValueError(f"Model parameters must be non-negative: {', '.join(negative)}")

