"""Shared instantaneous flux definitions for CABOMB V1.1."""

from __future__ import annotations

from collections.abc import Mapping, Sequence
from dataclasses import dataclass
from typing import TYPE_CHECKING

import numpy as np

from cabomb_model.designs import CarrierDesign, get_design
from cabomb_model.model_modes import ModelModes, V1_1_MODES
from cabomb_model.scenarios import Scenario
from cabomb_model.states import STATE_INDEX

if TYPE_CHECKING:
    from cabomb_model.solver import SimulationResult


FLUX_ATTRIBUTE_MAP: dict[str, str] = {
    "liver_mps_uptake": "liver_mps_uptake",
    "spleen_mps_uptake": "spleen_mps_uptake",
    "liver_mps_release": "liver_mps_release",
    "spleen_mps_release": "spleen_mps_release",
    "liver_mps_degradation": "liver_mps_degradation",
    "spleen_mps_degradation": "spleen_mps_degradation",
    "tumor_entry": "tumor_entry",
    "tumor_return": "tumor_return",
    "ocn_binding": "ocn_binding",
    "ysa_binding_from_free": "ysa_binding_from_free",
    "ocn_to_ysa_transfer": "ocn_to_ysa_transfer",
    "internalization": "internalization",
    "ysa_internalization": "ysa_internalization",
    "nonspecific_internalization": "nonspecific_internalization",
    "intracellular_degradation": "intracellular_degradation",
    "cargo_release": "cargo_release",
}

CUMULATIVE_NAME_MAP: dict[str, str] = {
    name: f"cumulative_{name}" for name in FLUX_ATTRIBUTE_MAP
}

_NEGATIVE_FLUX_TOLERANCE = 1.0e-12
_NEGATIVE_STATE_TOLERANCE = 1.0e-8


@dataclass(frozen=True)
class InstantaneousFluxes:
    """Tracked non-negative process fluxes for one state vector."""

    liver_mps_uptake: float
    spleen_mps_uptake: float
    liver_mps_release: float
    spleen_mps_release: float
    liver_mps_degradation: float
    spleen_mps_degradation: float
    tumor_entry: float
    tumor_return: float
    ocn_binding: float
    ysa_binding_from_free: float
    ocn_to_ysa_transfer: float
    ysa_internalization: float
    nonspecific_internalization: float
    intracellular_degradation: float
    cargo_release: float

    @property
    def internalization(self) -> float:
        """Return total cellular internalization flux."""

        return self.ysa_internalization + self.nonspecific_internalization


def degradation_rate(pH: float, params: Mapping[str, float]) -> float:
    """Return acid-dependent intracellular degradation rate in h^-1."""

    pH_phys = params["pH_blood_phys"]
    pH_acid = params["pH_endosome_lysosome"]
    k_phys = params["k_deg_phys"]
    k_acid = params["k_deg_int_acid"]
    if pH <= pH_acid:
        return k_acid
    if pH >= pH_phys:
        return k_phys
    fraction_acid = (pH_phys - pH) / (pH_phys - pH_acid)
    return k_phys + fraction_acid * (k_acid - k_phys)


def compute_k_entry(params: Mapping[str, float], scenario: Scenario) -> float:
    """Compute the effective blood-to-tumor entry rate in h^-1."""

    pcap_factor = scenario.pcap_tumor_cm_s / params["Pcap_tumor_ref"]
    ifp_factor = (1.0 + params["ifp_entry_sensitivity"] * params["IFP_ref"]) / (
        1.0 + params["ifp_entry_sensitivity"] * scenario.ifp_mmHg
    )
    return params["k_entry"] * scenario.epr_factor * pcap_factor * ifp_factor


def instantaneous_fluxes(
    y: Sequence[float],
    params: Mapping[str, float],
    design: str | CarrierDesign,
    scenario: Scenario,
    modes: ModelModes = V1_1_MODES,
) -> InstantaneousFluxes:
    """Evaluate tracked process fluxes for one state and scenario."""

    state = np.asarray(y, dtype=float)
    carrier_design = get_design(design)
    i_ocn = carrier_design.i_ocn
    i_ysa = carrier_design.i_ysa
    blood = float(state[STATE_INDEX["A_blood"]])
    liver = float(state[STATE_INDEX["A_liver"]])
    spleen = float(state[STATE_INDEX["A_spleen"]])
    tumor_free = float(state[STATE_INDEX["A_tumor_free"]])
    ocn_bound = float(state[STATE_INDEX["A_ocn_bound"]])
    ysa_bound = float(state[STATE_INDEX["A_ysa_bound"]])
    internalized = float(state[STATE_INDEX["A_internalized"]])
    ysa_rate = params["k_ysa_on"] * params["R_EphA2_eff"]
    degradation = degradation_rate(scenario.intracellular_pH, params)
    degradation_flux = degradation * internalized

    if modes.mps_model == "saturable_v1_2":
        concentration = blood / params["Vblood"]
        liver_uptake = (
            scenario.mps_factor
            * params["Vmax_liver"]
            * concentration
            / (params["Km_liver"] + concentration)
        )
        spleen_uptake = (
            scenario.mps_factor
            * params["Vmax_spleen"]
            * concentration
            / (params["Km_spleen"] + concentration)
        )
        liver_release = params["k_liver_release"] * liver
        spleen_release = params["k_spleen_release"] * spleen
    else:
        liver_uptake = scenario.mps_factor * params["k_liver_uptake"] * blood
        spleen_uptake = scenario.mps_factor * params["k_spleen_uptake"] * blood
        liver_release = 0.0
        spleen_release = 0.0

    if modes.targeting_model == "accessibility_v1_1":
        ysa_binding_from_free = i_ysa * ysa_rate * tumor_free
        ocn_to_ysa_transfer = (
            i_ocn
            * i_ysa
            * float(params.get("alpha_ocn_cell_access", 0.0))
            * ysa_rate
            * ocn_bound
        )
    else:
        direct_ysa_allowed = i_ysa * (1 - i_ocn)
        ysa_binding_from_free = direct_ysa_allowed * ysa_rate * tumor_free
        ocn_to_ysa_transfer = (
            i_ocn
            * i_ysa
            * params["k_cross"]
            * params["R_EphA2_eff"]
            * ocn_bound
        )

    return InstantaneousFluxes(
        liver_mps_uptake=liver_uptake,
        spleen_mps_uptake=spleen_uptake,
        liver_mps_release=liver_release,
        spleen_mps_release=spleen_release,
        liver_mps_degradation=params["k_liver_elim"] * liver,
        spleen_mps_degradation=params["k_spleen_elim"] * spleen,
        tumor_entry=compute_k_entry(params, scenario) * blood,
        tumor_return=params["k_return"] * tumor_free,
        ocn_binding=i_ocn * params["k_ocn_on"] * tumor_free,
        ysa_binding_from_free=ysa_binding_from_free,
        ocn_to_ysa_transfer=ocn_to_ysa_transfer,
        ysa_internalization=i_ysa * params["k_int"] * ysa_bound,
        nonspecific_internalization=params["k_int_ns"] * tumor_free,
        intracellular_degradation=degradation_flux,
        cargo_release=params["Y_cargo"] * degradation_flux,
    )


def _validate_result_grid(result: "SimulationResult") -> None:
    """Validate the saved trajectory grid used for postprocessing."""

    if result.t.ndim != 1 or result.y.ndim != 2:
        raise ValueError("Simulation time and state trajectories must be arrays")
    if result.y.shape[1] != result.t.size:
        raise ValueError("Simulation time and state trajectory lengths must match")
    if result.t.size == 0 or not np.all(np.isfinite(result.t)):
        raise ValueError("Simulation time axis must be finite and nonempty")
    if result.t.size > 1 and not np.all(np.diff(result.t) > 0.0):
        raise ValueError("Simulation time axis must be strictly increasing")


def flux_trajectories(
    result: "SimulationResult", params: Mapping[str, float]
) -> dict[str, np.ndarray]:
    """Evaluate tracked process rates on a saved simulation trajectory.

    Tiny negative rates caused by accepted state-solver roundoff are clipped to
    zero. A larger negative rate is rejected rather than silently converting an
    invalid simulation into a monotone cumulative trajectory.
    """

    _validate_result_grid(result)
    values = {
        name: np.zeros_like(result.t, dtype=float) for name in FLUX_ATTRIBUTE_MAP
    }
    for index in range(result.t.size):
        state = result.y[:, index]
        if not np.all(np.isfinite(state)):
            raise ValueError(f"State is non-finite at saved index {index}")
        minimum_state = float(np.min(state))
        if minimum_state < -_NEGATIVE_STATE_TOLERANCE:
            raise ValueError(
                "A state is negative beyond numerical tolerance "
                f"at saved index {index}: {minimum_state}"
            )
        postprocess_state = np.maximum(state, 0.0)
        flux = instantaneous_fluxes(
            postprocess_state, params, result.design, result.scenario, result.modes
        )
        for name, attribute in FLUX_ATTRIBUTE_MAP.items():
            value = float(getattr(flux, attribute))
            if not np.isfinite(value):
                raise ValueError(f"Flux {name} is non-finite at saved index {index}")
            if value < -_NEGATIVE_FLUX_TOLERANCE:
                raise ValueError(
                    f"Flux {name} is negative beyond numerical tolerance "
                    f"at saved index {index}: {value}"
                )
            values[name][index] = max(0.0, value)
    return values


def _cumulative_trapezoid(values: np.ndarray, time: np.ndarray) -> np.ndarray:
    """Return composite-trapezoid cumulative integral with a zero origin."""

    increments = 0.5 * (values[1:] + values[:-1]) * np.diff(time)
    return np.concatenate((np.array([0.0]), np.cumsum(increments)))


def cumulative_flux_trajectories(
    result: "SimulationResult", params: Mapping[str, float]
) -> dict[str, np.ndarray]:
    """Integrate tracked process rates after solving the ten-state ODE.

    Values are gross cumulative, process-specific normalized amounts,
    calculated only by trapezoidal postprocessing on the saved simulation time
    grid (hours). Carrier pathways, including all MPS metrics, use normalized
    carrier amount; ``cargo_release`` uses normalized payload amount.
    """

    fluxes = flux_trajectories(result, params)
    return {
        CUMULATIVE_NAME_MAP[name]: _cumulative_trapezoid(values, result.t)
        for name, values in fluxes.items()
    }


def cumulative_flux_metrics(
    result: "SimulationResult", params: Mapping[str, float]
) -> dict[str, float]:
    """Return final gross cumulative process amounts for one simulation."""

    return {
        name: float(values[-1])
        for name, values in cumulative_flux_trajectories(result, params).items()
    }
