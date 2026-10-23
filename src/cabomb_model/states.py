"""State definitions for the frozen CABOMB V1 ODE system."""

from __future__ import annotations

STATES: tuple[str, ...] = (
    "A_blood",
    "A_liver",
    "A_spleen",
    "A_tumor_free",
    "A_ocn_bound",
    "A_ysa_bound",
    "A_internalized",
    "A_elim_deg",
    "A_cargo",
    "A_cargo_lost",
)

CARRIER_STATES: tuple[str, ...] = (
    "A_blood",
    "A_liver",
    "A_spleen",
    "A_tumor_free",
    "A_ocn_bound",
    "A_ysa_bound",
    "A_internalized",
    "A_elim_deg",
)

ACTIVE_CARRIER_STATES: tuple[str, ...] = (
    "A_blood",
    "A_liver",
    "A_spleen",
    "A_tumor_free",
    "A_ocn_bound",
    "A_ysa_bound",
    "A_internalized",
)

TUMOR_ASSOCIATED_STATES: tuple[str, ...] = (
    "A_tumor_free",
    "A_ocn_bound",
    "A_ysa_bound",
    "A_internalized",
)

STATE_INDEX: dict[str, int] = {name: idx for idx, name in enumerate(STATES)}


def zeros_state() -> list[float]:
    """Return a zero-filled state vector in canonical order."""

    return [0.0 for _ in STATES]


def initial_state(dose: float) -> list[float]:
    """Return the IV bolus initial condition."""

    state = zeros_state()
    state[STATE_INDEX["A_blood"]] = float(dose)
    return state

