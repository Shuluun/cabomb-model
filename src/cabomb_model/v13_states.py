"""State definitions for the optional V1.3 explicit-circulation mode.

This module is intentionally separate from ``states.py`` so the frozen
V1.1/V1.2 ten-state path is never touched. V1.3 keeps the same tissue/organ
states and replaces the single lumped ``A_blood`` with nine vascular states.
"""

from __future__ import annotations

V13_STATES: tuple[str, ...] = (
    "A_central_venous",
    "A_right_heart",
    "A_lung_vascular",
    "A_left_heart",
    "A_systemic_arterial",
    "A_liver_vascular",
    "A_spleen_vascular",
    "A_tumor_vascular",
    "A_rest_vascular",
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

V13_STATE_INDEX: dict[str, int] = {name: index for index, name in enumerate(V13_STATES)}

V13_VASCULAR_STATES: tuple[str, ...] = V13_STATES[:9]

V13_CARRIER_STATES: tuple[str, ...] = V13_STATES[:16]

V13_ACTIVE_CARRIER_STATES: tuple[str, ...] = V13_STATES[:15]

V13_TUMOR_ASSOCIATED_STATES: tuple[str, ...] = (
    "A_tumor_free",
    "A_ocn_bound",
    "A_ysa_bound",
    "A_internalized",
)


def v13_initial_state(dose: float) -> list[float]:
    """Return the V1.3 initial state with the whole IV bolus in central venous."""

    state = [0.0] * len(V13_STATES)
    state[V13_STATE_INDEX["A_central_venous"]] = float(dose)
    return state
