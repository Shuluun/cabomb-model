"""Conditional single-batch supply; not cytosolic free or retained total calcium."""

from __future__ import annotations

import math

from cabomb_model.calcium_equivalent import calcium_equivalent_um, required_caco3_mass_pg


def equivalent_from_q(q_pg_per_pl: float, dissolved_fraction: float = 1.0) -> float:
    """Reuse stoichiometry at a 1 pL unit reference, not a typical cell volume."""
    return calcium_equivalent_um(q_pg_per_pl, dissolved_fraction, 1.0)


def required_q(target_um: float, dissolved_fraction: float = 1.0) -> float:
    """Return pure CaCO3 pg per pL required for a conditional equivalent target."""
    return required_caco3_mass_pg(target_um, dissolved_fraction, 1.0)


def _nonnegative_finite(value: float, name: str) -> None:
    if not math.isfinite(value) or value < 0:
        raise ValueError(f"{name} must be finite and nonnegative.")


def normalized_supply(x: float) -> dict[str, float]:
    """Return cumulative fraction, remaining fraction and dF/dx for x=t/t_half."""
    _nonnegative_finite(x, "normalized_time")
    exponent = -math.log(2) * x
    remaining = math.exp(exponent)
    return {
        "x_half_times": x,
        "dissolved_fraction": -math.expm1(exponent),
        "remaining_fraction": remaining,
        "normalized_rate_per_half_time": math.log(2) * remaining,
    }


def single_batch_supply(
    mass0_pg: float,
    volume_pl: float,
    time: float,
    half_time: float,
    *,
    time_unit: str,
    half_time_unit: str,
) -> dict[str, float | str]:
    """Calculate a single time point with explicit s/min/h units and no default rate.

    All pure CaCO3 is already internalized at t=0. Rate units follow time_unit.
    An input half time is a condition, not an inferred or validated sample value.
    """
    seconds = {"s": 1.0, "min": 60.0, "h": 3600.0}
    if time_unit not in seconds or half_time_unit not in seconds:
        raise ValueError("Use explicit time units: s, min or h; normalized x has a separate API.")
    _nonnegative_finite(time, "time")
    if not math.isfinite(half_time) or half_time <= 0:
        raise ValueError("half_time must be finite and positive.")
    complete = calcium_equivalent_um(mass0_pg, 1.0, volume_pl)
    half_in_unit = half_time * (seconds[half_time_unit] / seconds[time_unit])
    if not math.isfinite(half_in_unit) or half_in_unit <= 0:
        raise ValueError("Converted half_time is outside finite numeric range.")
    state = normalized_supply(time / half_in_unit)
    fraction = state["dissolved_fraction"]
    concentration = calcium_equivalent_um(mass0_pg, fraction, volume_pl)
    rate = (complete / half_in_unit) * state["normalized_rate_per_half_time"]
    if not math.isfinite(rate):
        raise ValueError("Supply rate is outside finite numeric range.")
    return {
        "time": time, "time_unit": time_unit,
        "half_time": half_time, "half_time_unit": half_time_unit,
        **state,
        "mass0_pg": mass0_pg,
        "remaining_caco3_pg": mass0_pg * state["remaining_fraction"],
        "dissolved_caco3_pg": mass0_pg * fraction,
        "complete_equivalent_um": complete,
        "cumulative_supply_equivalent_um": concentration,
        "supply_rate_um_per_time_unit": rate,
        "rate_unit": f"uM/{time_unit}",
    }
