"""Physical accounting and explicit time-unit contracts, without sample assumptions."""

import math

import numpy as np
import pytest

from cabomb_model.calcium_equivalent import calcium_equivalent_um
from cabomb_model.calcium_supply_time import equivalent_from_q, normalized_supply, required_q, single_batch_supply


@pytest.mark.parametrize("target,q", [(10, 0.00100087), (100, 0.0100087), (1000, 0.100087)])
def test_independent_budget_and_inverse(target, q):
    assert required_q(target) == pytest.approx(q, rel=1e-12)
    for fraction in (0.1, 0.5, 1):
        assert required_q(target, fraction) == pytest.approx(q / fraction)
        assert equivalent_from_q(required_q(target, fraction), fraction) == pytest.approx(target)


@pytest.mark.parametrize("mass,volume", [(0.1, 2), (0.2, 4), (0.05, 1), (1, 20)])
def test_same_q_does_not_require_a_cell_size(mass, volume):
    assert calcium_equivalent_um(mass, 0.5, volume) == pytest.approx(equivalent_from_q(0.05, 0.5))


@pytest.mark.parametrize("x,fraction", [(0, 0), (1, 0.5), (2, 0.75), (3, 0.875), (4, 0.9375)])
def test_half_time_checkpoints(x, fraction):
    result = normalized_supply(x)
    assert result["dissolved_fraction"] == pytest.approx(fraction, abs=1e-15)
    assert result["remaining_fraction"] + result["dissolved_fraction"] == pytest.approx(1)
    assert result["normalized_rate_per_half_time"] == pytest.approx(math.log(2) * (1-fraction))


def test_zero_tiny_and_long_time_monotonicity():
    results = [normalized_supply(float(x)) for x in np.linspace(0, 100, 501)]
    fractions = [r["dissolved_fraction"] for r in results]
    assert fractions[0] == 0
    assert np.all(np.diff(fractions) >= 0)
    assert all(0 <= f <= 1 for f in fractions)
    tiny = normalized_supply(1e-15)["dissolved_fraction"]
    assert tiny == pytest.approx(math.log(2)*1e-15, rel=1e-14, abs=0)
    assert normalized_supply(1e308)["remaining_fraction"] == 0


@pytest.mark.parametrize("time", [0, 0.1, 1, 2, 4, 100])
def test_mass_ledger_and_static_agreement(time):
    # Purely mathematical unit test, not a selected biological half-time scenario.
    state = single_batch_supply(0.1, 2, time, 1, time_unit="h", half_time_unit="h")
    assert state["remaining_caco3_pg"] + state["dissolved_caco3_pg"] == pytest.approx(0.1)
    assert state["cumulative_supply_equivalent_um"] == pytest.approx(calcium_equivalent_um(0.1, state["dissolved_fraction"], 2))
    assert state["cumulative_supply_equivalent_um"] <= state["complete_equivalent_um"]


def test_time_units_and_rate_units_are_consistent():
    hour = single_batch_supply(0.1, 2, 1, 1, time_unit="h", half_time_unit="h")
    sec = single_batch_supply(0.1, 2, 3600, 1, time_unit="s", half_time_unit="h")
    minute = single_batch_supply(0.1, 2, 60, 3600, time_unit="min", half_time_unit="s")
    assert hour["dissolved_fraction"] == sec["dissolved_fraction"] == minute["dissolved_fraction"] == 0.5
    assert sec["supply_rate_um_per_time_unit"] * 3600 == pytest.approx(hour["supply_rate_um_per_time_unit"])
    assert minute["supply_rate_um_per_time_unit"] * 60 == pytest.approx(hour["supply_rate_um_per_time_unit"])
    assert sec["rate_unit"] == "uM/s"
    epsilon = 1e-5
    upper = single_batch_supply(0.1, 2, 1+epsilon, 1, time_unit="h", half_time_unit="h")
    lower = single_batch_supply(0.1, 2, 1-epsilon, 1, time_unit="h", half_time_unit="h")
    derivative = (upper["cumulative_supply_equivalent_um"]-lower["cumulative_supply_equivalent_um"])/(2*epsilon)
    assert derivative == pytest.approx(hour["supply_rate_um_per_time_unit"], rel=1e-9)


def test_zero_mass_and_zero_static_dissolution():
    state = single_batch_supply(0, 2, 1, 2, time_unit="min", half_time_unit="min")
    assert state["cumulative_supply_equivalent_um"] == state["supply_rate_um_per_time_unit"] == 0
    assert equivalent_from_q(1, 0) == 0
    assert equivalent_from_q(0, 1) == 0
    with pytest.raises(ValueError):
        required_q(10, 0)


@pytest.mark.parametrize("value", [-1, float("nan"), float("inf")])
def test_invalid_normalized_time_and_q(value):
    with pytest.raises(ValueError):
        normalized_supply(value)
    with pytest.raises(ValueError):
        equivalent_from_q(value)


@pytest.mark.parametrize("half", [0, -1, float("nan"), float("inf")])
def test_invalid_half_time(half):
    with pytest.raises(ValueError):
        single_batch_supply(1, 1, 0, half, time_unit="h", half_time_unit="h")


@pytest.mark.parametrize("kwargs", [dict(time=-1), dict(time=float("inf")), dict(volume_pl=0),
                                    dict(mass0_pg=-1), dict(time_unit="normalized"), dict(half_time_unit="hours")])
def test_invalid_batch_inputs(kwargs):
    args = dict(mass0_pg=1, volume_pl=1, time=1, half_time=1, time_unit="h", half_time_unit="h")
    args.update(kwargs)
    with pytest.raises(ValueError):
        single_batch_supply(**args)
