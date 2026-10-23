from dataclasses import replace

import numpy as np
import pytest

from cabomb_model.fluxes import (
    cumulative_flux_metrics,
    cumulative_flux_trajectories,
    instantaneous_fluxes,
)
from cabomb_model.metrics import max_abs_cargo_residual, max_abs_carrier_residual
from cabomb_model.model import rhs
from cabomb_model.model_modes import V1_1_MODES, V1_2_MODES
from cabomb_model.parameters import default_model_parameters
from cabomb_model.scenarios import make_scenario
from cabomb_model.solver import simulate
from cabomb_model.states import STATE_INDEX, STATES


MPS_CUMULATIVE_NAMES = {
    "cumulative_liver_mps_uptake",
    "cumulative_spleen_mps_uptake",
    "cumulative_liver_mps_release",
    "cumulative_spleen_mps_release",
    "cumulative_liver_mps_degradation",
    "cumulative_spleen_mps_degradation",
}


def flux_at_blood(params, blood, scenario=None):
    y = np.zeros(len(STATES))
    y[STATE_INDEX["A_blood"]] = blood
    return instantaneous_fluxes(
        y, params, "ocn-ysa", scenario or make_scenario(params), V1_2_MODES
    )


def test_saturable_uptake_is_zero_monotone_and_bounded():
    params = default_model_parameters()
    values = [flux_at_blood(params, blood) for blood in (0.0, 0.1, 1.0, 1000.0)]
    liver = [value.liver_mps_uptake for value in values]
    spleen = [value.spleen_mps_uptake for value in values]
    assert liver[0] == spleen[0] == 0.0
    assert liver == sorted(liver)
    assert spleen == sorted(spleen)
    assert liver[-1] <= params["Vmax_liver"]
    assert spleen[-1] <= params["Vmax_spleen"]


def test_low_concentration_slope_matches_documented_linear_mapping():
    params = default_model_parameters()
    blood = 1.0e-9 * params["Vblood"]
    flux = flux_at_blood(params, blood)
    assert flux.liver_mps_uptake == pytest.approx(
        params["k_liver_uptake"] * blood, rel=1.0e-8
    )
    assert flux.spleen_mps_uptake == pytest.approx(
        params["k_spleen_uptake"] * blood, rel=1.0e-8
    )


def test_release_and_degradation_flux_destinations_close_both_ledgers():
    params = default_model_parameters()
    result = simulate(params, modes=V1_2_MODES)
    assert result.success
    assert max_abs_carrier_residual(result, params) < 1.0e-7
    assert max_abs_cargo_residual(result, params) < 1.0e-7


def test_linear_mode_remains_exact_v11_even_with_new_parameters():
    params = default_model_parameters()
    base = simulate(params, modes=V1_1_MODES)
    changed = dict(
        params,
        Vmax_liver=999.0,
        Vmax_spleen=999.0,
        Km_liver=999.0,
        Km_spleen=999.0,
        k_liver_release=9.0,
        k_spleen_release=9.0,
    )
    inactive = simulate(changed, modes=V1_1_MODES)
    assert np.allclose(base.y, inactive.y, rtol=0.0, atol=0.0)


@pytest.mark.parametrize(
    ("organ_state", "release_parameter"),
    (
        ("A_liver", "k_liver_release"),
        ("A_spleen", "k_spleen_release"),
    ),
)
def test_mps_release_returns_intact_carrier_without_terminal_ledger_flux(
    organ_state, release_parameter
):
    params = default_model_parameters()
    params[release_parameter] = 0.25
    params["k_liver_elim"] = 0.0
    params["k_spleen_elim"] = 0.0
    y = np.zeros(len(STATES))
    y[STATE_INDEX[organ_state]] = 0.4

    derivative = rhs(0.0, y, params, modes=V1_2_MODES)
    release = params[release_parameter] * y[STATE_INDEX[organ_state]]

    assert derivative[STATE_INDEX["A_blood"]] == pytest.approx(release)
    assert derivative[STATE_INDEX[organ_state]] == pytest.approx(-release)
    assert derivative[STATE_INDEX["A_elim_deg"]] == 0.0
    assert derivative[STATE_INDEX["A_cargo_lost"]] == 0.0


@pytest.mark.parametrize(
    ("organ_state", "elim_parameter", "release_parameter"),
    (
        ("A_liver", "k_liver_elim", "k_liver_release"),
        ("A_spleen", "k_spleen_elim", "k_spleen_release"),
    ),
)
def test_mps_degradation_has_terminal_ledger_destination_without_blood_return(
    organ_state, elim_parameter, release_parameter
):
    params = default_model_parameters()
    params[release_parameter] = 0.0
    params[elim_parameter] = 0.25
    y = np.zeros(len(STATES))
    y[STATE_INDEX[organ_state]] = 0.4

    derivative = rhs(0.0, y, params, modes=V1_2_MODES)
    degradation = params[elim_parameter] * y[STATE_INDEX[organ_state]]

    assert derivative[STATE_INDEX["A_blood"]] == 0.0
    assert derivative[STATE_INDEX[organ_state]] == pytest.approx(-degradation)
    assert derivative[STATE_INDEX["A_elim_deg"]] == pytest.approx(degradation)
    assert derivative[STATE_INDEX["A_cargo_lost"]] == pytest.approx(
        params["Y_cargo"] * degradation
    )


def test_saturable_mps_scenario_factor_scales_only_uptake_and_high_limit():
    params = default_model_parameters()
    baseline_scenario = make_scenario(params)
    scaled_scenario = replace(baseline_scenario, name="scaled-mps", mps_factor=1.7)
    blood = 0.5
    baseline = flux_at_blood(params, blood, baseline_scenario)
    scaled = flux_at_blood(params, blood, scaled_scenario)

    assert scaled.liver_mps_uptake == pytest.approx(1.7 * baseline.liver_mps_uptake)
    assert scaled.spleen_mps_uptake == pytest.approx(1.7 * baseline.spleen_mps_uptake)

    y = np.zeros(len(STATES))
    y[STATE_INDEX["A_liver"]] = 0.4
    y[STATE_INDEX["A_spleen"]] = 0.3
    baseline_organ = instantaneous_fluxes(
        y, params, "ocn-ysa", baseline_scenario, V1_2_MODES
    )
    scaled_organ = instantaneous_fluxes(
        y, params, "ocn-ysa", scaled_scenario, V1_2_MODES
    )
    assert scaled_organ.liver_mps_release == baseline_organ.liver_mps_release
    assert scaled_organ.spleen_mps_release == baseline_organ.spleen_mps_release
    assert (
        scaled_organ.liver_mps_degradation
        == baseline_organ.liver_mps_degradation
    )
    assert (
        scaled_organ.spleen_mps_degradation
        == baseline_organ.spleen_mps_degradation
    )

    high = flux_at_blood(params, 1.0e15 * params["Vblood"], scaled_scenario)
    assert high.liver_mps_uptake == pytest.approx(
        scaled_scenario.mps_factor * params["Vmax_liver"], rel=1.0e-12
    )
    assert high.spleen_mps_uptake == pytest.approx(
        scaled_scenario.mps_factor * params["Vmax_spleen"], rel=1.0e-12
    )


def test_v12_mps_cumulative_fluxes_are_present_and_monotone():
    params = default_model_parameters()
    result = simulate(params, modes=V1_2_MODES)

    cumulative = cumulative_flux_trajectories(result, params)

    assert MPS_CUMULATIVE_NAMES <= set(cumulative)
    for name in MPS_CUMULATIVE_NAMES:
        assert cumulative[name][0] == 0.0
        assert np.min(np.diff(cumulative[name])) >= -1.0e-12


def test_v12_mps_endpoints_converge_on_saved_time_grid():
    params = default_model_parameters()
    endpoints = {
        points: cumulative_flux_metrics(
            simulate(params, modes=V1_2_MODES, n_timepoints=points), params
        )
        for points in (241, 481, 961)
    }

    for name in MPS_CUMULATIVE_NAMES:
        coarse_gap = abs(endpoints[481][name] - endpoints[241][name])
        fine_gap = abs(endpoints[961][name] - endpoints[481][name])
        assert fine_gap <= coarse_gap + 1.0e-12
        assert endpoints[481][name] == pytest.approx(
            endpoints[961][name], rel=1.0e-3, abs=1.0e-8
        )


def test_v11_mps_cumulative_release_is_zero():
    params = default_model_parameters()
    metrics = cumulative_flux_metrics(simulate(params, modes=V1_1_MODES), params)

    assert metrics["cumulative_liver_mps_release"] == 0.0
    assert metrics["cumulative_spleen_mps_release"] == 0.0
