from dataclasses import replace
from pathlib import Path

import numpy as np
import pandas as pd
import pytest

from cabomb_model.designs import DESIGNS
from cabomb_model.fluxes import (
    cumulative_flux_metrics,
    cumulative_flux_trajectories,
    flux_trajectories,
    instantaneous_fluxes,
)
from cabomb_model.metrics import (
    max_abs_cargo_residual,
    max_abs_carrier_residual,
    minimum_state_value,
    summary_metrics,
)
from cabomb_model.model_modes import V1_1_MODES
from cabomb_model.parameters import default_model_parameters
from cabomb_model.scenarios import make_scenario
from cabomb_model.solver import simulate
from cabomb_model.states import CARRIER_STATES, STATE_INDEX, STATES


def test_ocn_to_ysa_flux_requires_both_switches_and_positive_alpha():
    params = default_model_parameters()
    params["alpha_ocn_cell_access"] = 0.5
    scenario = make_scenario(params)
    y = np.zeros(len(STATES))
    y[STATE_INDEX["A_ocn_bound"]] = 0.2

    values = {
        design: instantaneous_fluxes(y, params, design, scenario).ocn_to_ysa_transfer
        for design in DESIGNS
    }

    assert values["ocn-ysa"] > 0.0
    assert values["unmodified"] == 0.0
    assert values["ocn"] == 0.0
    assert values["ysa"] == 0.0


def test_alpha_does_not_change_tumor_entry_flux():
    params = default_model_parameters()
    scenario = make_scenario(params)
    y = np.zeros(len(STATES))
    y[STATE_INDEX["A_blood"]] = 1.0

    entries = []
    for alpha in (0.0, 0.5, 1.0):
        run = dict(params, alpha_ocn_cell_access=alpha)
        entries.append(instantaneous_fluxes(y, run, "ocn-ysa", scenario).tumor_entry)

    assert entries[0] == entries[1] == entries[2]


def test_alpha_zero_reproduces_v1_state_trajectories_and_metrics():
    params = default_model_parameters()
    params["alpha_ocn_cell_access"] = 0.0
    trajectory_reference = pd.read_csv(
        Path("tests/fixtures/v1_state_trajectories_01dba30.csv")
    )
    summary_reference = pd.read_csv(
        Path("tests/fixtures/v1_summary_metrics_01dba30.csv")
    ).set_index("design")

    for design in DESIGNS:
        result = simulate(params, design=design, modes=V1_1_MODES)
        expected = trajectory_reference.loc[
            trajectory_reference["design"] == design, list(result.state_names)
        ].to_numpy(dtype=float).T
        assert result.y.shape[0] == 10
        assert np.allclose(result.y, expected, rtol=1e-8, atol=1e-10)
        actual = summary_metrics(result, params)
        for metric in (
            "max_tumor_fraction",
            "tumor_auc",
            "internalized_fraction_final",
            "released_cargo_fraction_final",
            "elim_deg_fraction_final",
        ):
            assert actual[metric] == pytest.approx(
                float(summary_reference.loc[design, metric]), rel=1e-8, abs=1e-10
            )


def test_new_transfer_is_equal_and_opposite_in_carrier_equations():
    from cabomb_model.model import rhs

    params = default_model_parameters()
    params["alpha_ocn_cell_access"] = 1.0
    scenario = make_scenario(params)
    y = np.zeros(len(STATES))
    y[STATE_INDEX["A_ocn_bound"]] = 0.2
    base = dict(params, alpha_ocn_cell_access=0.0)

    delta = rhs(1.0, y, params, "ocn-ysa", scenario) - rhs(
        1.0, y, base, "ocn-ysa", scenario
    )

    assert delta[STATE_INDEX["A_ocn_bound"]] < 0.0
    assert delta[STATE_INDEX["A_ysa_bound"]] == pytest.approx(
        -delta[STATE_INDEX["A_ocn_bound"]]
    )
    carrier_delta = sum(delta[STATE_INDEX[name]] for name in CARRIER_STATES)
    assert carrier_delta == pytest.approx(0.0, abs=1e-15)


def test_positive_alpha_only_changes_combination_trajectory():
    params = default_model_parameters()
    for design in DESIGNS:
        alpha_zero = simulate(dict(params, alpha_ocn_cell_access=0.0), design=design)
        alpha_one_params = dict(params, alpha_ocn_cell_access=1.0)
        alpha_one = simulate(alpha_one_params, design=design)
        if design == "ocn-ysa":
            assert not np.allclose(alpha_zero.y, alpha_one.y, rtol=0.0, atol=1e-12)
            assert np.max(alpha_one.state("A_ysa_bound")) > np.max(
                alpha_zero.state("A_ysa_bound")
            )
        else:
            assert np.allclose(alpha_zero.y, alpha_one.y, rtol=0.0, atol=1e-12)
        assert minimum_state_value(alpha_one) >= -1.0e-8
        assert max_abs_carrier_residual(alpha_one, alpha_one_params) < 1.0e-7
        assert max_abs_cargo_residual(alpha_one, alpha_one_params) < 1.0e-7


def test_cumulative_fluxes_start_at_zero_and_are_monotone():
    params = default_model_parameters()
    params["alpha_ocn_cell_access"] = 0.75
    result = simulate(params, design="ocn-ysa")

    cumulative = cumulative_flux_trajectories(result, params)

    assert set(cumulative) == {
        "cumulative_liver_mps_uptake",
        "cumulative_spleen_mps_uptake",
        "cumulative_liver_mps_release",
        "cumulative_spleen_mps_release",
        "cumulative_liver_mps_degradation",
        "cumulative_spleen_mps_degradation",
        "cumulative_tumor_entry",
        "cumulative_tumor_return",
        "cumulative_ocn_binding",
        "cumulative_ysa_binding_from_free",
        "cumulative_ocn_to_ysa_transfer",
        "cumulative_internalization",
        "cumulative_ysa_internalization",
        "cumulative_nonspecific_internalization",
        "cumulative_intracellular_degradation",
        "cumulative_cargo_release",
    }
    for name, values in cumulative.items():
        assert values.shape == result.t.shape, name
        assert values[0] == 0.0, name
        assert np.min(np.diff(values)) >= -1.0e-12, name


def test_internalization_components_and_cargo_release_identity():
    params = default_model_parameters()
    result = simulate(params, design="ocn-ysa")

    cumulative = cumulative_flux_trajectories(result, params)

    assert np.allclose(
        cumulative["cumulative_internalization"],
        cumulative["cumulative_ysa_internalization"]
        + cumulative["cumulative_nonspecific_internalization"],
        rtol=0.0,
        atol=1.0e-12,
    )
    assert np.allclose(
        cumulative["cumulative_cargo_release"],
        params["Y_cargo"]
        * cumulative["cumulative_intracellular_degradation"],
        rtol=1.0e-12,
        atol=1.0e-12,
    )


def test_alpha_transfer_cumulative_metric_has_expected_zero_and_positive_cases():
    params = default_model_parameters()
    zero_params = dict(params, alpha_ocn_cell_access=0.0)
    positive_params = dict(params, alpha_ocn_cell_access=0.5)
    zero = simulate(zero_params, design="ocn-ysa")
    positive = simulate(positive_params, design="ocn-ysa")

    assert (
        cumulative_flux_metrics(zero, zero_params)[
            "cumulative_ocn_to_ysa_transfer"
        ]
        == 0.0
    )
    assert (
        cumulative_flux_metrics(positive, positive_params)[
            "cumulative_ocn_to_ysa_transfer"
        ]
        > 0.0
    )


def test_flux_trajectories_clip_only_tiny_negative_state_roundoff():
    params = default_model_parameters()
    result = simulate(params, design="ocn-ysa")
    perturbed_y = result.y.copy()
    perturbed_y[STATE_INDEX["A_tumor_free"], 0] = -1.0e-14
    perturbed = replace(result, y=perturbed_y)

    fluxes = flux_trajectories(perturbed, params)

    assert all(values[0] >= 0.0 for values in fluxes.values())


def test_flux_trajectories_use_accepted_state_tolerance_for_roundoff():
    params = default_model_parameters()
    params["alpha_ocn_cell_access"] = 1.0
    result = simulate(params, design="ocn-ysa")
    perturbed_y = result.y.copy()
    perturbed_y[STATE_INDEX["A_ocn_bound"], 1] = -5.0e-9
    accepted = replace(result, y=perturbed_y)

    fluxes = flux_trajectories(accepted, params)

    assert fluxes["ocn_to_ysa_transfer"][1] == 0.0
    assert accepted.y[STATE_INDEX["A_ocn_bound"], 1] == -5.0e-9


def test_flux_trajectories_reject_state_below_numerical_tolerance():
    params = default_model_parameters()
    params["alpha_ocn_cell_access"] = 1.0
    result = simulate(params, design="ocn-ysa")
    perturbed_y = result.y.copy()
    perturbed_y[STATE_INDEX["A_ocn_bound"], 1] = -2.0e-8
    invalid = replace(result, y=perturbed_y)

    with pytest.raises(ValueError, match="state is negative beyond numerical tolerance"):
        flux_trajectories(invalid, params)


def test_flux_trajectories_do_not_hide_unrelated_negative_rate():
    params = default_model_parameters()
    params["alpha_ocn_cell_access"] = 1.0
    result = simulate(params, design="ocn-ysa")
    perturbed_y = result.y[:, 1:2].copy()
    perturbed_y[STATE_INDEX["A_ocn_bound"], 0] = -5.0e-9
    accepted_roundoff = replace(result, t=result.t[1:2], y=perturbed_y)
    invalid_params = dict(params, k_return=-1.0)

    with pytest.raises(ValueError, match="Flux tumor_return is negative"):
        flux_trajectories(accepted_roundoff, invalid_params)


def test_cumulative_fluxes_reject_nonincreasing_time_axis():
    params = default_model_parameters()
    result = simulate(params, design="ocn-ysa")
    invalid = replace(result, t=result.t.copy())
    invalid.t[2] = invalid.t[1]

    with pytest.raises(ValueError, match="strictly increasing"):
        cumulative_flux_trajectories(invalid, params)


@pytest.mark.parametrize("alpha", [0.0, 1.0])
def test_cumulative_metrics_converge_with_saved_time_grid(alpha):
    params = default_model_parameters()
    params["alpha_ocn_cell_access"] = alpha
    endpoints = {}
    for points in (241, 481, 961):
        result = simulate(params, design="ocn-ysa", n_timepoints=points)
        endpoints[points] = cumulative_flux_metrics(result, params)

    for metric in endpoints[961]:
        coarse_gap = abs(endpoints[481][metric] - endpoints[241][metric])
        fine_gap = abs(endpoints[961][metric] - endpoints[481][metric])
        assert fine_gap <= coarse_gap + 1.0e-12, metric
        assert endpoints[481][metric] == pytest.approx(
            endpoints[961][metric], rel=1.0e-3, abs=1.0e-8
        ), metric
