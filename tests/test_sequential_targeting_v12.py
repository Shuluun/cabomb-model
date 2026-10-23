import numpy as np
import pytest

from cabomb_model.designs import DESIGNS
from cabomb_model.fluxes import cumulative_flux_metrics, instantaneous_fluxes
from cabomb_model.model import compute_k_entry, rhs
from cabomb_model.model_modes import V1_1_MODES, V1_2_MODES
from cabomb_model.parameters import default_model_parameters
from cabomb_model.scenarios import make_scenario
from cabomb_model.solver import simulate
from cabomb_model.states import CARRIER_STATES, STATE_INDEX, STATES


def test_v12_combo_has_no_direct_free_ysa_binding_but_ysa_only_does():
    params = default_model_parameters()
    scenario = make_scenario(params)
    y = np.zeros(len(STATES))
    y[STATE_INDEX["A_tumor_free"]] = 0.2

    combo = instantaneous_fluxes(y, params, "ocn-ysa", scenario, V1_2_MODES)
    ysa = instantaneous_fluxes(y, params, "ysa", scenario, V1_2_MODES)

    assert combo.ysa_binding_from_free == 0.0
    assert ysa.ysa_binding_from_free > 0.0


def test_v12_k_cross_requires_both_design_switches():
    params = default_model_parameters()
    scenario = make_scenario(params)
    y = np.zeros(len(STATES))
    y[STATE_INDEX["A_ocn_bound"]] = 0.2

    values = {
        design: instantaneous_fluxes(
            y, params, design, scenario, V1_2_MODES
        ).ocn_to_ysa_transfer
        for design in DESIGNS
    }

    assert values["ocn-ysa"] > 0.0
    assert values["unmodified"] == values["ocn"] == values["ysa"] == 0.0


def test_alpha_and_k_cross_are_mode_isolated():
    params = default_model_parameters()
    v11_base = simulate(params, modes=V1_1_MODES)
    v11_changed_cross = simulate(dict(params, k_cross=0.6), modes=V1_1_MODES)
    assert np.allclose(v11_base.y, v11_changed_cross.y, rtol=0.0, atol=0.0)

    v12_base = simulate(params, modes=V1_2_MODES)
    v12_changed_alpha = simulate(
        dict(params, alpha_ocn_cell_access=1.0), modes=V1_2_MODES
    )
    assert np.allclose(v12_base.y, v12_changed_alpha.y, rtol=0.0, atol=0.0)


def test_sequential_transfer_is_equal_and_opposite():
    params = default_model_parameters()
    scenario = make_scenario(params)
    y = np.zeros(len(STATES))
    y[STATE_INDEX["A_ocn_bound"]] = 0.2
    zero = dict(params, k_cross=0.0)

    delta = rhs(1.0, y, params, "ocn-ysa", scenario, V1_2_MODES) - rhs(
        1.0, y, zero, "ocn-ysa", scenario, V1_2_MODES
    )

    assert delta[STATE_INDEX["A_ocn_bound"]] < 0.0
    assert delta[STATE_INDEX["A_ysa_bound"]] == pytest.approx(
        -delta[STATE_INDEX["A_ocn_bound"]]
    )
    assert sum(delta[STATE_INDEX[name]] for name in CARRIER_STATES) == pytest.approx(
        0.0
    )


def test_targeting_modes_do_not_change_entry_formula_at_same_blood_state():
    params = default_model_parameters()
    scenario = make_scenario(params)

    assert compute_k_entry(params, scenario) == compute_k_entry(
        dict(params, alpha_ocn_cell_access=1.0, k_cross=0.6), scenario
    )


def test_k_cross_zero_and_positive_cumulative_cases():
    params = default_model_parameters()
    zero = dict(params, k_cross=0.0)

    assert (
        cumulative_flux_metrics(simulate(zero, modes=V1_2_MODES), zero)[
            "cumulative_ocn_to_ysa_transfer"
        ]
        == 0.0
    )
    assert (
        cumulative_flux_metrics(simulate(params, modes=V1_2_MODES), params)[
            "cumulative_ocn_to_ysa_transfer"
        ]
        > 0.0
    )
