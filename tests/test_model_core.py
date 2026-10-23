import numpy as np
import pytest

from cabomb_model.designs import DESIGNS
from cabomb_model.metrics import (
    max_abs_cargo_residual,
    max_abs_carrier_residual,
    minimum_state_value,
    summary_metrics,
)
from cabomb_model.model import compute_k_entry, degradation_rate, rhs
from cabomb_model.parameters import default_model_parameters
from cabomb_model.scenarios import make_scenario
from cabomb_model.solver import simulate
from cabomb_model.states import STATE_INDEX, initial_state


@pytest.mark.parametrize("value", [-0.01, 1.01, np.nan, np.inf])
def test_alpha_ocn_cell_access_rejects_out_of_domain_values(value):
    params = default_model_parameters()
    params["alpha_ocn_cell_access"] = value
    with pytest.raises(ValueError, match="alpha_ocn_cell_access"):
        simulate(params)


def test_all_required_designs_are_present():
    assert set(DESIGNS) == {"unmodified", "ocn", "ysa", "ocn-ysa"}
    assert DESIGNS["unmodified"].i_ocn == 0
    assert DESIGNS["unmodified"].i_ysa == 0
    assert DESIGNS["ocn"].i_ocn == 1
    assert DESIGNS["ocn"].i_ysa == 0
    assert DESIGNS["ysa"].i_ocn == 0
    assert DESIGNS["ysa"].i_ysa == 1
    assert DESIGNS["ocn-ysa"].i_ocn == 1
    assert DESIGNS["ocn-ysa"].i_ysa == 1


def test_zero_dose_produces_zero_states():
    params = default_model_parameters()
    params["D0"] = 0.0
    result = simulate(params)
    assert result.success
    assert np.max(np.abs(result.y)) == 0.0


def test_states_are_nonnegative_and_mass_ledgers_close():
    params = default_model_parameters()
    result = simulate(params)
    assert result.success
    assert minimum_state_value(result) >= -1.0e-8
    assert max_abs_carrier_residual(result, params) < 1.0e-7
    assert max_abs_cargo_residual(result, params) < 1.0e-7


def test_removing_ocn_removes_ocn_bound_state():
    params = default_model_parameters()
    result = simulate(params, design="ysa")
    assert np.max(np.abs(result.state("A_ocn_bound"))) < 1.0e-10


def test_removing_ysa_removes_ysa_bound_state_and_ysa_internalization_gain():
    params = default_model_parameters()
    ocn_only = simulate(params, design="ocn")
    combo = simulate(params, design="ocn-ysa")
    assert np.max(np.abs(ocn_only.state("A_ysa_bound"))) < 1.0e-10
    assert combo.state("A_internalized")[-1] > ocn_only.state("A_internalized")[-1]


def test_internalization_rate_does_not_change_tumor_entry_rate():
    params = default_model_parameters()
    scenario = make_scenario(params)
    base_entry = compute_k_entry(params, scenario)
    params_high_int = dict(params)
    params_high_int["k_int"] = params["k_int"] * 5.0
    assert compute_k_entry(params_high_int, scenario) == base_entry


def test_carrier_design_switches_do_not_change_initial_tumor_delivery():
    params = default_model_parameters()
    scenario = make_scenario(params)
    y0 = initial_state(params["D0"])
    initial_delivery = {
        design: rhs(0.0, y0, params, design=design, scenario=scenario)[
            STATE_INDEX["A_tumor_free"]
        ]
        for design in DESIGNS
    }
    assert len(set(np.round(value, 15) for value in initial_delivery.values())) == 1


def test_epr_scenarios_monotonically_change_tumor_entry():
    params = default_model_parameters()
    low = compute_k_entry(params, make_scenario(params, epr="low"))
    medium = compute_k_entry(params, make_scenario(params, epr="medium"))
    high = compute_k_entry(params, make_scenario(params, epr="high"))
    assert low < medium < high


def test_higher_ifp_or_lower_permeability_reduces_entry():
    params = default_model_parameters()
    medium = make_scenario(params, ifp_mmHg=20.0, pcap_tumor_cm_s=params["Pcap_tumor"])
    high_ifp = make_scenario(params, ifp_mmHg=45.0, pcap_tumor_cm_s=params["Pcap_tumor"])
    low_perm = make_scenario(params, ifp_mmHg=20.0, pcap_tumor_cm_s=params["Pcap_tumor"] * 0.5)
    assert compute_k_entry(params, high_ifp) < compute_k_entry(params, medium)
    assert compute_k_entry(params, low_perm) < compute_k_entry(params, medium)


def test_acidic_ph_degrades_faster_than_physiological_ph():
    params = default_model_parameters()
    acidic = degradation_rate(5.5, params)
    physiological = degradation_rate(7.4, params)
    assert acidic > physiological


def test_repeated_deterministic_runs_are_identical():
    params = default_model_parameters()
    first = simulate(params)
    second = simulate(params)
    assert np.array_equal(first.t, second.t)
    assert np.allclose(first.y, second.y, rtol=0.0, atol=1.0e-12)


def test_solver_stability_between_reasonable_configurations():
    params = default_model_parameters()
    primary = simulate(params, method="LSODA", rtol=1.0e-8, atol=1.0e-10)
    comparison = simulate(params, method="RK45", rtol=1.0e-7, atol=1.0e-9)
    primary_summary = summary_metrics(primary, params)
    comparison_summary = summary_metrics(comparison, params)
    assert primary.success
    assert comparison.success
    assert np.array_equal(primary.t, comparison.t)
    assert np.allclose(primary.y, comparison.y, rtol=1.0e-4, atol=1.0e-7)
    assert minimum_state_value(comparison) >= -1.0e-8
    assert max_abs_carrier_residual(comparison, params) < 1.0e-7
    assert max_abs_cargo_residual(comparison, params) < 1.0e-7
    assert np.isclose(
        primary_summary["tumor_auc"],
        comparison_summary["tumor_auc"],
        rtol=5.0e-3,
        atol=1.0e-6,
    )


def test_summary_preserves_v1_metrics_and_adds_cumulative_endpoints():
    params = default_model_parameters()
    result = simulate(params, design="ocn-ysa")

    summary = summary_metrics(result, params)

    assert summary["released_cargo_fraction_final"] == pytest.approx(
        result.state("A_cargo")[-1] / params["Y_cargo"] / params["D0"]
    )
    for name in (
        "cumulative_tumor_entry",
        "cumulative_tumor_return",
        "cumulative_ocn_binding",
        "cumulative_ysa_binding_from_free",
        "cumulative_ocn_to_ysa_transfer",
        "cumulative_internalization",
        "cumulative_intracellular_degradation",
        "cumulative_cargo_release",
    ):
        assert name in summary
        assert float(summary[name]) >= 0.0
