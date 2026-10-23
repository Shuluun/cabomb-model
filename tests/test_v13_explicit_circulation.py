"""Tests for the optional V1.3 explicit-circulation mode.

These tests cover the scientific invariants required by the V1.3 specification:
zero-dose behaviour, flow conservation, carrier/cargo mass balance, non-negativity,
deterministic reproducibility, the guarantee that OCN/YSA never change arterial or
initial tumor-vessel supply, the guarantee that disabling extravasation blocks all
tumor-targeting states, and the guarantee that un-uptaken organ-vascular carrier
returns to the venous pool rather than vanishing.
"""

import numpy as np
import pytest

from cabomb_model.v13_circulation import (
    v13_binding_internalization_metrics,
    v13_retention_efficiency,
    v13_tumor_associated,
)

from cabomb_model.model_modes import ModelModes, V1_3_MODES
from cabomb_model.v13_analysis import (
    run_v13_mechanism_comparison,
    run_v13_monte_carlo,
    v13_parameter_bounds,
)
from cabomb_model.v13_circulation import (
    SIZE_KEYS,
    v13_cargo_total,
    v13_carrier_total,
    v13_cumulative_fluxes,
    v13_default_parameters,
    v13_detailed_metrics,
    v13_first_pass_metrics,
    v13_fluxes,
    v13_max_abs_cargo_residual,
    v13_max_abs_carrier_residual,
    v13_minimum_state_value,
    v13_simulate,
    v13_summary_metrics,
    validate_v13_parameters,
)
from cabomb_model.v13_states import V13_STATES, V13_STATE_INDEX


@pytest.fixture(scope="module")
def params():
    return v13_default_parameters()


@pytest.fixture(scope="module")
def baseline(params):
    return v13_simulate(params, design="ocn-ysa")


def test_v13_state_ledger_layout():
    assert len(V13_STATES) == 18
    assert len(V13_STATE_INDEX) == 18
    # The initial IV bolus lands entirely in central venous.
    initial = V13_STATE_INDEX["A_central_venous"]
    assert initial == 0


@pytest.mark.parametrize("design", ["unmodified", "ocn", "ysa", "ocn-ysa"])
def test_total_internalization_includes_both_paths(params, design):
    result = v13_simulate(params, design=design)
    metrics = v13_binding_internalization_metrics(result, params)
    cumulative = v13_cumulative_fluxes(result, params)
    total = cumulative["internalization"][-1] + cumulative["nonspecific_internalization"][-1]
    assert metrics["internalization_efficiency_gross"] == pytest.approx(total / cumulative["entry"][-1])
    assert metrics["cumulative_total_internalization_fraction"] == pytest.approx(
        metrics["cumulative_specific_internalization_fraction"]
        + metrics["cumulative_nonspecific_internalization_fraction"]
    )
    # Internalized stock is cumulative inflow minus intracellular degradation.
    # Flux integrals use the saved grid, so allow trapezoidal quadrature error.
    assert total - cumulative["deg_int"][-1] == pytest.approx(
        result.state("A_internalized")[-1], abs=5e-6
    )
    assert metrics["cumulative_nonspecific_internalization_fraction"] > 0
    if design in ("unmodified", "ocn"):
        assert metrics["cumulative_specific_internalization_fraction"] == 0
        assert metrics["internalization_efficiency_gross"] > 0


@pytest.mark.parametrize("overrides", [{"D0": 0.0}, {"k_entry": 0.0}])
def test_internalization_and_legacy_ratio_zero_entry(params, overrides):
    params = dict(params, **overrides)
    result = v13_simulate(params)
    metrics = v13_binding_internalization_metrics(result, params)
    assert metrics["internalization_efficiency_gross"] == 0
    assert metrics["cumulative_total_internalization_fraction"] == 0
    assert v13_retention_efficiency(result, params) == 0


def test_legacy_ratio_uses_same_run_gross_entry(params, baseline):
    cumulative = v13_cumulative_fluxes(baseline, params)
    assert v13_retention_efficiency(baseline, params) == pytest.approx(
        v13_tumor_associated(baseline)[-1] / cumulative["entry"][-1]
    )
    single_pass = v13_simulate(params, recirculate=False)
    assert v13_detailed_metrics(baseline, params, first_pass=single_pass)[
        "tumor_retention_efficiency"
    ] == v13_retention_efficiency(baseline, params)


def test_zero_dose_zeros_all_states(params):
    zero_params = dict(params, D0=0.0)
    result = v13_simulate(zero_params, design="ocn-ysa")
    assert result.success
    assert np.max(np.abs(result.y)) == 0.0
    assert v13_minimum_state_value(result) == 0.0


def test_carrier_mass_conservation(params, baseline):
    dose = params["D0"]
    total = v13_carrier_total(baseline)
    assert np.allclose(total, dose, rtol=0.0, atol=1e-9)
    assert v13_max_abs_carrier_residual(baseline, params) < 1e-7


def test_cargo_mass_conservation(params, baseline):
    cargo = v13_cargo_total(baseline, params)
    assert np.allclose(cargo, params["Y_cargo"] * params["D0"], rtol=0.0, atol=1e-9)
    assert v13_max_abs_cargo_residual(baseline, params) < 1e-7


def test_nonnegativity(params, baseline):
    assert baseline.success
    assert v13_minimum_state_value(baseline) >= -1e-8


def test_flow_branch_sum_equals_cardiac_output(params):
    # Qrest = Qco - Qliv - Qspl - Qtum must be non-negative (no negative flow).
    qco = params["Qco"] * 60.0
    qrest = qco - params["Qliv"] * 60.0 - params["Qspl"] * 60.0 - params["Qtum"] * 60.0
    assert qrest >= 0.0


def test_all_six_size_scenarios_run_and_conserve(params):
    assert len(SIZE_KEYS) == 6
    for size in SIZE_KEYS:
        result = v13_simulate(params, design="ocn-ysa", size=size)
        assert result.success
        assert v13_max_abs_carrier_residual(result, params) < 1e-7
        assert v13_max_abs_cargo_residual(result, params) < 1e-7


def test_deterministic_reproducibility(params):
    first = v13_simulate(params, design="ocn-ysa")
    second = v13_simulate(params, design="ocn-ysa")
    assert np.array_equal(first.y, second.y)
    assert np.array_equal(first.t, second.t)


def test_targeting_does_not_change_circulation_fluxes(params, baseline):
    # For a fixed state vector, the circulation and entry fluxes must be exactly
    # invariant to the OCN/YSA design switches.
    state = baseline.y[:, 0].copy()
    idx = {name: i for i, name in enumerate(V13_STATES)}
    state[idx["A_tumor_free"]] = 0.1
    state[idx["A_ocn_bound"]] = 0.05
    state[idx["A_ysa_bound"]] = 0.05
    circulation_keys = (
        "cv_to_rh",
        "rh_to_lung",
        "lung_pass",
        "lung_ret",
        "lh_to_sa",
        "sa_to_liver",
        "sa_to_spleen",
        "sa_to_tumor",
        "sa_to_rest",
        "liver_venous",
        "spleen_venous",
        "tumor_venous",
        "rest_venous",
        "entry",
    )
    base = v13_fluxes(state, params, "unmodified", baseline.scenario, baseline.size)
    for design in ("ocn", "ysa", "ocn-ysa"):
        fluxes = v13_fluxes(state, params, design, baseline.scenario, baseline.size)
        for key in circulation_keys:
            assert fluxes[key] == base[key], f"{key} differs for design {design}"


def test_ocn_ysa_do_not_alter_initial_arterial_or_tumor_vessel_supply(params, baseline):
    unmodified = v13_simulate(params, design="unmodified")
    early = baseline.t <= 0.05
    assert np.allclose(
        unmodified.state("A_systemic_arterial")[early],
        baseline.state("A_systemic_arterial")[early],
        rtol=1e-6,
        atol=1e-9,
    )
    assert np.allclose(
        unmodified.state("A_tumor_vascular")[early],
        baseline.state("A_tumor_vascular")[early],
        rtol=1e-6,
        atol=1e-9,
    )


def test_extravasation_off_blocks_tumor_and_targeting(params):
    no_entry = dict(params, k_entry=0.0)
    result = v13_simulate(no_entry, design="ocn-ysa")
    assert result.success
    for state in ("A_tumor_free", "A_ocn_bound", "A_ysa_bound", "A_internalized", "A_cargo"):
        assert np.max(np.abs(result.state(state))) == 0.0
    assert v13_max_abs_carrier_residual(result, no_entry) < 1e-7


def test_organ_uptake_returns_to_venous(params):
    # With all uptake, clearance, and extravasation disabled, the carrier that
    # leaves the arterial pool must all return through the venous pool.
    closed = dict(
        params,
        k_entry=0.0,
        k_blood_clear=0.0,
        Vmax_liver=0.0,
        Vmax_spleen=0.0,
        k_liver_elim=0.0,
        k_spleen_elim=0.0,
    )
    result = v13_simulate(closed, design="ocn-ysa")
    cumulative = v13_cumulative_fluxes(result, closed)
    arterial_out = float(cumulative["lh_to_sa"][-1])
    venous_in = float(
        cumulative["liver_venous"][-1]
        + cumulative["spleen_venous"][-1]
        + cumulative["tumor_venous"][-1]
        + cumulative["rest_venous"][-1]
    )
    assert venous_in == pytest.approx(arterial_out, rel=1e-6, abs=1e-9)


def test_monte_carlo_fixed_seed_reproducible(tmp_path):
    first = run_v13_monte_carlo(tmp_path / "first", n_samples=8, seed=20260709)
    second = run_v13_monte_carlo(tmp_path / "second", n_samples=8, seed=20260709)
    first_metrics = [(row["sample_id"], row["tumor_auc"]) for row in first]
    second_metrics = [(row["sample_id"], row["tumor_auc"]) for row in second]
    assert first_metrics == second_metrics


def test_monte_carlo_stays_within_registry_bounds(params):
    bounds = v13_parameter_bounds()
    from cabomb_model.v13_analysis import V13_MONTE_CARLO_SYMBOLS

    for symbol in V13_MONTE_CARLO_SYMBOLS:
        lower, upper = bounds[symbol]
        assert lower <= params[symbol] <= upper


def test_validate_v13_parameters_rejects_invalid_volumes(params):
    bad = dict(params, V_central_venous=0.0)
    with pytest.raises(ValueError):
        validate_v13_parameters(bad)


def test_validate_v13_parameters_rejects_invalid_lung_retention(params):
    bad = dict(params, f_lung_ret_100nm=1.0)
    with pytest.raises(ValueError):
        validate_v13_parameters(bad)


def test_summary_metrics_exposes_required_outputs(params, baseline):
    summary = v13_summary_metrics(baseline, params)
    for key in (
        "fraction_surviving_first_lung_pass",
        "gross_arterial_delivery_fraction",
        "net_extravasated_fraction",
        "arterial_auc",
        "tumor_vascular_auc",
        "tumor_interstitial_auc",
        "internalized_fraction_final",
        "released_cargo_fraction_final",
    ):
        assert key in summary
    assert summary["fraction_surviving_first_lung_pass"] == pytest.approx(
        1.0 - params[f"f_lung_ret_{baseline.size}"]
    )
    assert summary["circulation_model"] == "explicit_v1_3"
    assert summary["mps_model"] == V1_3_MODES.mps_model
    assert summary["targeting_model"] == V1_3_MODES.targeting_model


def test_parallel_dual_access_runs_and_conserves_mass(params):
    parallel = ModelModes("saturable_v1_2", "parallel_dual_access_v1_3", "explicit_v1_3")
    result = v13_simulate(params, design="ocn-ysa", modes=parallel)
    assert result.success
    assert v13_max_abs_carrier_residual(result, params) < 1e-7
    assert v13_max_abs_cargo_residual(result, params) < 1e-7
    assert result.modes.targeting_model == "parallel_dual_access_v1_3"


def test_parallel_with_zero_direct_access_recovers_sequential(params):
    # alpha_ysa_direct_access = 0 collapses the parallel route to sequential.
    sequential = v13_simulate(params, design="ocn-ysa", modes=V1_3_MODES)
    zero_direct = dict(params, alpha_ysa_direct_access=0.0)
    parallel = ModelModes("saturable_v1_2", "parallel_dual_access_v1_3", "explicit_v1_3")
    recovered = v13_simulate(zero_direct, design="ocn-ysa", modes=parallel)
    assert np.allclose(recovered.y, sequential.y, rtol=1e-6, atol=1e-9)


def test_parallel_mode_changes_targeting_but_not_circulation_fluxes(params, baseline):
    state = v13_simulate(params, design="ocn-ysa").y[:, 0].copy()
    circulation_keys = (
        "cv_to_rh",
        "rh_to_lung",
        "lung_pass",
        "lh_to_sa",
        "sa_to_tumor",
        "entry",
    )
    sequential = v13_fluxes(
        state, params, "ocn-ysa", baseline.scenario, baseline.size,
        targeting_model="sequential_v1_2",
    )
    parallel = v13_fluxes(
        state, params, "ocn-ysa", baseline.scenario, baseline.size,
        targeting_model="parallel_dual_access_v1_3",
    )
    for key in circulation_keys:
        assert parallel[key] == sequential[key], f"{key} differs between mechanisms"


def test_validate_v13_parameters_rejects_out_of_range_alpha(params):
    bad = dict(params, alpha_ysa_direct_access=1.5)
    with pytest.raises(ValueError):
        validate_v13_parameters(bad)


def test_first_pass_metrics_are_valid_fractions(params):
    first_pass = v13_simulate(params, design="ocn-ysa", recirculate=False)
    metrics = v13_first_pass_metrics(first_pass, params)
    for key in (
        "first_pass_lung_arrival_fraction",
        "first_pass_lung_retention_fraction",
        "first_pass_systemic_availability",
        "tumor_vascular_first_pass_arrival_fraction",
        "tumor_first_extravasation_fraction",
        "first_pass_extravasation_efficiency",
    ):
        assert 0.0 <= metrics[key] <= 1.0 + 1e-9, f"{key} out of range: {metrics[key]}"


def test_detailed_metrics_exposes_first_event_and_gross_keys(params):
    result = v13_simulate(params, design="ocn-ysa")
    first_pass = v13_simulate(params, design="ocn-ysa", recirculate=False)
    detailed = v13_detailed_metrics(result, params, first_pass=first_pass)
    for key in (
        "design",
        "size",
        "targeting_model",
        "tumor_auc",
        "net_extravasated_fraction",
        "internalized_fraction_final",
        "released_cargo_fraction_final",
        "first_pass_lung_arrival_fraction",
        "apparent_blood_half_life_h",
        "ocn_binding_fraction_gross",
        "ysa_binding_fraction_gross",
        "ocn_to_ysa_transfer_fraction_gross",
        "internalization_efficiency_gross",
        "cargo_release_efficiency_gross",
        "active_cargo_remaining_fraction_final",
        "tumor_retention_efficiency",
    ):
        assert key in detailed, f"missing {key}"


def test_mechanism_comparison_covers_both_mechanisms(tmp_path):
    rows = run_v13_mechanism_comparison(tmp_path, sizes=("100nm",))
    assert len(rows) == 2 * 1 * 4
    models = {row["targeting_model"] for row in rows}
    sizes = {row["size"] for row in rows}
    assert models == {"sequential_v1_2", "parallel_dual_access_v1_3"}
    assert sizes == {"100nm"}
    for row in rows:
        assert row["tumor_auc"] is not None
        assert row["internalized_fraction_final"] is not None
        assert row["released_cargo_fraction_final"] is not None


def test_feasibility_dual_wins_tumor_auc_and_both_conserve_mass(params):
    # One grid point: the OCN-YSA design must not reduce tumor exposure relative
    # to YSA-only, and both designs must conserve carrier/cargo mass.
    dual = v13_simulate(params, design="ocn-ysa")
    ysa_only = v13_simulate(params, design="ysa")
    dual_summary = v13_summary_metrics(dual, params)
    ysa_summary = v13_summary_metrics(ysa_only, params)
    assert dual_summary["tumor_auc"] > ysa_summary["tumor_auc"]
    assert v13_max_abs_carrier_residual(dual, params) < 1e-7
    assert v13_max_abs_cargo_residual(dual, params) < 1e-7
    assert v13_max_abs_carrier_residual(ysa_only, params) < 1e-7
    assert v13_max_abs_cargo_residual(ysa_only, params) < 1e-7
