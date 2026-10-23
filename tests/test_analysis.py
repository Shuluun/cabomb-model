import csv
from dataclasses import replace
from pathlib import Path
import re

import numpy as np
import pytest

import cabomb_model.analysis as analysis
from cabomb_model.analysis import (
    MONTE_CARLO_SYMBOLS,
    TARGET_METRICS,
    run_all_analyses,
    run_baseline,
    run_local_sensitivity,
    run_monte_carlo,
    run_ocn_accessibility_sweep,
    summarize_paired_design_advantage,
)
from cabomb_model.designs import DESIGNS
from cabomb_model.metrics import summary_metrics
from cabomb_model.model_modes import ModelModes, V1_1_MODES, V1_2_MODES
from cabomb_model.parameters import default_model_parameters
from cabomb_model.solver import simulate


SUMMARY_TABLES = (
    "baseline_summary.csv",
    "design_comparison.csv",
    "scenario_comparison.csv",
    "monte_carlo.csv",
    "ocn_accessibility_sweep.csv",
)


def test_baseline_writes_expected_outputs(tmp_path: Path):
    outputs = run_baseline(tmp_path)
    assert outputs["trajectory_path"].exists()
    assert outputs["summary_path"].exists()
    assert outputs["figure_path"].exists()


def test_monte_carlo_is_reproducible(tmp_path: Path):
    first = run_monte_carlo(tmp_path / "a", n_samples=8, seed=123)
    second = run_monte_carlo(tmp_path / "b", n_samples=8, seed=123)
    assert first == second
    for filename in (
        "monte_carlo.csv",
        "monte_carlo_design_advantage_summary.csv",
    ):
        assert (tmp_path / "a" / "tables" / filename).read_bytes() == (
            tmp_path / "b" / "tables" / filename
        ).read_bytes()


def test_monte_carlo_pairs_all_designs_with_identical_parameters(tmp_path: Path):
    rows = run_monte_carlo(tmp_path, n_samples=6, seed=123)
    assert len(rows) == 24
    for sample_id in range(6):
        sample = [row for row in rows if int(row["sample_id"]) == sample_id]
        assert [row["design"] for row in sample] == list(DESIGNS)
        for symbol in MONTE_CARLO_SYMBOLS:
            assert len({float(row[symbol]) for row in sample}) == 1
    alpha_values = {
        float(row["alpha_ocn_cell_access"])
        for row in rows
        if row["design"] == "ocn-ysa"
    }
    assert all(0.0 <= value <= 1.0 for value in alpha_values)
    assert any(0.0 < value < 1.0 for value in alpha_values)


def test_paired_summary_uses_strict_wins_and_bounded_fractions():
    rows = [
        {
            "sample_id": 0,
            "design": "unmodified",
            "tumor_auc": 0.0,
            "cumulative_cargo_release": 0.0,
        },
        {
            "sample_id": 0,
            "design": "ocn",
            "tumor_auc": 0.0,
            "cumulative_cargo_release": 0.0,
        },
        {
            "sample_id": 0,
            "design": "ysa",
            "tumor_auc": 1.0,
            "cumulative_cargo_release": 2.0,
        },
        {
            "sample_id": 0,
            "design": "ocn-ysa",
            "tumor_auc": 2.0,
            "cumulative_cargo_release": 3.0,
        },
        {
            "sample_id": 1,
            "design": "unmodified",
            "tumor_auc": 0.0,
            "cumulative_cargo_release": 0.0,
        },
        {
            "sample_id": 1,
            "design": "ocn",
            "tumor_auc": 0.0,
            "cumulative_cargo_release": 0.0,
        },
        {
            "sample_id": 1,
            "design": "ysa",
            "tumor_auc": 2.0,
            "cumulative_cargo_release": 2.0,
        },
        {
            "sample_id": 1,
            "design": "ocn-ysa",
            "tumor_auc": 2.0,
            "cumulative_cargo_release": 1.0,
        },
    ]
    for row in rows:
        row.update(
            mps_model="linear_v1", targeting_model="accessibility_v1_1"
        )
    summary = summarize_paired_design_advantage(rows)
    assert summary["sample_count"] == 2
    assert summary["cargo_release_win_count"] == 1
    assert summary["tumor_auc_win_count"] == 1
    assert summary["both_win_count"] == 1
    assert summary["cargo_release_win_fraction"] == 0.5
    assert summary["tumor_auc_win_fraction"] == 0.5
    assert summary["both_win_fraction"] == 0.5
    assert summary["mps_model"] == "linear_v1"
    assert summary["targeting_model"] == "accessibility_v1_1"
    assert "not probabilities of clinical success" in str(summary["interpretation"])


def test_paired_summary_rejects_duplicate_or_incomplete_samples():
    complete = [
        {
            "sample_id": 0,
            "design": design,
            "tumor_auc": 1.0,
            "cumulative_cargo_release": 1.0,
            "mps_model": "linear_v1",
            "targeting_model": "accessibility_v1_1",
        }
        for design in DESIGNS
    ]
    with pytest.raises(ValueError, match="Duplicate design"):
        summarize_paired_design_advantage([*complete, dict(complete[-1])])
    with pytest.raises(ValueError, match="all four designs"):
        summarize_paired_design_advantage(complete[:-1])
    with pytest.raises(ValueError, match="at least one sample"):
        summarize_paired_design_advantage([])


def test_paired_summary_rejects_nonfinite_comparison_metrics():
    rows = [
        {
            "sample_id": 0,
            "design": design,
            "tumor_auc": np.nan if design == "ocn-ysa" else 1.0,
            "cumulative_cargo_release": 1.0,
            "mps_model": "linear_v1",
            "targeting_model": "accessibility_v1_1",
        }
        for design in DESIGNS
    ]

    with pytest.raises(ValueError, match="finite comparison metrics"):
        summarize_paired_design_advantage(rows)


def test_local_sensitivity_covers_three_targets_and_alpha(tmp_path: Path):
    rows = run_local_sensitivity(tmp_path)
    assert {row["target_metric"] for row in rows} == set(TARGET_METRICS)
    alpha_rows = [
        row for row in rows if row["parameter_symbol"] == "alpha_ocn_cell_access"
    ]
    assert len(alpha_rows) == 3
    assert all(row["sensitivity_method"] == "one_sided_boundary" for row in alpha_rows)
    assert all(float(row["low_value"]) == 0.0 for row in alpha_rows)
    assert all(float(row["high_value"]) == 0.1 for row in alpha_rows)


def test_local_sensitivity_uses_generic_metric_columns(tmp_path: Path):
    row = run_local_sensitivity(tmp_path)[0]
    assert {
        "baseline_metric",
        "low_metric",
        "high_metric",
        "normalized_sensitivity",
        "abs_normalized_sensitivity",
    } <= set(row)


@pytest.mark.parametrize("target_metric", TARGET_METRICS)
def test_local_sensitivity_maps_each_requested_target_metric(
    tmp_path: Path, target_metric: str
):
    params = default_model_parameters()
    expected = float(summary_metrics(simulate(params), params)[target_metric])
    rows = run_local_sensitivity(tmp_path, target_metrics=(target_metric,))
    assert len(rows) == 13
    assert all(row["target_metric"] == target_metric for row in rows)
    assert all(float(row["baseline_metric"]) == pytest.approx(expected) for row in rows)
    for row in rows:
        for column in (
            "baseline_metric",
            "low_metric",
            "high_metric",
            "finite_difference_slope",
            "normalized_sensitivity",
            "abs_normalized_sensitivity",
        ):
            value = row[column]
            assert value is None or np.isfinite(float(value))


def test_local_sensitivity_is_deterministic_and_alpha_changes_combo_targets(
    tmp_path: Path,
):
    targets = ("cumulative_internalization", "cumulative_cargo_release")
    first = run_local_sensitivity(tmp_path / "a", target_metrics=targets)
    second = run_local_sensitivity(tmp_path / "b", target_metrics=targets)
    assert first == second
    alpha_rows = [
        row for row in first if row["parameter_symbol"] == "alpha_ocn_cell_access"
    ]
    assert all(float(row["high_metric"]) != float(row["low_metric"]) for row in alpha_rows)


def test_full_analysis_writes_tables(tmp_path: Path):
    run_all_analyses(tmp_path)
    expected = [
        "baseline_timecourse.csv",
        "baseline_summary.csv",
        "design_comparison.csv",
        "scenario_comparison.csv",
        "local_sensitivity.csv",
        "monte_carlo.csv",
        "monte_carlo_design_advantage_summary.csv",
        "ocn_accessibility_sweep.csv",
    ]
    for filename in expected:
        assert (tmp_path / "tables" / filename).exists()
    assert (tmp_path / "figures" / "parameters_master_preview.png").exists()
    assert (tmp_path / "figures" / "cumulative_cargo_release_vs_alpha.svg").exists()
    assert (
        tmp_path / "figures" / "tumor_auc_vs_cumulative_cargo_release.svg"
    ).exists()


def test_ocn_accessibility_figures_explain_scientific_interpretation(tmp_path: Path):
    run_ocn_accessibility_sweep(tmp_path)
    release_svg = (
        tmp_path / "figures" / "cumulative_cargo_release_vs_alpha.svg"
    ).read_text(encoding="utf-8")
    tradeoff_svg = (
        tmp_path / "figures" / "tumor_auc_vs_cumulative_cargo_release.svg"
    ).read_text(encoding="utf-8")

    assert "effective YSA-EphA2 access retained by OCN-bound carrier" in release_svg
    assert "may shift carrier toward internalization and release" in release_svg
    assert "benefit is not guaranteed" in release_svg
    assert "Tumor AUC measures exposure/retention" in tradeoff_svg
    assert "cumulative release measures cell-processing throughput" in tradeoff_svg
    assert "not efficacy or clinical outcomes" in tradeoff_svg
    for svg in (release_svg, tradeoff_svg):
        assert all(line == line.rstrip(" \t") for line in svg.splitlines())

    def text_y(svg: str, label: str) -> float:
        match = re.search(
            rf"<!-- {re.escape(label)} -->\s*"
            r'<g transform="translate\([^ ]+ ([0-9.]+)\)',
            svg,
        )
        assert match is not None
        return float(match.group(1))

    release_xlabel_y = text_y(
        release_svg,
        "OCN-bound cellular accessibility, alpha (dimensionless)",
    )
    release_caption_y = text_y(
        release_svg,
        "Alpha is the effective YSA-EphA2 access retained by OCN-bound carrier.",
    )
    tradeoff_xlabel_y = text_y(
        tradeoff_svg,
        "Tumor-associated carrier AUC (normalized carrier h)",
    )
    tradeoff_caption_y = text_y(
        tradeoff_svg,
        "Tumor AUC measures exposure/retention; cumulative release measures cell-processing throughput.",
    )
    assert release_caption_y - release_xlabel_y >= 18.0
    assert tradeoff_caption_y - tradeoff_xlabel_y >= 18.0


def test_ocn_accessibility_sweep_has_twenty_rows_and_only_combo_changes(tmp_path: Path):
    rows = run_ocn_accessibility_sweep(tmp_path)
    assert len(rows) == 20
    assert {float(row["alpha_ocn_cell_access"]) for row in rows} == {
        0.0,
        0.25,
        0.5,
        0.75,
        1.0,
    }
    assert {str(row["design"]) for row in rows} == {
        "unmodified",
        "ocn",
        "ysa",
        "ocn-ysa",
    }
    for design in ("unmodified", "ocn", "ysa"):
        selected = [row for row in rows if row["design"] == design]
        for metric in (
            "tumor_auc",
            "cumulative_internalization",
            "cumulative_cargo_release",
            "cumulative_ocn_to_ysa_transfer",
        ):
            values = np.asarray([float(row[metric]) for row in selected])
            assert np.allclose(values, values[0], rtol=0.0, atol=1e-12), metric

    combo = sorted(
        (row for row in rows if row["design"] == "ocn-ysa"),
        key=lambda row: float(row["alpha_ocn_cell_access"]),
    )
    assert float(combo[0]["cumulative_ocn_to_ysa_transfer"]) == 0.0
    assert float(combo[-1]["cumulative_ocn_to_ysa_transfer"]) > 0.0
    assert not np.allclose(
        [float(row["cumulative_cargo_release"]) for row in combo],
        float(combo[0]["cumulative_cargo_release"]),
        rtol=0.0,
        atol=1e-12,
    )


def test_ocn_accessibility_sweep_writes_required_table_and_figures(tmp_path: Path):
    rows = run_ocn_accessibility_sweep(tmp_path)
    table = tmp_path / "tables" / "ocn_accessibility_sweep.csv"
    release_figure = tmp_path / "figures" / "cumulative_cargo_release_vs_alpha.svg"
    tradeoff_figure = (
        tmp_path / "figures" / "tumor_auc_vs_cumulative_cargo_release.svg"
    )
    assert table.exists()
    assert release_figure.exists()
    assert tradeoff_figure.exists()

    with table.open(newline="", encoding="utf-8") as handle:
        saved_rows = list(csv.DictReader(handle))
    assert len(saved_rows) == len(rows)
    assert {
        "alpha_ocn_cell_access",
        "design",
        "scenario",
        "max_tumor_fraction",
        "tumor_auc",
        "cumulative_internalization",
        "cumulative_cargo_release",
        "cumulative_ocn_to_ysa_transfer",
        "success",
        "max_abs_carrier_residual",
        "max_abs_cargo_residual",
        "minimum_state_value",
    }.issubset(saved_rows[0])

    release_svg = release_figure.read_text(encoding="utf-8")
    tradeoff_svg = tradeoff_figure.read_text(encoding="utf-8")
    assert not any(line.endswith((" ", "\t")) for line in release_svg.splitlines())
    assert not any(line.endswith((" ", "\t")) for line in tradeoff_svg.splitlines())
    assert "Cumulative cargo release (normalized cargo)" in release_svg
    assert "Tumor-associated carrier AUC (normalized carrier h)" in tradeoff_svg
    assert "Exploratory normalized outputs; not clinical predictions." in release_svg
    assert "alpha=1" in tradeoff_svg


def test_full_analysis_summary_tables_satisfy_invariants(tmp_path: Path):
    run_all_analyses(tmp_path)
    for filename in SUMMARY_TABLES:
        with (tmp_path / "tables" / filename).open(newline="", encoding="utf-8") as handle:
            rows = list(csv.DictReader(handle))
        assert rows, filename
        for row in rows:
            assert row["success"] == "True"
            assert float(row["minimum_state_value"]) >= analysis.NONNEGATIVE_TOLERANCE
            assert (
                float(row["max_abs_carrier_residual"])
                <= analysis.MASS_RESIDUAL_TOLERANCE
            )
            assert (
                float(row["max_abs_cargo_residual"])
                <= analysis.MASS_RESIDUAL_TOLERANCE
            )


def test_baseline_rejects_failed_solver_result(monkeypatch, tmp_path: Path):
    params = default_model_parameters()
    failed = replace(simulate(params), success=False, message="forced failure")

    def fake_simulate(*args, **kwargs):
        return failed

    monkeypatch.setattr(analysis, "simulate", fake_simulate)
    with pytest.raises(RuntimeError, match="solver failed"):
        run_baseline(tmp_path)


def test_summary_and_design_output_record_modes(tmp_path: Path):
    rows = analysis.run_design_comparison(tmp_path, modes=V1_2_MODES)
    assert len(rows) == 4
    assert {row["mps_model"] for row in rows} == {"saturable_v1_2"}
    assert {row["targeting_model"] for row in rows} == {"sequential_v1_2"}


def test_v12_monte_carlo_pairs_modes_and_systemic_parameters(tmp_path: Path):
    rows = analysis.run_monte_carlo(
        tmp_path, n_samples=4, seed=123, modes=V1_2_MODES
    )
    assert len(rows) == 16
    for sample_id in range(4):
        sample = [row for row in rows if int(row["sample_id"]) == sample_id]
        assert {row["mps_model"] for row in sample} == {"saturable_v1_2"}
        assert {row["targeting_model"] for row in sample} == {"sequential_v1_2"}
        for symbol in analysis.V12_MONTE_CARLO_SYMBOLS:
            assert len({float(row[symbol]) for row in sample}) == 1


def test_alpha_sweep_remains_v11_specific(tmp_path: Path):
    with pytest.raises(ValueError, match="accessibility_v1_1"):
        analysis.run_ocn_accessibility_sweep(tmp_path, modes=V1_2_MODES)


def test_default_analysis_modes_remain_v11(tmp_path: Path):
    rows = analysis.run_design_comparison(tmp_path, modes=V1_1_MODES)
    assert {row["mps_model"] for row in rows} == {"linear_v1"}
    assert {row["targeting_model"] for row in rows} == {"accessibility_v1_1"}


def test_v12_local_sensitivity_uses_only_v12_parameter_set(tmp_path: Path):
    rows = analysis.run_local_sensitivity(
        tmp_path, target_metrics=("tumor_auc",), modes=V1_2_MODES
    )
    assert {row["parameter_symbol"] for row in rows} == set(
        analysis.V12_SENSITIVITY_SYMBOLS
    )
    assert {row["mps_model"] for row in rows} == {"saturable_v1_2"}
    assert {row["targeting_model"] for row in rows} == {"sequential_v1_2"}


@pytest.mark.parametrize(
    ("modes", "included", "excluded"),
    (
        (
            ModelModes("saturable_v1_2", "accessibility_v1_1"),
            {"Vmax_liver", "Km_spleen", "k_liver_release", "alpha_ocn_cell_access"},
            {"k_cross", "k_liver_uptake", "k_spleen_uptake"},
        ),
        (
            ModelModes("linear_v1", "sequential_v1_2"),
            {"k_liver_uptake", "k_spleen_uptake", "k_cross"},
            {"alpha_ocn_cell_access", "Vmax_liver", "Km_liver"},
        ),
    ),
)
def test_mixed_modes_use_only_active_parameter_sets(
    tmp_path: Path,
    modes: ModelModes,
    included: set[str],
    excluded: set[str],
):
    rows = analysis.run_local_sensitivity(
        tmp_path, target_metrics=("tumor_auc",), modes=modes
    )
    symbols = {str(row["parameter_symbol"]) for row in rows}
    assert included <= symbols
    assert symbols.isdisjoint(excluded)


def test_paired_advantage_summary_records_and_requires_one_mode_pair():
    rows = [
        {
            "sample_id": 0,
            "design": design,
            "tumor_auc": 1.0,
            "cumulative_cargo_release": 1.0,
            "mps_model": "saturable_v1_2",
            "targeting_model": "sequential_v1_2",
        }
        for design in DESIGNS
    ]
    summary = summarize_paired_design_advantage(rows)
    assert summary["mps_model"] == "saturable_v1_2"
    assert summary["targeting_model"] == "sequential_v1_2"

    mixed = [dict(row) for row in rows]
    mixed[-1]["targeting_model"] = "accessibility_v1_1"
    with pytest.raises(ValueError, match="one consistent mode pair"):
        summarize_paired_design_advantage(mixed)
