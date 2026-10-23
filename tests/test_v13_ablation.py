"""Integration checks for the scoped comparison and its saved evidence."""

import csv
import hashlib
import json
from pathlib import Path

import numpy as np
import pytest

from cabomb_model.v13_analysis import run_v13_design_comparison
from cabomb_model.v13_circulation import v13_default_parameters, v13_simulate, v13_detailed_metrics


def read_csv(path):
    with path.open(encoding="utf-8", newline="") as handle:
        return list(csv.DictReader(handle))


def test_ablation_common_inputs_metrics_and_provenance(tmp_path):
    rows = run_v13_design_comparison(tmp_path)
    assert [row["design"] for row in rows] == ["unmodified", "ocn", "ysa", "ocn-ysa"]
    manifest = json.loads((tmp_path / "run_manifest.json").read_text(encoding="utf-8"))
    assert manifest["size"] == "100nm"
    assert manifest["modes"]["targeting"] == "sequential_v1_2"
    assert manifest["design_switches"]["unmodified"] == {"i_ocn": 0, "i_ysa": 0}
    grids = []
    for row in rows:
        name = row["design"]
        trajectory = read_csv(tmp_path / "tables" / f"{name}_timecourse.csv")
        cumulative = read_csv(tmp_path / "tables" / f"{name}_cumulative_fluxes.csv")
        grids.append([float(point["time_h"]) for point in trajectory])
        assert float(trajectory[0]["A_central_venous"]) == manifest["parameters"]["D0"]
        assert row["time_final_h"] == 72
        assert row["max_abs_carrier_residual"] < 1e-7
        assert row["max_abs_cargo_residual"] < 1e-7
        assert row["minimum_state_value"] >= -1e-8
        dose = manifest["parameters"]["D0"]
        expected = (float(cumulative[-1]["internalization"]) + float(cumulative[-1]["nonspecific_internalization"])) / dose
        assert row["cumulative_total_internalization_fraction"] == pytest.approx(expected)
        assert float(cumulative[-1]["deg_int"]) / dose == pytest.approx(row["cargo_release_efficiency_gross"])
    assert all(grid == grids[0] for grid in grids)
    assert manifest["time_grid_h"] == grids[0]
    for name, digest in manifest["output_sha256"].items():
        assert hashlib.sha256((tmp_path / name).read_bytes()).hexdigest() == digest
    summary = (tmp_path / "summary_zh.md").read_text(encoding="utf-8")
    assert "相对变化" in summary and "累计非特异内化" in summary
    assert "实验重复或置信区间" in summary
    # Compare postprocessed CSVs with direct integration on the same backend.
    # Frozen historical replay is reported separately by check_solver_agreement.py.
    reference = v13_simulate(v13_default_parameters(), design="ocn-ysa", size="100nm").as_rows()
    current = read_csv(tmp_path / "tables/ocn-ysa_timecourse.csv")
    assert len(reference) == len(current)
    for key in reference[0]:
        if key == "time_h" or key.startswith("A_"):
            assert np.allclose([float(r[key]) for r in reference], [float(r[key]) for r in current], rtol=1e-9, atol=1e-11)
        else:
            assert [str(r[key]) for r in reference] == [r[key] for r in current]


def test_v13_report_total_components_and_legacy_retention_exclusion(tmp_path, monkeypatch):
    import cabomb_model.v13_reporting as reporting

    params = v13_default_parameters()
    result = v13_simulate(params, design="unmodified")
    single = v13_simulate(params, design="unmodified", recirculate=False)
    detailed = v13_detailed_metrics(result, params, first_pass=single)
    payload = {name: [] for name in ("design_comparison", "size_comparison", "scenario_comparison", "local_sensitivity", "monte_carlo", "mass_balance_audit")}
    payload.update(baseline={"result": result, "params": params}, baseline_detailed=[detailed])
    monkeypatch.setattr(reporting, "run_v13_all_analyses", lambda output_dir: payload)
    path = reporting.build_v13_markdown_report(tmp_path, tmp_path / "report.md")
    text = path.read_text(encoding="utf-8")
    assert "Tumor retention efficiency" not in text
    total_line = next(line for line in text.splitlines() if "Total internalization / gross entry" in line)
    assert f"{detailed['nonspecific_internalization_fraction_gross']:.6g}" in total_line
    assert "Specific internalization / gross entry | 0 |" in text
    assert "cumulative lung retention over all passes" in text
