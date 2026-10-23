import csv
import hashlib
from pathlib import Path

from cabomb_model.cli import main


def test_validate_registry_cli_passes():
    assert main(["validate-registry"]) == 0


def test_run_baseline_cli_passes(tmp_path: Path):
    assert main(["run-baseline", "--output-dir", str(tmp_path)]) == 0
    assert (tmp_path / "tables" / "baseline_summary.csv").exists()


def _tree_snapshot(path: Path) -> dict[str, tuple[int, int]]:
    if not path.exists():
        return {}
    return {
        str(child.relative_to(path)): (child.stat().st_size, child.stat().st_mtime_ns)
        for child in path.rglob("*")
        if child.is_file()
    }


def _sha256(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def test_run_analysis_cli_regenerates_parameter_registry_preview(tmp_path: Path):
    preview = tmp_path / "figures" / "parameters_master_preview.png"
    preview.parent.mkdir(parents=True)
    preview.write_bytes(b"stale preview")
    preview.unlink()
    inputs_before = _tree_snapshot(Path("inputs"))
    raw_before = _tree_snapshot(Path("data/raw"))

    assert main(["run-analysis", "--output-dir", str(tmp_path)]) == 0
    assert preview.exists()
    assert preview.stat().st_size > 0
    first_hash = _sha256(preview)

    assert main(["run-analysis", "--output-dir", str(tmp_path)]) == 0
    assert _sha256(preview) == first_hash
    assert _tree_snapshot(Path("inputs")) == inputs_before
    assert _tree_snapshot(Path("data/raw")) == raw_before


def test_v12_baseline_cli_writes_versioned_mode_columns(tmp_path: Path):
    assert main([
        "run-baseline",
        "--model-version", "v1.2",
        "--output-dir", str(tmp_path),
    ]) == 0
    with (tmp_path / "tables" / "baseline_summary.csv").open(
        newline="", encoding="utf-8"
    ) as handle:
        rows = list(csv.DictReader(handle))
    assert rows[0]["mps_model"] == "saturable_v1_2"
    assert rows[0]["targeting_model"] == "sequential_v1_2"
    with (tmp_path / "tables" / "baseline_timecourse.csv").open(
        newline="", encoding="utf-8"
    ) as handle:
        timecourse_rows = list(csv.DictReader(handle))
    assert timecourse_rows[0]["mps_model"] == "saturable_v1_2"
    assert timecourse_rows[0]["targeting_model"] == "sequential_v1_2"


def test_v12_analysis_cli_writes_versioned_mode_columns(tmp_path: Path):
    assert main([
        "run-analysis",
        "--model-version", "v1.2",
        "--output-dir", str(tmp_path),
    ]) == 0
    with (tmp_path / "tables" / "design_comparison.csv").open(
        newline="", encoding="utf-8"
    ) as handle:
        rows = list(csv.DictReader(handle))
    assert {row["mps_model"] for row in rows} == {"saturable_v1_2"}
    assert {row["targeting_model"] for row in rows} == {"sequential_v1_2"}


def test_v12_build_report_cli_uses_requested_model_modes(tmp_path: Path):
    report_path = tmp_path / "model_report_v1_2.md"
    html_path = tmp_path / "model_report_v1_2.html"

    assert main([
        "build-report",
        "--model-version", "v1.2",
        "--output-dir", str(tmp_path / "outputs"),
        "--report-path", str(report_path),
        "--html-path", str(html_path),
    ]) == 0

    report = report_path.read_text(encoding="utf-8")
    assert "CABOMB V1.2" in report
    assert "sequential structure assumption" in report
    assert "OCN accessibility sweep" not in report


def test_v13_baseline_cli_writes_explicit_circulation_columns(tmp_path: Path):
    assert main([
        "run-baseline",
        "--model-version", "v1.3",
        "--output-dir", str(tmp_path),
    ]) == 0
    with (tmp_path / "tables" / "baseline_summary.csv").open(
        newline="", encoding="utf-8"
    ) as handle:
        rows = list(csv.DictReader(handle))
    assert rows[0]["circulation_model"] == "explicit_v1_3"
    assert rows[0]["mps_model"] == "saturable_v1_2"
    assert rows[0]["targeting_model"] == "sequential_v1_2"
