from pathlib import Path

import pytest

import cabomb_model.reporting as reporting_module
from cabomb_model.reporting import (
    _markdown_to_html,
    build_html_report,
    build_markdown_report,
)
from cabomb_model.model_modes import ModelModes, V1_2_MODES


def _stub_analysis_results() -> dict[str, object]:
    return {
        "baseline": None,
        "design_comparison": [],
        "scenario_comparison": [],
        "ocn_accessibility_sweep": [],
        "local_sensitivity": [],
        "monte_carlo_design_advantage_summary": None,
    }


def _write_stubbed_report(
    tmp_path: Path,
    monkeypatch: pytest.MonkeyPatch,
    modes: ModelModes,
) -> str:
    monkeypatch.setattr(
        reporting_module,
        "run_all_analyses",
        lambda output_dir, *, modes: _stub_analysis_results(),
    )
    return reporting_module.build_markdown_report(
        tmp_path / "outputs",
        tmp_path / f"{modes.mps_model}-{modes.targeting_model}.md",
        modes=modes,
    ).read_text(encoding="utf-8")


def test_report_generation_is_deterministic(tmp_path: Path):
    report = tmp_path / "model_report.md"
    first = build_markdown_report(tmp_path / "outputs", report).read_text(encoding="utf-8")
    second = build_markdown_report(tmp_path / "outputs", report).read_text(encoding="utf-8")
    assert first == second
    assert "Generated:" not in second
    assert "wall-clock timestamp" in second
    assert "parameters_master_preview.png" in second


def test_html_report_renders_tables(tmp_path: Path):
    markdown = build_markdown_report(tmp_path / "outputs", tmp_path / "model_report.md")
    html = build_html_report(markdown, tmp_path / "model_report.html").read_text(
        encoding="utf-8"
    )
    assert "<table>" in html
    assert "<pre>" not in html
    assert "| Rank | Parameter" not in html


def test_html_report_renders_sensitivity_targets_as_level_three_headings():
    html = _markdown_to_html(
        "## Multi-target local sensitivity\n\n"
        "### `tumor_auc`\n\n"
        "### `cumulative_internalization`\n\n"
        "### `cumulative_cargo_release`\n"
    )

    assert "<h3><code>tumor_auc</code></h3>" in html
    assert "<h3><code>cumulative_internalization</code></h3>" in html
    assert "<h3><code>cumulative_cargo_release</code></h3>" in html
    assert "<p>### " not in html


def test_v11_report_distinguishes_metric_classes_and_probability_boundary(
    tmp_path: Path,
):
    report_path = build_markdown_report(
        tmp_path / "outputs", tmp_path / "model_report.md"
    )
    report = report_path.read_text(encoding="utf-8")

    assert "CABOMB V1.1" in report
    assert "## State stock, exposure AUC, and cumulative flux" in report
    assert "State stock" in report
    assert "Exposure AUC" in report
    assert "Cumulative flux" in report
    assert "alpha_ocn_cell_access" in report
    assert "A_dual_bound" in report
    assert "not probabilities of clinical success" in report
    assert "Peak tumor-associated fraction" in report
    assert "Peak tumor stock" not in report

    html = build_html_report(
        report_path, tmp_path / "model_report.html"
    ).read_text(encoding="utf-8")
    assert "CABOMB V1.1" in html
    assert "State stock" in html
    assert "not probabilities of clinical success" in html
    assert "Peak tumor-associated fraction" in html
    assert "Peak tumor stock" not in html


def test_streamlit_declares_exploratory_alpha_and_metric_classes():
    source = (Path(__file__).parents[1] / "app" / "streamlit_app.py").read_text(
        encoding="utf-8"
    )

    assert 'params["alpha_ocn_cell_access"] = st.slider(' in source
    assert '"Peak tumor-associated fraction"' in source
    assert '"Peak tumor state stock (normalized carrier)"' not in source
    assert "min_value=0.0" in source
    assert "max_value=1.0" in source
    assert "value=0.0" in source
    assert "step=0.05" in source
    assert "Phenomenological and unvalidated" in source
    assert "State stock is amount present at one time" in source
    assert "not probabilities of clinical success" in source


def test_streamlit_exposes_modes_and_only_active_mechanism_controls():
    source = (Path(__file__).parents[1] / "app" / "streamlit_app.py").read_text(
        encoding="utf-8"
    )

    assert "V1.1 accessibility" in source
    assert "V1.2 sequential" in source
    assert "Vmax_liver" in source
    assert "Km_liver" in source
    assert "k_liver_release" in source
    assert "k_cross" in source
    assert "not a measured probability" in source


def test_streamlit_displays_explicit_model_mode_identifiers():
    source = (Path(__file__).parents[1] / "app" / "streamlit_app.py").read_text(
        encoding="utf-8"
    )

    assert 'f"Active MPS model: {modes.mps_model}"' in source
    assert 'f"Active targeting model: {modes.targeting_model}"' in source
    assert '"mps_model",' in source
    assert '"targeting_model",' in source


@pytest.mark.parametrize(
    ("modes", "mps_heading", "targeting_heading"),
    [
        (
            ModelModes("linear_v1", "accessibility_v1_1"),
            "## Linear MPS mechanism",
            "## Accessibility targeting mechanism",
        ),
        (
            ModelModes("linear_v1", "sequential_v1_2"),
            "## Linear MPS mechanism",
            "## Sequential targeting mechanism",
        ),
        (
            ModelModes("saturable_v1_2", "accessibility_v1_1"),
            "## Saturable MPS mechanism",
            "## Accessibility targeting mechanism",
        ),
        (
            ModelModes("saturable_v1_2", "sequential_v1_2"),
            "## Saturable MPS mechanism",
            "## Sequential targeting mechanism",
        ),
    ],
)
def test_report_describes_mps_and_targeting_modes_orthogonally(
    tmp_path: Path,
    monkeypatch: pytest.MonkeyPatch,
    modes: ModelModes,
    mps_heading: str,
    targeting_heading: str,
):
    report = _write_stubbed_report(tmp_path, monkeypatch, modes)

    assert f"`mps_model` | `{modes.mps_model}`" in report
    assert f"`targeting_model` | `{modes.targeting_model}`" in report
    assert mps_heading in report
    assert targeting_heading in report
    assert ("Michaelis-Menten" in report) == (
        modes.mps_model == "saturable_v1_2"
    )
    assert ("## Accessibility targeting mechanism" in report) == (
        modes.targeting_model == "accessibility_v1_1"
    )
    assert ("## Sequential targeting mechanism" in report) == (
        modes.targeting_model == "sequential_v1_2"
    )


@pytest.mark.parametrize(
    ("targeting_model", "required_priority", "forbidden_term"),
    [
        (
            "accessibility_v1_1",
            "`alpha_ocn_cell_access`",
            "strict OCN-to-YSA sequence assumption",
        ),
        (
            "sequential_v1_2",
            "`k_cross`",
            "alpha",
        ),
    ],
)
def test_wet_lab_priority_follows_targeting_mode(
    tmp_path: Path,
    monkeypatch: pytest.MonkeyPatch,
    targeting_model: str,
    required_priority: str,
    forbidden_term: str,
):
    report = _write_stubbed_report(
        tmp_path,
        monkeypatch,
        ModelModes("linear_v1", targeting_model),
    )
    wet_lab_section = report.split("## Wet-lab implications", 1)[1].split(
        "## Interpretation boundary", 1
    )[0]

    assert required_priority in wet_lab_section
    assert forbidden_term not in wet_lab_section


@pytest.mark.parametrize(
    "modes",
    [
        ModelModes("linear_v1", "sequential_v1_2"),
        ModelModes("saturable_v1_2", "sequential_v1_2"),
    ],
)
def test_nonempty_sensitivity_boundary_names_both_selected_modes(
    tmp_path: Path,
    monkeypatch: pytest.MonkeyPatch,
    modes: ModelModes,
):
    results = _stub_analysis_results()
    results["local_sensitivity"] = [
        {
            "target_metric": "tumor_auc",
            "parameter_symbol": "k_cross",
            "sensitivity_method": "central_relative_20pct",
            "abs_normalized_sensitivity": 0.25,
        }
    ]
    monkeypatch.setattr(
        reporting_module,
        "run_all_analyses",
        lambda output_dir, *, modes: results,
    )

    report = reporting_module.build_markdown_report(
        tmp_path / "outputs",
        tmp_path / f"sensitivity-{modes.mps_model}-{modes.targeting_model}.md",
        modes=modes,
    ).read_text(encoding="utf-8")

    assert (
        "Local sensitivities cover only active parameters for the selected "
        f"`mps_model={modes.mps_model}` and "
        f"`targeting_model={modes.targeting_model}` mechanisms."
    ) in report
    if modes.mps_model == "linear_v1":
        assert "active saturable-MPS and sequential-targeting parameters" not in report


def test_v12_equations_document_rhs_and_four_design_paths():
    equations = (Path(__file__).parents[1] / "docs" / "equations.md").read_text(
        encoding="utf-8"
    )
    normalized_equations = " ".join(equations.split())

    for rhs_term in (
        "- J_liver_MPS - J_spleen_MPS",
        "+ J_liver_release + J_spleen_release",
        "J_liver_MPS - J_liver_release - k_liver_elim A_liver",
        "J_spleen_MPS - J_spleen_release - k_spleen_elim A_spleen",
        "+ k_liver_elim A_liver + k_spleen_elim A_spleen",
    ):
        assert rhs_term in equations
    assert (
        "J_liver_release and J_spleen_release do not enter "
        "`A_elim_deg` or `A_cargo_lost`"
    ) in normalized_equations
    assert (
        "J_ysa_from_free = I_YSA * (1 - I_OCN) * k_ysa_on "
        "* R_EphA2_eff * A_tumor_free"
    ) in equations
    assert "V1 non-targeting terms" not in equations
    assert "For the OCN-YSA design, `J_ysa_from_free = 0`" in equations
    for design in ("unmodified", "OCN only", "YSA only", "OCN-YSA"):
        assert f"| {design} |" in equations


def test_report_states_unmodeled_organ_passage_boundary(
    tmp_path: Path,
    monkeypatch: pytest.MonkeyPatch,
):
    report = _write_stubbed_report(tmp_path, monkeypatch, V1_2_MODES)

    assert (
        "Unmodeled organ passage is not automatically terminal elimination"
        in report
    )


def test_report_generated_output_paths_follow_requested_output_directory(
    tmp_path: Path,
    monkeypatch: pytest.MonkeyPatch,
):
    output_dir = tmp_path / "outputs" / "v1_2"
    monkeypatch.setattr(
        reporting_module,
        "run_all_analyses",
        lambda output_dir, *, modes: _stub_analysis_results(),
    )
    report = build_markdown_report(
        output_dir,
        tmp_path / "model_report_v1_2.md",
        modes=V1_2_MODES,
    ).read_text(encoding="utf-8")

    expected_root = output_dir.as_posix()
    assert f"`{expected_root}/tables/baseline_timecourse.csv`" in report
    assert f"`{expected_root}/figures/baseline_timecourse.svg`" in report
    assert "`outputs/tables/baseline_timecourse.csv`" not in report


def test_v12_report_states_all_scientific_boundaries(tmp_path: Path):
    report = build_markdown_report(
        tmp_path / "outputs",
        tmp_path / "model_report_v1_2.md",
        modes=V1_2_MODES,
    ).read_text(encoding="utf-8")

    for statement in (
        "not a full PBPK model",
        "not an absolute mass concentration",
        "one candidate mechanism",
        "sequential structure assumption",
        "not a dynamic receptor model",
        "not probabilities of clinical success",
        "Unmodeled organ passage is not automatically terminal elimination",
    ):
        assert statement in report
