from pathlib import Path
from xml.etree import ElementTree
from zipfile import ZipFile

import pytest

from cabomb_model.model_modes import V1_2_MODES
from cabomb_model.parameters import (
    default_model_parameters,
    load_registry,
    validate_registry,
)
from cabomb_model.solver import simulate


XLSX_PATH = (
    Path(__file__).resolve().parents[1]
    / "data"
    / "processed"
    / "parameters_master.xlsx"
)
SPREADSHEET_NAMESPACE = "http://schemas.openxmlformats.org/spreadsheetml/2006/main"
V12_SYMBOLS = {
    "Vmax_liver",
    "Km_liver",
    "k_liver_release",
    "Vmax_spleen",
    "Km_spleen",
    "k_spleen_release",
    "k_cross",
}
V12_APPROVED_SCAFFOLD = {
    "Vmax_liver": (
        2.08,
        0.52,
        6.24,
        "normalized carrier / h",
    ),
    "Km_liver": (
        1.0,
        0.2,
        5.0,
        "normalized carrier / apparent L",
    ),
    "k_liver_release": (0.01, 0.0, 0.10, "h^-1"),
    "Vmax_spleen": (
        0.78,
        0.104,
        2.60,
        "normalized carrier / h",
    ),
    "Km_spleen": (
        1.0,
        0.2,
        5.0,
        "normalized carrier / apparent L",
    ),
    "k_spleen_release": (0.01, 0.0, 0.10, "h^-1"),
    "k_cross": (0.10, 0.0, 0.60, "h^-1"),
}
V12_SOURCE_VALUE_TEXT = (
    "The number is an approved, transparent V1.2 implementation assumption; "
    "no measured CABOMB-specific value was supplied."
)


def test_parameter_registry_has_required_metadata_and_executable_values():
    validation = validate_registry(load_registry())
    assert validation.ok, "\n".join(validation.errors)


def test_alpha_ocn_cell_access_is_explicitly_exploratory():
    rows = {row["parameter_symbol"]: row for row in load_registry()}
    alpha = rows["alpha_ocn_cell_access"]
    assert float(alpha["baseline_value"]) == 0.0
    assert float(alpha["lower_bound"]) == 0.0
    assert float(alpha["upper_bound"]) == 1.0
    assert alpha["unit"] == "dimensionless"
    assert "exploratory" in alpha["value_type"].lower()
    assert "not" in alpha["notes"].lower()


def test_alpha_ocn_cell_access_cites_the_approved_v11_specification():
    rows = {row["parameter_symbol"]: row for row in load_registry()}
    alpha = rows["alpha_ocn_cell_access"]
    assert alpha["source_type"] == "model implementation assumption"
    assert alpha["source_file"] == (
        "docs/superpowers/specs/"
        "2026-07-17-cabomb-v1.1-metrics-ocn-access-design.md"
    )
    assert alpha["page_table_section"] == "Section 4. New exploratory parameter"
    assert "docs/open_questions.md" not in alpha["source_file"]
    assert "phenomenological" in alpha["source_value_text"].lower()


def test_xlsx_registry_numeric_columns_are_stored_as_numbers():
    with ZipFile(XLSX_PATH) as workbook:
        sheet = ElementTree.fromstring(workbook.read("xl/worksheets/sheet1.xml"))

    cells = {
        cell.attrib["r"]: cell
        for cell in sheet.iter(f"{{{SPREADSHEET_NAMESPACE}}}c")
    }
    alpha_row = next(
        reference[1:]
        for reference, cell in cells.items()
        if reference.startswith("B")
        and cell.findtext(f"{{{SPREADSHEET_NAMESPACE}}}v")
        == "alpha_ocn_cell_access"
    )

    numeric_cells = [
        cell
        for reference, cell in cells.items()
        if reference[0] in "DEF" and reference[1:] != "1"
    ]
    assert numeric_cells
    assert all(
        cell.attrib.get("t") not in {"inlineStr", "s", "str"}
        for cell in numeric_cells
    )

    for column, expected in zip("DEF", (0.0, 0.0, 1.0), strict=True):
        cell = cells[f"{column}{alpha_row}"]
        assert float(cell.findtext(f"{{{SPREADSHEET_NAMESPACE}}}v")) == expected


def test_xlsx_registry_matches_csv_for_every_field():
    csv_rows = load_registry()
    with ZipFile(XLSX_PATH) as workbook:
        sheet = ElementTree.fromstring(workbook.read("xl/worksheets/sheet1.xml"))

    cells = {
        cell.attrib["r"]: cell.findtext(f"{{{SPREADSHEET_NAMESPACE}}}v") or ""
        for cell in sheet.iter(f"{{{SPREADSHEET_NAMESPACE}}}c")
    }
    columns = "ABCDEFGHIJKLMNOPQR"
    headers = [cells[f"{column}1"] for column in columns]
    assert headers == list(csv_rows[0])

    assert len(csv_rows) == 108
    for row_number, csv_row in enumerate(csv_rows, start=2):
        for column, field in zip(columns, headers, strict=True):
            workbook_value = cells.get(f"{column}{row_number}", "")
            csv_value = csv_row[field]
            if column in "DEF" and csv_value:
                assert float(workbook_value) == float(csv_value)
            else:
                assert workbook_value == csv_value

    alpha = next(row for row in csv_rows if row["parameter_symbol"] == "alpha_ocn_cell_access")
    alpha_row = next(
        row_number
        for row_number in range(2, len(csv_rows) + 2)
        if cells[f"B{row_number}"] == "alpha_ocn_cell_access"
    )
    assert cells[f"I{alpha_row}"] == alpha["source_file"]
    assert cells[f"J{alpha_row}"] == alpha["page_table_section"]


def test_v12_registry_rows_are_explicit_exploratory_assumptions():
    rows = {row["parameter_symbol"]: row for row in load_registry()}
    assert V12_SYMBOLS <= set(rows)
    for symbol in V12_SYMBOLS:
        row = rows[symbol]
        assert "exploratory" in row["value_type"].lower()
        assert row["confidence_level"].lower() == "low"
        assert "2026-07-22-cabomb-v1.2" in row["source_file"]
        assert "not" in row["notes"].lower()


def test_v12_registry_rows_match_the_approved_numeric_scaffold_exactly():
    rows = {row["parameter_symbol"]: row for row in load_registry()}
    for symbol, (baseline, lower, upper, unit) in V12_APPROVED_SCAFFOLD.items():
        row = rows[symbol]
        assert float(row["baseline_value"]) == baseline
        assert float(row["lower_bound"]) == lower
        assert float(row["upper_bound"]) == upper
        assert row["unit"] == unit
        assert row["source_value_text"] == V12_SOURCE_VALUE_TEXT


@pytest.mark.parametrize("symbol", ["Km_liver", "Km_spleen"])
def test_v12_rejects_nonpositive_km(symbol):
    params = default_model_parameters()
    params[symbol] = 0.0
    with pytest.raises(
        ValueError,
        match=rf"^{symbol} must be finite and positive\.$",
    ):
        simulate(params, modes=V1_2_MODES)


@pytest.mark.parametrize(
    "symbol",
    ["Vmax_liver", "k_liver_release", "Vmax_spleen", "k_spleen_release", "k_cross"],
)
def test_v12_rejects_negative_active_rates(symbol):
    params = default_model_parameters()
    params[symbol] = -1.0
    with pytest.raises(
        ValueError,
        match=rf"^{symbol} must be finite and nonnegative\.$",
    ):
        simulate(params, modes=V1_2_MODES)


@pytest.mark.parametrize("symbol", ["Vmax_liver", "k_cross"])
def test_v12_reports_missing_mode_active_parameters(symbol):
    params = default_model_parameters()
    del params[symbol]

    with pytest.raises(
        ValueError,
        match=rf"^{symbol} is required by the active model modes\.$",
    ):
        simulate(params, modes=V1_2_MODES)


def test_alpha_ysa_direct_access_is_explicitly_exploratory():
    rows = {row["parameter_symbol"]: row for row in load_registry()}
    alpha = rows["alpha_ysa_direct_access"]
    assert float(alpha["baseline_value"]) == 0.5
    assert float(alpha["lower_bound"]) == 0.0
    assert float(alpha["upper_bound"]) == 1.0
    assert alpha["unit"] == "dimensionless"
    assert "exploratory" in alpha["value_type"].lower()
    assert alpha["confidence_level"].lower() == "low"
    assert "not" in alpha["notes"].lower()


def test_v13_six_size_factor_rows_are_present_and_obsolete_symbols_removed():
    rows = {row["parameter_symbol"]: row for row in load_registry()}
    for size in ("100nm", "200nm", "500nm", "800nm", "2um", "3um"):
        for prefix in ("f_lung_ret_", "size_mps_", "size_extra_", "size_int_"):
            assert prefix + size in rows, f"missing {prefix}{size}"
    for obsolete in ("f_lung_ret_2_3um", "size_mps_2_3um", "size_extra_2_3um", "size_int_2_3um"):
        assert obsolete not in rows

