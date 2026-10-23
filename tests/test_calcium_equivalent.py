"""Unit, domain and CLI checks for the independent stoichiometric calculator."""

import csv
import json

import pytest

from cabomb_model.calcium_equivalent import (
    calcium_equivalent_um,
    main,
    required_caco3_mass_pg,
)


def test_independent_si_hand_calculation_and_inverse():
    # 0.1 pg = 1e-13 g; half dissolved; 2 pL = 2e-12 L.
    expected_molar = (5e-14 / 100.087) / 2e-12
    concentration = calcium_equivalent_um(0.1, 0.5, 2.0)
    assert concentration == pytest.approx(249.78268906051736, rel=1e-12)
    assert concentration / 1e6 == pytest.approx(expected_molar, rel=1e-12)
    assert required_caco3_mass_pg(concentration, 0.5, 2) == pytest.approx(0.1)


@pytest.mark.parametrize("mass,fraction,volume", [(0, 0.5, 2), (0.1, 0, 2), (0, 0, 2)])
def test_zero_forward(mass, fraction, volume):
    assert calcium_equivalent_um(mass, fraction, volume) == 0


def test_scaling_and_legal_upper_fraction():
    base = calcium_equivalent_um(0.1, 0.5, 2)
    assert calcium_equivalent_um(0.2, 0.5, 2) == pytest.approx(2 * base)
    assert calcium_equivalent_um(0.1, 1, 2) == pytest.approx(2 * base)
    assert calcium_equivalent_um(0.1, 0.5, 4) == pytest.approx(base / 2)
    assert required_caco3_mass_pg(0, 1, 2) == 0
    assert required_caco3_mass_pg(base, 0.5, 4) == pytest.approx(0.2)


@pytest.mark.parametrize("args", [(-1, 0.5, 2), (1, -0.1, 2), (1, 1.1, 2),
                                  (1, 0.5, 0), (1, 0.5, -1),
                                  (float("nan"), 0.5, 2), (1, float("inf"), 2),
                                  (1, 0.5, float("nan")), (float("inf"), 0.5, 2),
                                  (1, 0.5, float("inf"))])
@pytest.mark.parametrize("function", [calcium_equivalent_um, required_caco3_mass_pg])
def test_illegal_inputs(function, args):
    with pytest.raises(ValueError):
        function(*args)


@pytest.mark.parametrize("target,message", [(0, "non-unique"), (1, "unattainable")])
def test_inverse_zero_dissolution(target, message):
    with pytest.raises(ValueError, match=message):
        required_caco3_mass_pg(target, 0, 2)


def test_cli_output_roundtrip_and_traceability(tmp_path):
    forward = tmp_path / "forward"
    assert main(["--mass-pg", "0.1", "--dissolved-fraction", "0.5", "--volume-pl", "2", "--output-dir", str(forward)]) == 0
    payload = json.loads((forward / "run_inputs.json").read_text(encoding="utf-8"))
    assert payload["calcium_equivalent_um"] == pytest.approx(249.78268906051736)
    assert len(payload["source_sha256"]) == 64
    assert "not free" in payload["interpretation"]
    inverse = tmp_path / "inverse"
    main(["--target-um", str(payload["calcium_equivalent_um"]), "--dissolved-fraction", "0.5", "--volume-pl", "2", "--output-dir", str(inverse)])
    rows = list(csv.DictReader((inverse / "input_registry.csv").read_text(encoding="utf-8").splitlines()))
    mass = next(row for row in rows if row["parameter_symbol"] == "m")
    assert float(mass["baseline_value"]) == pytest.approx(0.1)
    assert mass["value_type"] == "derived_conditional_mass"


def test_cli_invalid_input_creates_no_results(tmp_path):
    out = tmp_path / "invalid"
    with pytest.raises(SystemExit) as exc:
        main(["--mass-pg", "0.1", "--dissolved-fraction", "0.5", "--volume-pl", "0", "--output-dir", str(out)])
    assert exc.value.code == 2
    assert not out.exists()
