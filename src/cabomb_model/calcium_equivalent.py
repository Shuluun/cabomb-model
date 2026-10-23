"""Conditional calcium supply equivalent, independent of the delivery ODE.

Inputs refer to pure CaCO3 mass per cell, dissolved fraction and cytosolic volume.
The result is not free cytosolic calcium or a validated transfer/dose prediction.
"""

from __future__ import annotations

import argparse
import csv
import hashlib
import json
import math
from pathlib import Path

from cabomb_model.io import write_csv

from cabomb_model.resource_paths import resource_path

REGISTRY_PATH = resource_path("calcium_equivalent_parameters.csv")


def molar_mass_caco3() -> float:
    """Read the source-traceable CaCO3 molar mass (g/mol)."""
    with REGISTRY_PATH.open(encoding="utf-8", newline="") as handle:
        rows = {row["parameter_symbol"]: row for row in csv.DictReader(handle)}
    value = float(rows["M_CaCO3"]["baseline_value"])
    if not math.isfinite(value) or value <= 0:
        raise ValueError("M_CaCO3 must be finite and positive.")
    return value


def _validate(value: float, name: str, *, positive: bool = False) -> None:
    if not math.isfinite(value) or value < 0 or (positive and value == 0):
        raise ValueError(f"{name} must be finite and {'positive' if positive else 'nonnegative'}.")


def _validate_fraction_volume(dissolved_fraction: float, volume_pl: float) -> None:
    _validate(dissolved_fraction, "dissolved_fraction")
    if dissolved_fraction > 1:
        raise ValueError("dissolved_fraction must be in [0, 1].")
    _validate(volume_pl, "volume_pl", positive=True)


def calcium_equivalent_um(mass_pg: float, dissolved_fraction: float, volume_pl: float) -> float:
    """Return micromolar equivalent; pg/pL cancels the two 1e-12 SI factors."""
    _validate(mass_pg, "mass_pg")
    _validate_fraction_volume(dissolved_fraction, volume_pl)
    if mass_pg == 0 or dissolved_fraction == 0:
        return 0.0
    value = (mass_pg / volume_pl) * dissolved_fraction / molar_mass_caco3() * 1e6
    if not math.isfinite(value):
        raise ValueError("Calcium equivalent exceeds the finite numeric range.")
    return value


def required_caco3_mass_pg(target_um: float, dissolved_fraction: float, volume_pl: float) -> float:
    """Return required pure CaCO3 pg/cell, defined only for nonzero dissolution."""
    _validate(target_um, "target_um")
    _validate_fraction_volume(dissolved_fraction, volume_pl)
    if dissolved_fraction == 0:
        reason = "unattainable" if target_um > 0 else "non-unique"
        raise ValueError(f"Inverse mass is {reason} when dissolved_fraction is zero.")
    value = (target_um / 1e6) * molar_mass_caco3() * volume_pl / dissolved_fraction
    if not math.isfinite(value):
        raise ValueError("Required mass exceeds the finite numeric range.")
    return value


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    mode = parser.add_mutually_exclusive_group(required=True)
    mode.add_argument("--mass-pg", type=float, help="Pure CaCO3 uptake mass per cell, pg")
    mode.add_argument("--target-um", type=float, help="Conditional equivalent target, micromolar; not a treatment target")
    parser.add_argument("--dissolved-fraction", required=True, type=float)
    parser.add_argument("--volume-pl", required=True, type=float, help="Cytosolic volume per cell, pL")
    parser.add_argument("--input-note", default="User-specified conditional scenario; not measured experimental inputs")
    parser.add_argument("--output-dir", required=True, type=Path)
    args = parser.parse_args(argv)
    try:
        if args.mass_pg is not None:
            mass = args.mass_pg
            concentration = calcium_equivalent_um(mass, args.dissolved_fraction, args.volume_pl)
            direction = "forward"
        else:
            concentration = args.target_um
            mass = required_caco3_mass_pg(concentration, args.dissolved_fraction, args.volume_pl)
            direction = "inverse"
    except ValueError as exc:
        parser.error(str(exc))
    # Validate before creating any result files.
    row = {
        "direction": direction, "pure_caco3_mass_pg_per_cell": mass,
        "dissolved_fraction": args.dissolved_fraction, "cytosolic_volume_pl": args.volume_pl,
        "molar_mass_g_per_mol": molar_mass_caco3(), "calcium_equivalent_um": concentration,
        "input_note": args.input_note,
        "interpretation": "conditional calcium supply equivalent; not free cytosolic Ca2+",
    }
    output = args.output_dir
    write_csv(output / "calcium_equivalent.csv", [row])
    registry = list(csv.DictReader(REGISTRY_PATH.read_text(encoding="utf-8").splitlines()))
    values = {"m": mass, "f": args.dissolved_fraction, "V": args.volume_pl}
    for record in registry:
        if record["parameter_symbol"] in values:
            record["baseline_value"] = str(values[record["parameter_symbol"]])
            record["source_value_text"] = record["baseline_value"] + " " + record["unit"]
            record["status"] = "conditional_run_input"
            record["value_type"] = "user_supplied_condition" if direction == "forward" or record["parameter_symbol"] != "m" else "derived_conditional_mass"
            record["source_file"] = "command line; run_inputs.json"
            record["notes"] = args.input_note + "; bounds are domain constraints, not uncertainty intervals"
    write_csv(output / "input_registry.csv", registry)
    manifest = {**row, "registry_sha256": hashlib.sha256(REGISTRY_PATH.read_bytes()).hexdigest(),
                "source_sha256": hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
                "formula": "C_eq[uM] = 1e6 * m[pg] * f / (M[g/mol] * V[pL])"}
    (output / "run_inputs.json").write_text(json.dumps(manifest, ensure_ascii=False, indent=2, allow_nan=False) + "\n", encoding="utf-8")
    (output / "summary_zh.md").write_text(
        "# 每细胞供钙当量：条件换算\n\n"
        f"纯 CaCO₃ {mass:.9g} pg/细胞，溶解比例 {args.dissolved_fraction:.9g}，"
        f"胞质体积 {args.volume_pl:.9g} pL，折算供钙当量 **{concentration:.9g} µM**。\n\n"
        f"输入说明：{args.input_note}。\n\n"
        "1 mol CaCO₃ 提供 1 mol Ca；摩尔质量100.087 g/mol来自NIST Chemistry WebBook，"
        "具体来源及输入计量基础见 input_registry.csv。\n\n"
        "此处把全部已溶解载体钙按胞质体积折算；不证明跨区室转移或真实游离Ca²⁺浓度，"
        "不提供细胞死亡阈值或治疗剂量。输入是纯CaCO₃质量，不能直接用整颗粒质量、"
        "元素钙质量或归一化A_cargo替代。\n\n"
        "解析敏感性：其他量固定时，质量或溶解比例加倍，当量加倍（比例须≤1）；"
        "体积加倍，当量减半。这是公式比例关系，不是实验置信区间。\n",
        encoding="utf-8",
    )
    print(json.dumps(row, ensure_ascii=False))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
