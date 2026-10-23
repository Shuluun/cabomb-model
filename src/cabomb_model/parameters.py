"""Parameter registry loading and executable parameter validation."""

from __future__ import annotations

import csv
import math
from dataclasses import dataclass
from pathlib import Path
from typing import Iterable, Mapping

from cabomb_model.model_modes import ModelModes, V1_1_MODES

REQUIRED_REGISTRY_COLUMNS: tuple[str, ...] = (
    "status",
    "parameter_symbol",
    "definition",
    "baseline_value",
    "lower_bound",
    "upper_bound",
    "unit",
    "source_type",
    "source_file",
    "page_table_section",
    "species",
    "nanoparticle_material",
    "nanoparticle_size",
    "administration_route",
    "value_type",
    "confidence_level",
    "notes",
    "source_value_text",
)

NUMERIC_MODEL_SYMBOLS: tuple[str, ...] = (
    "D0",
    "Vblood",
    "Vliv",
    "Vspl",
    "Vtum",
    "Pcap_tumor",
    "Pcap_tumor_ref",
    "P_EPR",
    "IFP_tumor_project",
    "IFP_ref",
    "ifp_entry_sensitivity",
    "pH_internal",
    "pH_blood_phys",
    "pH_endosome_lysosome",
    "k_entry",
    "k_return",
    "k_liver_uptake",
    "k_spleen_uptake",
    "k_liver_elim",
    "k_spleen_elim",
    "k_blood_clear",
    "k_ocn_on",
    "k_ocn_off",
    "k_ocn_loss",
    "alpha_ocn_cell_access",
    "k_ysa_on",
    "k_ysa_off",
    "R_EphA2_eff",
    "k_int",
    "k_int_ns",
    "k_deg_int_acid",
    "k_deg_phys",
    "k_tumor_deg",
    "Y_cargo",
    "k_cargo_loss",
    "epr_factor_low",
    "epr_factor_medium",
    "epr_factor_high",
    "mps_factor_low",
    "mps_factor_medium",
    "mps_factor_high",
    "Vmax_liver",
    "Km_liver",
    "k_liver_release",
    "Vmax_spleen",
    "Km_spleen",
    "k_spleen_release",
    "k_cross",
)
SATURABLE_MPS_V12_SYMBOLS: tuple[str, ...] = (
    "Vmax_liver",
    "Km_liver",
    "k_liver_release",
    "Vmax_spleen",
    "Km_spleen",
    "k_spleen_release",
)
SEQUENTIAL_TARGETING_V12_SYMBOLS: tuple[str, ...] = ("k_cross",)

# V1.3 explicit-circulation symbols. These are loaded by the V1.3 code path
# only; they are absent from the frozen V1.1/V1.2 NUMERIC_MODEL_SYMBOLS set.
V13_FLOW_SYMBOLS: tuple[str, ...] = ("Qco", "Qliv", "Qspl", "Qtum")
V13_VOLUME_SYMBOLS: tuple[str, ...] = (
    "V_central_venous",
    "V_right_heart",
    "V_lung_vascular",
    "V_left_heart",
    "V_systemic_arterial",
    "V_liver_vascular",
    "V_spleen_vascular",
    "V_tumor_vascular",
    "V_rest_vascular",
)
V13_LUNG_RET_SYMBOLS: tuple[str, ...] = (
    "f_lung_ret_100nm",
    "f_lung_ret_200nm",
    "f_lung_ret_500nm",
    "f_lung_ret_800nm",
    "f_lung_ret_2um",
    "f_lung_ret_3um",
)
V13_SIZE_MPS_SYMBOLS: tuple[str, ...] = (
    "size_mps_100nm",
    "size_mps_200nm",
    "size_mps_500nm",
    "size_mps_800nm",
    "size_mps_2um",
    "size_mps_3um",
)
V13_SIZE_EXTRA_SYMBOLS: tuple[str, ...] = (
    "size_extra_100nm",
    "size_extra_200nm",
    "size_extra_500nm",
    "size_extra_800nm",
    "size_extra_2um",
    "size_extra_3um",
)
V13_SIZE_INT_SYMBOLS: tuple[str, ...] = (
    "size_int_100nm",
    "size_int_200nm",
    "size_int_500nm",
    "size_int_800nm",
    "size_int_2um",
    "size_int_3um",
)
V13_TARGETING_SYMBOLS: tuple[str, ...] = ("alpha_ysa_direct_access",)
V13_CIRCULATION_SYMBOLS: tuple[str, ...] = (
    V13_VOLUME_SYMBOLS
    + V13_LUNG_RET_SYMBOLS
    + V13_SIZE_MPS_SYMBOLS
    + V13_SIZE_EXTRA_SYMBOLS
    + V13_SIZE_INT_SYMBOLS
)
V13_EXECUTABLE_SYMBOLS: tuple[str, ...] = (
    V13_FLOW_SYMBOLS + V13_CIRCULATION_SYMBOLS + V13_TARGETING_SYMBOLS
)


@dataclass(frozen=True)
class RegistryValidation:
    """Structured registry validation result."""

    ok: bool
    errors: tuple[str, ...]
    warnings: tuple[str, ...]


def repo_root() -> Path:
    """Find the repository root from this installed source tree."""

    from cabomb_model.resource_paths import source_root
    return source_root() or Path(__file__).resolve().parent


def default_registry_path() -> Path:
    """Return the default processed parameter registry path."""

    from cabomb_model.resource_paths import selected_registry_path
    return selected_registry_path()


def load_registry(path: str | Path | None = None) -> list[dict[str, str]]:
    """Load parameter registry rows as dictionaries."""

    registry_path = Path(path) if path else default_registry_path()
    with registry_path.open("r", encoding="utf-8-sig", newline="") as handle:
        reader = csv.DictReader(handle)
        return [dict(row) for row in reader]


def _parse_float(value: str, symbol: str, field: str) -> float:
    if value is None or str(value).strip() == "":
        raise ValueError(f"{symbol} has empty {field}.")
    try:
        parsed = float(str(value).strip())
    except ValueError as exc:
        raise ValueError(f"{symbol} has non-numeric {field}: {value!r}.") from exc
    return parsed


def validate_registry(
    records: Iterable[Mapping[str, str]],
    *,
    numeric_symbols: Iterable[str] = NUMERIC_MODEL_SYMBOLS,
) -> RegistryValidation:
    """Validate required metadata and numeric bounds for registry parameters."""

    rows = list(records)
    errors: list[str] = []
    warnings: list[str] = []
    if not rows:
        errors.append("Parameter registry has no rows.")
        return RegistryValidation(False, tuple(errors), tuple(warnings))

    missing_columns = [col for col in REQUIRED_REGISTRY_COLUMNS if col not in rows[0]]
    if missing_columns:
        errors.append(f"Missing columns: {', '.join(missing_columns)}")

    by_symbol: dict[str, Mapping[str, str]] = {}
    for idx, row in enumerate(rows, start=2):
        symbol = (row.get("parameter_symbol") or "").strip()
        if not symbol:
            errors.append(f"Row {idx} is missing parameter_symbol.")
            continue
        if symbol in by_symbol:
            errors.append(f"Duplicate parameter_symbol {symbol!r}.")
        by_symbol[symbol] = row
        for field in REQUIRED_REGISTRY_COLUMNS:
            if str(row.get(field, "")).strip() == "":
                errors.append(f"{symbol} is missing required metadata field {field}.")
        try:
            baseline = _parse_float(row.get("baseline_value", ""), symbol, "baseline_value")
            lower = _parse_float(row.get("lower_bound", ""), symbol, "lower_bound")
            upper = _parse_float(row.get("upper_bound", ""), symbol, "upper_bound")
        except ValueError as exc:
            errors.append(str(exc))
            continue
        if lower > baseline or baseline > upper:
            errors.append(
                f"{symbol} baseline {baseline} is outside bounds [{lower}, {upper}]."
            )
        if lower < 0 and symbol not in {"IFP_tumor_project"}:
            errors.append(f"{symbol} has negative lower bound {lower}.")
        value_type = row.get("value_type", "").lower()
        status = row.get("status", "").lower()
        if (
            "assumed" in value_type
            or "exploratory" in value_type
            or "exploratory" in status
        ):
            warnings.append(f"{symbol} is exploratory and must not be treated as measured.")

    for symbol in numeric_symbols:
        row = by_symbol.get(symbol)
        if row is None:
            errors.append(f"Executable parameter {symbol} is absent from registry.")

    return RegistryValidation(len(errors) == 0, tuple(errors), tuple(warnings))


def default_model_parameters(path: str | Path | None = None) -> dict[str, float]:
    """Return numeric baseline parameters for the executable ODE model."""

    records = load_registry(path)
    validation = validate_registry(records)
    if not validation.ok:
        joined = "\n".join(validation.errors)
        raise ValueError(f"Parameter registry is not executable:\n{joined}")
    by_symbol = {row["parameter_symbol"]: row for row in records}
    return {
        symbol: _parse_float(by_symbol[symbol]["baseline_value"], symbol, "baseline_value")
        for symbol in NUMERIC_MODEL_SYMBOLS
    }


def inactive_mode_parameter_symbols(modes: ModelModes) -> frozenset[str]:
    """Return V1.2-only parameters that are inactive for the selected modes."""

    inactive: set[str] = set()
    if modes.mps_model != "saturable_v1_2":
        inactive.update(SATURABLE_MPS_V12_SYMBOLS)
    if modes.targeting_model != "sequential_v1_2":
        inactive.update(SEQUENTIAL_TARGETING_V12_SYMBOLS)
    return frozenset(inactive)


def _required_mode_parameter(params: Mapping[str, float], symbol: str) -> float:
    try:
        return float(params[symbol])
    except KeyError as exc:
        raise ValueError(
            f"{symbol} is required by the active model modes."
        ) from exc


def validate_executable_parameters(
    params: Mapping[str, float],
    modes: ModelModes = V1_1_MODES,
) -> None:
    """Validate model-wide domains not expressible as non-negativity alone."""

    alpha = float(params.get("alpha_ocn_cell_access", 0.0))
    if not math.isfinite(alpha) or not 0.0 <= alpha <= 1.0:
        raise ValueError("alpha_ocn_cell_access must be finite and within [0, 1].")

    if modes.mps_model == "saturable_v1_2":
        for symbol in ("Km_liver", "Km_spleen", "Vblood"):
            value = _required_mode_parameter(params, symbol)
            if not math.isfinite(value) or value <= 0.0:
                raise ValueError(f"{symbol} must be finite and positive.")
        for symbol in SATURABLE_MPS_V12_SYMBOLS:
            if symbol.startswith("Km_"):
                continue
            value = _required_mode_parameter(params, symbol)
            if not math.isfinite(value) or value < 0.0:
                raise ValueError(f"{symbol} must be finite and nonnegative.")

    if modes.targeting_model == "sequential_v1_2":
        value = _required_mode_parameter(params, "k_cross")
        if not math.isfinite(value) or value < 0.0:
            raise ValueError("k_cross must be finite and nonnegative.")


def parameter_bounds(path: str | Path | None = None) -> dict[str, tuple[float, float]]:
    """Return lower and upper bounds for executable model parameters."""

    records = load_registry(path)
    validation = validate_registry(records)
    if not validation.ok:
        joined = "\n".join(validation.errors)
        raise ValueError(f"Parameter registry is not executable:\n{joined}")
    by_symbol = {row["parameter_symbol"]: row for row in records}
    return {
        symbol: (
            _parse_float(by_symbol[symbol]["lower_bound"], symbol, "lower_bound"),
            _parse_float(by_symbol[symbol]["upper_bound"], symbol, "upper_bound"),
        )
        for symbol in NUMERIC_MODEL_SYMBOLS
    }
