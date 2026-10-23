"""Carrier design switches for the frozen CABOMB V1 model."""

from __future__ import annotations

from dataclasses import dataclass


@dataclass(frozen=True)
class CarrierDesign:
    """Binary mechanism switches for one carrier design."""

    name: str
    i_ocn: int
    i_ysa: int
    description: str

    def __post_init__(self) -> None:
        if self.i_ocn not in (0, 1) or self.i_ysa not in (0, 1):
            raise ValueError("Design switches must be 0 or 1.")


DESIGNS: dict[str, CarrierDesign] = {
    "unmodified": CarrierDesign(
        name="unmodified",
        i_ocn=0,
        i_ysa=0,
        description="No OCN matrix retention and no YSA-EphA2 binding.",
    ),
    "ocn": CarrierDesign(
        name="ocn",
        i_ocn=1,
        i_ysa=0,
        description="OCN-mediated local matrix binding and retention only.",
    ),
    "ysa": CarrierDesign(
        name="ysa",
        i_ocn=0,
        i_ysa=1,
        description="YSA-EphA2 binding and receptor-mediated internalization only.",
    ),
    "ocn-ysa": CarrierDesign(
        name="ocn-ysa",
        i_ocn=1,
        i_ysa=1,
        description="Combined OCN retention and YSA internalization modules.",
    ),
}


def get_design(name: str | CarrierDesign) -> CarrierDesign:
    """Return a design by canonical name."""

    if isinstance(name, CarrierDesign):
        return name
    try:
        return DESIGNS[name]
    except KeyError as exc:
        valid = ", ".join(sorted(DESIGNS))
        raise ValueError(f"Unknown carrier design {name!r}. Valid designs: {valid}") from exc

