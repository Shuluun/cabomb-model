"""EPR, MPS, IFP, and pH scenarios for exploratory model runs."""

from __future__ import annotations

from dataclasses import dataclass
from typing import Mapping


@dataclass(frozen=True)
class Scenario:
    """A scenario modifies uncertain entry, clearance, and pH context."""

    name: str
    epr_label: str = "medium"
    mps_label: str = "medium"
    epr_factor: float = 1.0
    mps_factor: float = 1.0
    ifp_mmHg: float = 20.0
    pcap_tumor_cm_s: float = 3.0e-7
    intracellular_pH: float = 5.5


def make_scenario(
    params: Mapping[str, float],
    *,
    epr: str = "medium",
    mps: str = "medium",
    ifp_mmHg: float | None = None,
    pcap_tumor_cm_s: float | None = None,
    intracellular_pH: float | None = None,
) -> Scenario:
    """Build a scenario from registry-backed factors."""

    epr_key = f"epr_factor_{epr}"
    mps_key = f"mps_factor_{mps}"
    if epr_key not in params:
        raise ValueError(f"Unknown EPR scenario {epr!r}.")
    if mps_key not in params:
        raise ValueError(f"Unknown MPS scenario {mps!r}.")
    return Scenario(
        name=f"epr-{epr}_mps-{mps}",
        epr_label=epr,
        mps_label=mps,
        epr_factor=params[epr_key],
        mps_factor=params[mps_key],
        ifp_mmHg=params["IFP_tumor_project"] if ifp_mmHg is None else float(ifp_mmHg),
        pcap_tumor_cm_s=(
            params["Pcap_tumor"] if pcap_tumor_cm_s is None else float(pcap_tumor_cm_s)
        ),
        intracellular_pH=(
            params["pH_internal"] if intracellular_pH is None else float(intracellular_pH)
        ),
    )

