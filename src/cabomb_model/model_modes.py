"""Explicit model-version switches for CABOMB simulations."""

from __future__ import annotations

from dataclasses import dataclass


MPS_MODELS = ("linear_v1", "saturable_v1_2")
TARGETING_MODELS = ("accessibility_v1_1", "sequential_v1_2", "parallel_dual_access_v1_3")
CIRCULATION_MODELS = ("central_v1", "explicit_v1_3")


@dataclass(frozen=True)
class ModelModes:
    """Select MPS, tumor-targeting, and circulation model variants."""

    mps_model: str = "linear_v1"
    targeting_model: str = "accessibility_v1_1"
    circulation_model: str = "central_v1"

    def __post_init__(self) -> None:
        if self.mps_model not in MPS_MODELS:
            raise ValueError(f"Unknown MPS model {self.mps_model!r}.")
        if self.targeting_model not in TARGETING_MODELS:
            raise ValueError(f"Unknown targeting model {self.targeting_model!r}.")
        if self.circulation_model not in CIRCULATION_MODELS:
            raise ValueError(f"Unknown circulation model {self.circulation_model!r}.")


V1_1_MODES = ModelModes()
V1_2_MODES = ModelModes("saturable_v1_2", "sequential_v1_2")
V1_3_MODES = ModelModes("saturable_v1_2", "sequential_v1_2", "explicit_v1_3")
