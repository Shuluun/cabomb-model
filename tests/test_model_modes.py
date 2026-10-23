import numpy as np
import pytest

from cabomb_model.model_modes import ModelModes, V1_1_MODES, V1_2_MODES
from cabomb_model.parameters import default_model_parameters
from cabomb_model.solver import simulate


def test_mode_constants_are_explicit():
    assert V1_1_MODES == ModelModes("linear_v1", "accessibility_v1_1")
    assert V1_2_MODES == ModelModes("saturable_v1_2", "sequential_v1_2")


@pytest.mark.parametrize(
    ("mps_model", "targeting_model"),
    [("unknown", "accessibility_v1_1"), ("linear_v1", "unknown")],
)
def test_invalid_model_modes_are_rejected(mps_model, targeting_model):
    with pytest.raises(ValueError, match="model"):
        ModelModes(mps_model, targeting_model)


def test_omitted_modes_are_exact_v11_modes():
    params = default_model_parameters()
    implicit = simulate(params, design="ocn-ysa")
    explicit = simulate(params, design="ocn-ysa", modes=V1_1_MODES)
    assert implicit.modes == V1_1_MODES
    assert explicit.modes == V1_1_MODES
    assert np.array_equal(implicit.t, explicit.t)
    assert np.allclose(implicit.y, explicit.y, rtol=0.0, atol=0.0)


@pytest.mark.parametrize(
    "symbol",
    [
        "Vmax_liver",
        "Km_liver",
        "k_liver_release",
        "Vmax_spleen",
        "Km_spleen",
        "k_spleen_release",
        "k_cross",
    ],
)
def test_v11_ignores_inactive_v12_parameter_domains(symbol):
    params = default_model_parameters()
    params[symbol] = -1.0

    result = simulate(
        params,
        modes=V1_1_MODES,
        t_end_h=0.1,
        n_timepoints=3,
    )

    assert result.success
    assert result.modes == V1_1_MODES
