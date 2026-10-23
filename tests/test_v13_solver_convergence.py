"""Numerical accuracy checks separate from frozen solver-trajectory identity."""
import numpy as np

from cabomb_model.v13_circulation import (
    v13_default_parameters, v13_simulate,
    v13_max_abs_carrier_residual, v13_max_abs_cargo_residual,
)


def test_default_solution_converges_to_independent_tight_solvers():
    params=v13_default_parameters()
    default=v13_simulate(params,design='ocn-ysa',size='100nm')
    tight_lsoda=v13_simulate(params,design='ocn-ysa',size='100nm',rtol=1e-11,atol=1e-13)
    tight_radau=v13_simulate(params,design='ocn-ysa',size='100nm',method='Radau',rtol=1e-11,atol=1e-13)
    dose=params['D0']
    assert dose>0
    for result in [default,tight_lsoda,tight_radau]:
        assert result.success
        np.testing.assert_array_equal(result.t,default.t)
        assert np.isfinite(result.y).all()
        assert result.y.min()>=-1e-8*dose
        assert v13_max_abs_carrier_residual(result,params)<=1e-7*dose
        assert v13_max_abs_cargo_residual(result,params)<=1e-7*dose
    # Normalized whole-trajectory acceptance budgets, rather than a promise
    # that requested local solver tolerances are exact global-error bounds.
    assert np.max(np.abs(default.y-tight_radau.y))/dose<=1e-7
    assert np.max(np.abs(tight_lsoda.y-tight_radau.y))/dose<=1e-9
