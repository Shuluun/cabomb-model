"""Report numerical solver and frozen-trajectory agreement without changing tolerances."""
from pathlib import Path
from importlib.metadata import version
import csv
import json
import platform

import numpy as np

from cabomb_model.v13_circulation import (
    v13_default_parameters, v13_simulate,
    v13_max_abs_carrier_residual, v13_max_abs_cargo_residual,
)

ROOT=Path(__file__).resolve().parents[1]
with (ROOT/'tests/fixtures/v13_baseline_timecourse.csv').open(encoding='utf-8',newline='') as handle:
    rows=list(csv.DictReader(handle))
params=v13_default_parameters()
runs={
    'default_lsoda':v13_simulate(params,design='ocn-ysa',size='100nm'),
    'tight_lsoda':v13_simulate(params,design='ocn-ysa',size='100nm',rtol=1e-11,atol=1e-13),
    'tight_radau':v13_simulate(params,design='ocn-ysa',size='100nm',method='Radau',rtol=1e-11,atol=1e-13),
}
assert all(r.success for r in runs.values())
baseline=runs['default_lsoda']
assert np.array_equal(baseline.t,np.array([float(r['time_h']) for r in rows]))
reference=np.array([[float(row[state]) for row in rows] for state in baseline.state_names])

def compare(a,b):
    difference=np.abs(a-b)
    return {'maximum_absolute_difference':float(difference.max()),
            'historical_strict_allclose':bool(np.allclose(a,b,rtol=1e-9,atol=1e-11)),
            'states':{name:{'maximum_absolute_difference':float(difference[i].max()),
                    'time_of_max_h':float(baseline.t[np.argmax(difference[i])]),
                    'strict_allclose':bool(np.allclose(a[i],b[i],rtol=1e-9,atol=1e-11))}
                    for i,name in enumerate(baseline.state_names)}}

result={
    'platform':platform.platform(),'python':platform.python_version(),
    'versions':{name:version(name) for name in ['numpy','scipy','pandas']},
    'default_tolerances':{'rtol':baseline.rtol,'atol':baseline.atol},
    'historical_vs_default':compare(reference,baseline.y),
    'historical_vs_tight_lsoda':compare(reference,runs['tight_lsoda'].y),
    'default_vs_tight_lsoda':compare(baseline.y,runs['tight_lsoda'].y),
    'default_vs_tight_radau':compare(baseline.y,runs['tight_radau'].y),
    'tight_lsoda_vs_radau':compare(runs['tight_lsoda'].y,runs['tight_radau'].y),
    'mass_residuals':{name:{'carrier':v13_max_abs_carrier_residual(run,params),
                    'cargo':v13_max_abs_cargo_residual(run,params)} for name,run in runs.items()},
}
print('SOLVER_DIAGNOSTICS_JSON='+json.dumps(result,separators=(',',':'),allow_nan=False))
