"""Report numerical solver and frozen-trajectory agreement without changing tolerances."""
from pathlib import Path
from importlib.metadata import version
from importlib import import_module
import argparse
import csv
import hashlib
import json
import platform
import sys

import numpy as np

from cabomb_model.v13_circulation import (
    v13_default_parameters, v13_simulate,
    v13_max_abs_carrier_residual, v13_max_abs_cargo_residual,
)
from cabomb_model.parameters import default_registry_path

ROOT=Path(__file__).resolve().parents[1]
parser=argparse.ArgumentParser(description=__doc__)
parser.add_argument('--strict-history',action='store_true',
                    help='Also require the original rtol=1e-9/atol=1e-11 frozen-trajectory replay.')
args=parser.parse_args()
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
reference_info=json.loads((ROOT/'research/numerical_reference.json').read_text(encoding='utf-8'))
result['historical_source_identity']={
    name:hashlib.sha256((ROOT/name).read_bytes()).hexdigest()==expected
    for name,expected in reference_info['expected_source_sha256'].items()
}
result['historical_fixture_identity']=hashlib.sha256((ROOT/reference_info['fixture']).read_bytes()).hexdigest()==reference_info['fixture_sha256']
executed_identity={}
executed_paths={}
for name,expected in reference_info['expected_source_sha256'].items():
    if name.startswith('src/cabomb_model/'):
        module_name='cabomb_model.'+Path(name).stem
        actual_path=Path(import_module(module_name).__file__)
    else:
        actual_path=default_registry_path()
    executed_paths[name]=str(actual_path)
    executed_identity[name]=hashlib.sha256(actual_path.read_bytes()).hexdigest()==expected
result['executed_source_paths']=executed_paths
result['executed_source_identity']=executed_identity
dose=params['D0']
checks={
    'default_vs_tight_radau':result['default_vs_tight_radau']['maximum_absolute_difference']/dose<=1e-7,
    'tight_solver_agreement':result['tight_lsoda_vs_radau']['maximum_absolute_difference']/dose<=1e-9,
    'mass_balance':all(max(row.values())/dose<=1e-7 for row in result['mass_residuals'].values()),
    'nonnegativity':all(run.y.min()/dose>=-1e-8 for run in runs.values()),
}
result['numerical_acceptance_budgets_normalized']={'default_vs_tight_radau':1e-7,'tight_solver_agreement':1e-9}
result['numerical_checks']=checks
result['strict_historical_numerical_match']=result['historical_vs_default']['historical_strict_allclose']
result['strict_historical_replay_pass']=bool(
    result['strict_historical_numerical_match'] and result['historical_fixture_identity']
    and all(result['historical_source_identity'].values()) and all(executed_identity.values())
)
print('SOLVER_DIAGNOSTICS_JSON='+json.dumps(result,separators=(',',':'),allow_nan=False))
sys.exit(0 if all(checks.values()) and (not args.strict_history or result['strict_historical_replay_pass']) else 1)
