"""Fixed V1.3 input with age-since-internalization response scenarios."""
from __future__ import annotations

import dataclasses
import json
import math
import sys

import numpy as np
import pandas as pd
from scipy.integrate import solve_ivp

from pathlib import Path
import hashlib
from cabomb_release.methods import CONFIG
from cabomb_model.resource_paths import resource_path

def sha(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()

def write_json(path, payload):
    path.write_text(json.dumps(payload, ensure_ascii=False, indent=2, allow_nan=False)+"\n", encoding="utf-8")
from cabomb_model.fluxes import degradation_rate
from cabomb_model.model_modes import V1_3_MODES
from cabomb_model.scenarios import make_scenario
from cabomb_model.v13_circulation import v13_default_parameters, v13_fluxes, v13_rhs
from cabomb_model.v13_states import V13_STATE_INDEX, V13_STATES, v13_initial_state


def cdf(age: np.ndarray, beta: float, mean_h: float=4.0) -> np.ndarray:
    """Unit-platform Weibull cumulative response at a common mean wait."""
    age=np.asarray(age,dtype=float)
    if np.any(age<0): raise ValueError('negative age')
    tau=mean_h/math.gamma(1+1/beta)
    return -np.expm1(-(age/tau)**beta)


def cohort_response(times: np.ndarray, flux: np.ndarray, beta: float) -> tuple[np.ndarray,np.ndarray,np.ndarray]:
    """Trapezoid bin mass at midpoint age, never an injection-clock hazard."""
    times,flux=map(lambda a:np.asarray(a,dtype=float),(times,flux))
    if len(times)!=len(flux) or not np.isfinite(times).all() or not np.isfinite(flux).all() or np.any(np.diff(times)<=0) or np.any(flux<0):
        raise ValueError('invalid times or negative input')
    mass=.5*(flux[:-1]+flux[1:])*np.diff(times)
    midpoint=.5*(times[:-1]+times[1:])
    total=np.r_[0,np.cumsum(mass)]
    gross=np.zeros_like(times)
    for i in range(1,len(times)):
        gross[i]=np.dot(mass[:i],cdf(times[i]-midpoint[:i],beta))
    return gross,total-gross,total


def run_bridge(output_dir: str | Path) -> None:
    """Generate fixed model fluxes, independent response ledgers and checks."""
    output = Path(output_dir)
    output.mkdir(parents=True, exist_ok=False)
    params=v13_default_parameters()
    scenario=make_scenario(params,epr='medium',mps='medium')
    if scenario.intracellular_pH!=5.5 or degradation_rate(scenario.intracellular_pH,params)!=.25:
        raise RuntimeError('current default scenario no longer matches frozen reference')
    spec=CONFIG['bridge']
    # Output directory is explicit and must be new.
    write_json(output/'bridge_inputs.json',{'params':params,'scenario':dataclasses.asdict(scenario),
               'modes':dataclasses.asdict(V1_3_MODES),'size':'200nm','kernel_config':spec,
                'source_hashes':{
                   'bridge.py':sha(Path(__file__)),
                   'methods.py':sha(Path(__file__).with_name('methods.py')),
                   'frozen_config.json':sha(Path(__file__).parent/'resources/frozen_config.json'),
                   'parameters_master.csv':sha(resource_path('parameters_master.csv')),
                   **{'cabomb_model/'+p.name:sha(p) for p in sorted(resource_path('parameters_master.csv').parent.parent.glob('*.py'))}}})
    full_rows,metrics,checks,input_rows=[],[],[],[]
    for design in spec['designs']:
        n=len(V13_STATES)
        def augmented(t,z):
            y=z[:n]
            f=v13_fluxes(y,params,design,scenario,'200nm',V1_3_MODES.targeting_model)
            J=f['internalization']+f['nonspecific_internalization']
            return np.r_[v13_rhs(t,y,params,design,scenario,'200nm',V1_3_MODES),J,f['deg_int']]
        solution=solve_ivp(augmented,(0,72),np.r_[v13_initial_state(params['D0']),0,0],
                           method='LSODA',rtol=1e-9,atol=1e-11,dense_output=True)
        if not solution.success: raise RuntimeError(solution.message)
        results={}
        for dt in spec['quadrature_steps_h']:
            times=np.linspace(0,72,int(round(72/dt))+1)
            states=solution.sol(times)
            if not np.isfinite(states).all(): raise RuntimeError('nonfinite source state')
            fl=np.array([v13_fluxes(states[:n,i],params,design,scenario,'200nm',V1_3_MODES.targeting_model)['internalization']+
                         v13_fluxes(states[:n,i],params,design,scenario,'200nm',V1_3_MODES.targeting_model)['nonspecific_internalization'] for i in range(len(times))])
            minimum=float(fl.min())
            if not np.isfinite(fl).all(): raise RuntimeError('nonfinite source flux')
            negatives=int(np.sum(fl<0))
            if minimum < -1e-10: raise RuntimeError('negative flux beyond tolerance')
            fl=np.where(fl<0,0,fl)
            for beta in spec['beta_scenarios']:
                gross,pending,total=cohort_response(times,fl,beta)
                residual=float(np.max(np.abs(gross+pending-total)))
                identity=states[n]-states[V13_STATE_INDEX['A_internalized']]
                first_order_error=float(np.max(np.abs(gross-identity))) if beta==1 else None
                integral_error=float(np.max(np.abs(total-states[n])))
                reference_error=float(np.max(np.abs(identity-states[n+1])))
                check={'design':design,'dt_h':dt,'beta':beta,'solver_success':bool(solution.success),
                       'solver_nfev':solution.nfev,'minimum_source_state':float(states[:n].min()),
                       'minimum_raw_J':minimum,'negative_flux_clipped_count':negatives,
                       'minimum_pending':float(pending.min()),'minimum_gross':float(gross.min()),
                       'mass_ledger_max_abs_residual':residual,'input_integral_max_abs_error':integral_error,
                       'first_order_identity_max_abs_error':first_order_error,
                       'solver_identity_vs_gross_deg_max_abs_error':reference_error,
                       'gross_monotone':bool(np.all(np.diff(gross)>=-1e-12)),
                       'gross_not_above_input':bool(np.all(gross<=total+1e-10)),
                       'pass':bool(states[:n].min()>=-1e-10 and residual<=1e-10 and pending.min()>=-1e-10 and gross.min()>=0 and np.all(np.diff(gross)>=-1e-12) and np.all(gross<=total+1e-10) and integral_error<=1e-4 and (first_order_error is None or first_order_error<=1e-4) and reference_error<=1e-8)}
                checks.append(check)
                results[(dt,beta)]=(times,gross,pending,total)
                if dt==spec['quadrature_steps_h'][-1]:
                    for i,t in enumerate(times):
                        full_rows.append({'design':design,'time_h':t,'beta':beta,'mean_wait_h':4,
                                          'tau_h':4/math.gamma(1+1/beta),'gross_response':gross[i],
                                          'pending_response':pending[i],'cumulative_internalized_input':total[i]})
                    for t in spec['output_hours']:
                        i=int(round(t/dt))
                        metrics.append({'design':design,'time_h':t,'beta':beta,'gross_response':gross[i],
                                        'pending_response':pending[i],'cumulative_internalized_input':total[i]})
            if dt==spec['quadrature_steps_h'][-1]:
                for i,t in enumerate(times):
                    input_rows.append({'design':design,'time_h':t,'J_total':fl[i],
                                       'cumulative_input_solver':states[n,i],
                                       'A_internalized_original':states[V13_STATE_INDEX['A_internalized'],i],
                                       'cumulative_degradation_solver':states[n+1,i]})
        for beta in spec['beta_scenarios']:
            coarse,mid,fine=(results[(dt,beta)] for dt in spec['quadrature_steps_h'])
            delta_coarse=float(np.max(np.abs(coarse[1]-fine[1][::4])))
            delta_mid=float(np.max(np.abs(mid[1]-fine[1][::2])))
            checks.append({'design':design,'beta':beta,'check':'quadrature_convergence',
                           'coarse_vs_fine_max_abs':delta_coarse,'mid_vs_fine_max_abs':delta_mid,
                           'pass':bool(delta_mid<=1e-4)})
        print('Bridge complete:',design,flush=True)
    metrics=pd.DataFrame(metrics)
    metrics['rank_desc']=metrics.groupby(['time_h','beta']).gross_response.rank(method='min',ascending=False)
    metrics.to_csv(output/'bridge_metrics.csv',index=False,float_format='%.17g')
    pd.DataFrame(full_rows).to_csv(output/'bridge_curves.csv',index=False,float_format='%.17g')
    pd.DataFrame(input_rows).to_csv(output/'bridge_fixed_inputs.csv',index=False,float_format='%.17g')
    write_json(output/'bridge_checks.json',{'checks':checks,'all_pass':all(c['pass'] for c in checks),
               'warning':'Assumed unit response, no calcium fitted parameter transfer; not measured cargo or free calcium.'})
    print('Bridge all checks:',all(c['pass'] for c in checks),flush=True)


if __name__=='__main__':
    raise SystemExit("Use cabomb-release bridge --output-dir PATH")
