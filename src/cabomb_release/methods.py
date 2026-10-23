"""Empirical bounded percentage forecasts; no mechanistic parameter inference."""
from __future__ import annotations

import hashlib
import json
from pathlib import Path
from typing import Any

import numpy as np
from scipy.optimize import least_squares

CONFIG = json.loads((Path(__file__).resolve().parent / 'resources/frozen_config.json').read_text(encoding='utf-8'))
FITTED = CONFIG['fitted_methods']
NAIVE = CONFIG['naive_methods']
OLD = CONFIG['old_methods']


def decode(family: str, u: np.ndarray) -> np.ndarray:
    """Convert normalized transformed coordinates using the frozen registry."""
    u = np.asarray(u, dtype=float)
    values = []
    for index, name in enumerate(CONFIG['parameter_order'][family]):
        lo, hi = CONFIG['bounds'][name]
        z = u[..., index]
        mode = CONFIG['transforms'][name]
        if mode == 'log':
            value = np.exp(np.log(lo) + z * np.log(hi / lo))
        elif mode == 'log1p':
            value = np.expm1(np.log1p(lo) + z * (np.log1p(hi)-np.log1p(lo)))
        else:
            value = lo + z * (hi-lo)
        values.append(value)
    return np.stack(values, axis=-1)


def predict(family: str, x: np.ndarray, p: np.ndarray) -> np.ndarray:
    """Return a parameter-row by time-column prediction array."""
    x = np.asarray(x, dtype=float)
    p = np.atleast_2d(p)
    if np.any(x < 0):
        raise ValueError('negative time')
    if family in ('exp100', 'expA'):
        amplitude = 100 if family == 'exp100' else p[:, 1, None]
        return amplitude * (-np.expm1(-p[:, 0, None] * x))
    if family in ('weibull100', 'weibullA'):
        amplitude = 100 if family == 'weibull100' else p[:, 2, None]
        return amplitude * (-np.expm1(-np.power(x/p[:, 0, None], p[:, 1, None])))
    if family == 'capped_linear':
        return np.minimum(100, p[:, 0, None]*x)
    if family == 'higuchi_cap':
        return np.minimum(100, p[:, 0, None]*np.sqrt(x))
    if family == 'kp_cap':
        return np.minimum(100, p[:, 0, None]*np.power(x, p[:, 1, None]))
    if family == 'hc_cap':
        return 100*(1-np.maximum(0, 1-p[:, 0, None]*x)**3)
    raise ValueError(family)


def fit(family: str, x: np.ndarray, y: np.ndarray, curve_id: str,
        fold_suffix: str = '') -> dict[str, Any]:
    """Retain all 33 terminal optimizer states, including unsuccessful optima."""
    x, y = np.asarray(x, dtype=float), np.asarray(y, dtype=float)
    dim = len(CONFIG['parameter_order'][family])
    seed = (CONFIG['seed']+int(hashlib.sha256((curve_id+family+fold_suffix).encode()).hexdigest()[:8],16)) % 2**32
    rng = np.random.default_rng(seed)
    starts = np.vstack((np.full((1,dim),.5),rng.uniform(0,1,(32,dim))))
    runs = []
    best_jac = {}
    for start_index, start in enumerate(starts):
        try:
            result = least_squares(lambda u: predict(family,x,decode(family,u))[0]-y,
                                   start,bounds=(0,1),ftol=1e-11,xtol=1e-11,
                                   gtol=1e-11,max_nfev=20000)
            params = decode(family,result.x)
            rmse = float(np.sqrt(np.mean((predict(family,x,params)[0]-y)**2)))
            singular = np.linalg.svd(result.jac, compute_uv=False)
            rank = int(np.linalg.matrix_rank(result.jac))
            row = {'start_index':start_index,'initial_u':start.tolist(),
                   'terminal_u':result.x.tolist(),'parameters_scaled':params.tolist(),
                   'rmse_pp':rmse,'sse':float(rmse**2*len(y)),
                   'success':bool(result.success),'status':int(result.status),
                   'message':str(result.message),'nfev':int(result.nfev),
                   'optimality':float(result.optimality),
                   'bound_indices':np.flatnonzero(np.minimum(result.x,1-result.x)<=1e-5).tolist(),
                   'jacobian_rank':rank,'jacobian_singular_values':singular.tolist()}
            runs.append(row)
        except (ValueError,RuntimeError,FloatingPointError,np.linalg.LinAlgError) as exc:
            runs.append({'start_index':start_index,'initial_u':start.tolist(),
                         'success':False,'error':str(exc)})
    valid = [r for r in runs if r['success'] and np.isfinite(r.get('rmse_pp',np.nan))]
    info = {'family':family,'curve_id':curve_id,'fold_suffix':fold_suffix,'seed':int(seed),
            'n_positive':len(y),'n_parameters':dim,'nominal_n_minus_p':len(y)-dim,'runs':runs}
    if not valid:
        return {**info,'unresolved':True,'usable':False,'error':'no successful finite fit'}
    best = min(valid,key=lambda r:r['rmse_pp'])
    unresolved = any(not r['success'] and r.get('rmse_pp',np.inf)<best['rmse_pp']-1e-6 for r in runs)
    return {**info,**{k:v for k,v in best.items() if k!='initial_u'},
            'unresolved':unresolved,'usable':not unresolved,
            'n_failed_starts':sum(not r['success'] for r in runs)}


def naive(family: str, x: np.ndarray, y: np.ndarray, future_x: np.ndarray) -> np.ndarray:
    """Forecast with fixed persistence or last-three OLS slope, anchored last."""
    x,y,future_x=map(lambda z:np.asarray(z,dtype=float),(x,y,future_x))
    if family=='hold_last':
        return np.full_like(future_x,np.clip(y[-1],0,100))
    if family=='last3OLS_slope_anchored_last':
        a,b=x[-3:],y[-3:]
        slope=max(0,float(np.dot(a-a.mean(),b-b.mean())/np.dot(a-a.mean(),a-a.mean())))
        return np.clip(y[-1]+slope*(future_x-x[-1]),0,100)
    raise ValueError(family)


def rank_scores(scores: dict[str,float], candidates: list[str]) -> list[str]:
    """Resolve score ties in the frozen family order, within the frozen tolerance."""
    remaining=[f for f in candidates if np.isfinite(scores.get(f,np.nan))]
    ranked=[]
    while remaining:
        minimum=min(scores[f] for f in remaining)
        winner=next(f for f in remaining if scores[f]<=minimum+1e-8)
        ranked.append(winner)
        remaining.remove(winner)
    return ranked


def select_rolling(candidates: list[str], fold_errors: dict[str,list[float | None]],
                   full_fits: dict[str,dict], fallback: str | None,
                   n_folds: int) -> dict:
    """Use complete common folds only; retain ineligibility and fallback reasons."""
    scores,excluded={},{}
    for family in candidates:
        if n_folds==0:
            excluded[family]='no rolling folds; score undefined'
            continue
        errors=fold_errors.get(family,[])
        if len(errors)!=n_folds or any(e is None or not np.isfinite(e) for e in errors):
            excluded[family]='failed or unresolved common fold'
        else:
            scores[family]=float(np.sqrt(np.mean(np.square(errors)))) if n_folds else float('nan')
    ranking=rank_scores(scores,candidates)
    skipped=[]
    if n_folds>=2:
        for family in ranking:
            if family in NAIVE or full_fits.get(family,{}).get('usable',False):
                return {'selected_family':family,'fallback':False,'fallback_reason':None,
                        'n_folds':n_folds,'cv_scores_pp':scores,'excluded':excluded,
                        'ranking':ranking,'skipped_full_fit_failures':skipped}
            skipped.append(family)
    reason='fewer than two common folds' if n_folds<2 else 'no eligible usable full fit'
    return {'selected_family':fallback,'fallback':True,'fallback_reason':reason,
            'n_folds':n_folds,'cv_scores_pp':scores,'excluded':excluded,
            'ranking':ranking,'skipped_full_fit_failures':skipped,
            'unresolved':fallback is None}


def forecast(curve_id: str, time_h: np.ndarray, values: np.ndarray,
             future_time_h: np.ndarray, diagnostic_sink=None) -> dict:
    """Only early values enter this interface; late outcomes are unavailable."""
    time_h,values,future_time_h=map(lambda z:np.asarray(z,dtype=float),(time_h,values,future_time_h))
    if len(time_h)!=len(values) or time_h[0]!=0 or values[0]!=0 or not np.all(np.diff(time_h)>0):
        raise ValueError('requires structural zero and strictly increasing times')
    T=float(time_h[-1])
    x,y=time_h[1:]/T,values[1:]
    future_x=future_time_h/T
    if len(x)<4 or np.any(future_x<=1):
        raise ValueError('insufficient early values or invalid late time')
    full={}
    for family in FITTED:
        full[family]=fit(family,x,y,curve_id)
        if diagnostic_sink: diagnostic_sink(full[family])
    old_scores={f:full[f].get('rmse_pp',np.nan) for f in OLD}
    # Original pipeline requires all four valid optima for its selection.
    old_rank=rank_scores(old_scores,OLD) if all(full[f].get('usable',False) for f in OLD) else []
    old_selected=old_rank[0] if old_rank else None
    fold_errors={f:[] for f in FITTED+NAIVE}
    fold_rows=[]
    for n_prefix in range(4,len(x)):
        for family in FITTED+NAIVE:
            if family in NAIVE:
                predicted=float(naive(family,x[:n_prefix],y[:n_prefix],x[n_prefix:n_prefix+1])[0])
                valid=True
            else:
                fitted=fit(family,x[:n_prefix],y[:n_prefix],curve_id,f':prefix{n_prefix}')
                if diagnostic_sink: diagnostic_sink(fitted)
                valid=fitted.get('usable',False)
                predicted=float(predict(family,x[n_prefix:n_prefix+1],np.asarray(fitted['parameters_scaled']))[0,0]) if valid else None
            error=predicted-float(y[n_prefix]) if predicted is not None and np.isfinite(predicted) else None
            fold_errors[family].append(error)
            fold_rows.append({'curve_id':curve_id,'method':family,'n_prefix':n_prefix,
                              'heldout_time_h':float(time_h[n_prefix+1]),'observed_early_percent':float(y[n_prefix]),
                              'predicted_percent':predicted,'error_pp':error,'usable':bool(valid and error is not None)})
    selectors={'old4SSE':{'selected_family':old_selected,'fallback':False,'n_folds':0},
               'old4Rolling':select_rolling(OLD,fold_errors,full,old_selected,max(0,len(x)-4)),
               'expandedRolling':select_rolling(FITTED+NAIVE,fold_errors,full,old_selected,max(0,len(x)-4))}
    predictions={}
    for method in FITTED+NAIVE+list(selectors):
        family=selectors[method]['selected_family'] if method in selectors else method
        if family is None or (family in FITTED and not full[family].get('usable',False)):
            predictions[method]=[None]*len(future_x)
        elif family in NAIVE:
            predictions[method]=naive(family,x,y,future_x).tolist()
        else:
            predictions[method]=predict(family,future_x,np.asarray(full[family]['parameters_scaled']))[0].tolist()
    return {'curve_id':curve_id,'T_early_end_h':T,'full_fits':full,'fold_rows':fold_rows,
            'selectors':selectors,'predictions':predictions}
