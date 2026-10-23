"""Portable forecasting commands with separate fitting and late evaluation."""
from __future__ import annotations

import argparse
import hashlib
import json
from pathlib import Path

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

from cabomb_release import methods
from cabomb_release.bridge import run_bridge

KEY = ['study_id', 'curve_id', 'time_h']
ALL_METHODS = methods.FITTED + methods.NAIVE + ['old4SSE','old4Rolling','expandedRolling']

def sha(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()

def save_json(path: Path, payload: dict) -> None:
    path.write_text(json.dumps(payload,ensure_ascii=False,indent=2,allow_nan=False)+'\n',encoding='utf-8')

def read_series(path: Path, values: bool) -> pd.DataFrame:
    frame = pd.read_csv(path, dtype={'study_id':str,'curve_id':str})
    required = KEY + (['release_percent'] if values else [])
    if not set(required).issubset(frame.columns) or frame.empty:
        raise ValueError(f'{path.name}: required columns {required}, with at least one row')
    # A future-time file cannot carry observations into the forecasting step.
    if not values and set(frame.columns) != set(KEY):
        raise ValueError('future CSV must contain only study_id, curve_id, time_h')
    if frame[required].isna().any().any() or frame.duplicated(KEY).any():
        raise ValueError('missing values or duplicate study/curve/time keys')
    numeric = ['time_h'] + (['release_percent'] if values else [])
    for column in numeric:
        frame[column] = pd.to_numeric(frame[column],errors='raise')
        if not np.isfinite(frame[column]).all():
            raise ValueError('nonfinite numeric values')
    if (frame.time_h < 0).any():
        raise ValueError('negative time')
    if values and not frame.release_percent.between(0,100).all():
        raise ValueError('release_percent must be in [0,100]')
    return frame.sort_values(KEY).reset_index(drop=True)

def forecast_files(early_path: Path, future_path: Path, output: Path) -> None:
    early, future = read_series(early_path,True), read_series(future_path,False)
    groups = ['study_id','curve_id']
    if set(map(tuple,early[groups].drop_duplicates().values)) != set(map(tuple,future[groups].drop_duplicates().values)):
        raise ValueError('early and future curves must match exactly')
    # Validate the whole panel before any fitting or output creation.
    for identity, frame in early.groupby(groups,sort=True):
        times, values = frame.time_h.to_numpy(), frame.release_percent.to_numpy()
        query = future[(future.study_id==identity[0]) & (future.curve_id==identity[1])]
        if len(frame)<5 or times[0]!=0 or values[0]!=0 or np.any(np.diff(times)<=0) or np.any(query.time_h<=times[-1]):
            raise ValueError('each curve needs structural (0,0), four positive early times, and strictly later future times')
    output.mkdir(parents=True,exist_ok=False)
    rows, decisions, folds = [], [], []
    with (output/'optimization_diagnostics.jsonl').open('w',encoding='utf-8') as diagnostics:
        def sink(payload):
            diagnostics.write(json.dumps(payload,ensure_ascii=False,allow_nan=False)+'\n')
        for (study,curve), frame in early.groupby(groups,sort=True):
            query=future[(future.study_id==study)&(future.curve_id==curve)].time_h.to_numpy()
            result=methods.forecast(curve,frame.time_h.to_numpy(),frame.release_percent.to_numpy(),query,sink)
            for method, predictions in result['predictions'].items():
                for t,p in zip(query,predictions):
                    rows.append(dict(study_id=study,curve_id=curve,time_h=t,method=method,predicted_percent=p))
            decisions.append({'study_id':study,'curve_id':curve,'selectors':result['selectors']})
            folds.extend(dict(study_id=study,**row) for row in result['fold_rows'])
    pd.DataFrame(rows).to_csv(output/'predictions.csv',index=False,float_format='%.17g')
    pd.DataFrame(folds).to_csv(output/'rolling_folds.csv',index=False,float_format='%.17g')
    save_json(output/'decisions.json',{'curves':decisions})
    resource=Path(methods.__file__).parent/'resources/frozen_config.json'
    save_json(output/'forecast_lock.json',{
        'schema_version':1,'early_sha256':sha(early_path),'future_sha256':sha(future_path),
        'method_source_sha256':sha(Path(methods.__file__)),'protocol_sha256':sha(resource),
        'artifact_sha256':{name:sha(output/name) for name in [
            'predictions.csv','decisions.json','rolling_folds.csv','optimization_diagnostics.jsonl']},
        'late_observations_accessed':False,
        'meaning':'Provenance and integrity record; not cryptographic proof of blinding or an external preregistration.'})

def evaluate(run: Path, late_path: Path) -> None:
    lock=json.loads((run/'forecast_lock.json').read_text(encoding='utf-8'))
    for name,expected in lock['artifact_sha256'].items():
        if sha(run/name)!=expected:
            raise ValueError(f'forecast artifact changed: {name}')
    late=read_series(late_path,True)
    pred=pd.read_csv(run/'predictions.csv',dtype={'study_id':str,'curve_id':str})
    if pred.duplicated(KEY+['method']).any() or set(pred.method)!=set(ALL_METHODS):
        raise ValueError('duplicate or missing forecast methods')
    expected=set(map(tuple,late[KEY].values))
    for method,frame in pred.groupby('method'):
        if set(map(tuple,frame[KEY].values))!=expected:
            raise ValueError(f'late times do not exactly match {method} forecasts')
    joined=pred.merge(late[KEY+['release_percent']],on=KEY,validate='many_to_one')
    metrics=[]
    for (study,curve,method),frame in joined.groupby(['study_id','curve_id','method']):
        valid=np.isfinite(frame.predicted_percent).all()
        rmse=float(np.sqrt(np.mean((frame.predicted_percent-frame.release_percent)**2))) if valid else None
        metrics.append(dict(study_id=study,curve_id=curve,method=method,late_rmse_pp=rmse,
                            unresolved=not valid,failure_gt10pp_or_unresolved=not valid or rmse>10))
    out=run/'evaluation'
    out.mkdir(exist_ok=False)
    pd.DataFrame(metrics).to_csv(out/'curve_metrics.csv',index=False,float_format='%.17g')
    save_json(out/'evaluation_manifest.json',{'late_sha256':sha(late_path),'forecast_lock_sha256':sha(run/'forecast_lock.json'),
        'description':'Post-forecast scoring; no refitting or model selection on late outcomes.'})

def synthetic(output: Path) -> None:
    output.mkdir(parents=True,exist_ok=False)
    inputs=output/'inputs'
    inputs.mkdir()
    time=np.array([0.,1.,2.,4.,6.,8.,12.,24.,48.,72.])
    # Deterministic mathematical example, deliberately not experimental data.
    y=75*(-np.expm1(-(time/25)**.8))
    series=pd.DataFrame({'study_id':'synthetic','curve_id':'synthetic_weibull','time_h':time,'release_percent':y})
    early=series.iloc[:7]
    late=series.iloc[7:]
    early.to_csv(inputs/'early.csv',index=False)
    late[KEY].to_csv(inputs/'future_times.csv',index=False)
    late.to_csv(inputs/'late.csv',index=False)
    save_json(inputs/'generation.json',{'status':'synthetic, noise-free, not biological validation',
        'equation':'y(t)=75*(1-exp(-(t/25)^0.8))','time_unit':'hour','value_unit':'percent'})
    run=output/'forecast'
    forecast_files(inputs/'early.csv',inputs/'future_times.csv',run)
    evaluate(run,inputs/'late.csv')
    pred=pd.read_csv(run/'predictions.csv')
    fig,ax=plt.subplots(figsize=(7,4.3))
    ax.plot(time,y,'k.-',label='Synthetic reference')
    for method in ['old4SSE','expandedRolling','hold_last']:
        p=pred[pred.method==method]
        ax.plot(p.time_h,p.predicted_percent,'o--',label=method)
    ax.axvline(12,color='gray',linestyle=':',label='End of early window')
    ax.set(xlabel='Time (h)',ylabel='Cumulative release (%)',ylim=(0,100),
           title='Early-only forecasts on a noise-free synthetic curve')
    ax.legend(fontsize=8)
    fig.tight_layout()
    fig.savefig(output/'synthetic_forecast.png',dpi=170)
    plt.close(fig)
    (output/'caption.txt').write_text('Mathematically generated series. Late points are scored only after forecasts are saved. Successful recovery here does not establish performance on experiments.\n',encoding='utf-8')

def main(argv: list[str] | None=None) -> int:
    parser=argparse.ArgumentParser(prog='cabomb-release')
    commands=parser.add_subparsers(dest='command',required=True)
    for name in ['synthetic','bridge']:
        command=commands.add_parser(name)
        command.add_argument('--output-dir',type=Path,required=True)
    forecast=commands.add_parser('forecast')
    forecast.add_argument('--early',type=Path,required=True)
    forecast.add_argument('--future',type=Path,required=True)
    forecast.add_argument('--output-dir',type=Path,required=True)
    score=commands.add_parser('evaluate')
    score.add_argument('--run-dir',type=Path,required=True)
    score.add_argument('--late',type=Path,required=True)
    args=parser.parse_args(argv)
    if args.command=='synthetic': synthetic(args.output_dir)
    elif args.command=='bridge':
        run_bridge(args.output_dir)
        from cabomb_release.plotting import plot_bridge
        plot_bridge(args.output_dir)
    elif args.command=='forecast': forecast_files(args.early,args.future,args.output_dir)
    elif args.command=='evaluate': evaluate(args.run_dir,args.late)
    print('Completed:',args.command)
    return 0

if __name__=='__main__':
    raise SystemExit(main())
