"""Information boundaries and standalone resource behavior."""
import hashlib
import json
from pathlib import Path

import numpy as np
import pandas as pd
import pytest

from cabomb_model.parameters import default_registry_path
from cabomb_release import methods
from cabomb_release.cli import ALL_METHODS, evaluate, forecast_files, read_series


def test_registry_is_bundled_and_matches_readable_copy():
    root = Path(__file__).resolve().parents[1]
    bundled = default_registry_path()
    assert bundled.parent.name == 'resources'
    assert bundled.read_bytes() == (root/'data/processed/parameters_master.csv').read_bytes()


def test_future_file_rejects_embedded_outcomes(tmp_path):
    path=tmp_path/'future.csv'
    pd.DataFrame([dict(study_id='S',curve_id='C',time_h=3,release_percent=80)]).to_csv(path,index=False)
    with pytest.raises(ValueError,match='only'):
        read_series(path,False)


@pytest.mark.parametrize('bad',[np.nan,np.inf,-1.,101.])
def test_bad_observations_are_rejected_before_fitting(tmp_path,bad):
    path=tmp_path/'early.csv'
    pd.DataFrame(dict(study_id=['S'],curve_id=['C'],time_h=[0],release_percent=[bad])).to_csv(path,index=False)
    with pytest.raises(ValueError):
        read_series(path,True)


def write_locked_predictions(run, future_time=24):
    run.mkdir()
    pd.DataFrame([dict(study_id='S',curve_id='C',time_h=future_time,method=m,predicted_percent=20.)
                  for m in ALL_METHODS]).to_csv(run/'predictions.csv',index=False)
    digest=hashlib.sha256((run/'predictions.csv').read_bytes()).hexdigest()
    (run/'forecast_lock.json').write_text(json.dumps({'artifact_sha256':{'predictions.csv':digest}}))


def test_evaluation_does_not_refit_and_checks_exact_late_times(tmp_path,monkeypatch):
    def forbidden(*args,**kwargs):
        raise AssertionError('evaluation attempted refitting')
    monkeypatch.setattr(methods,'fit',forbidden)
    run=tmp_path/'run'
    write_locked_predictions(run)
    late=tmp_path/'late.csv'
    pd.DataFrame([dict(study_id='S',curve_id='C',time_h=25,release_percent=30.)]).to_csv(late,index=False)
    with pytest.raises(ValueError,match='exactly'):
        evaluate(run,late)
    assert not (run/'evaluation').exists()
    pd.DataFrame([dict(study_id='S',curve_id='C',time_h=24,release_percent=30.)]).to_csv(late,index=False)
    evaluate(run,late)
    frame=pd.read_csv(run/'evaluation/curve_metrics.csv')
    assert len(frame)==13
    np.testing.assert_allclose(frame.late_rmse_pp,10.)


def test_changed_forecasts_cannot_be_evaluated_as_original(tmp_path):
    run=tmp_path/'run'
    write_locked_predictions(run)
    with (run/'predictions.csv').open('a') as handle:
        handle.write('\n')
    with pytest.raises(ValueError,match='changed'):
        evaluate(run,tmp_path/'unused_late.csv')


def test_invalid_prefix_does_not_create_output_or_call_fitting(tmp_path,monkeypatch):
    def forbidden(*args,**kwargs):
        raise AssertionError('invalid input reached fitting')
    monkeypatch.setattr(methods,'forecast',forbidden)
    early=tmp_path/'early.csv'
    future=tmp_path/'future.csv'
    pd.DataFrame(dict(study_id=['S']*5,curve_id=['C']*5,time_h=[0,1,2,3,4],
                      release_percent=[1,2,3,4,5])).to_csv(early,index=False)
    pd.DataFrame([dict(study_id='S',curve_id='C',time_h=10)]).to_csv(future,index=False)
    with pytest.raises(ValueError,match='structural'):
        forecast_files(early,future,tmp_path/'output')
    assert not (tmp_path/'output').exists()
