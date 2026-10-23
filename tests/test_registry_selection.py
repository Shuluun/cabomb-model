"""User-selected registries must affect runs and remain traceable."""
import csv
import hashlib
import json

import pytest

from cabomb_model.cli import main
from cabomb_model.parameters import default_registry_path, load_registry
from cabomb_model.resource_paths import parameter_registry
from cabomb_model.v13_circulation import v13_default_parameters
from cabomb_model.v13_analysis import run_v13_design_comparison


def custom_registry(tmp_path):
    rows=load_registry()
    for row in rows:
        if row['parameter_symbol']=='D0':
            row.update(baseline_value='2',lower_bound='2',upper_bound='2',
                       source_value_text='2 normalized carrier; software test condition')
    path=tmp_path/'custom.csv'
    with path.open('w',encoding='utf-8',newline='') as handle:
        writer=csv.DictWriter(handle,fieldnames=list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)
    return path


def test_custom_registry_controls_cli_and_is_snapshotted(tmp_path):
    baseline=v13_default_parameters()['D0']
    path=custom_registry(tmp_path)
    output=tmp_path/'run'
    assert main(['run-baseline','--model-version','v1.3','--registry',str(path),
                 '--output-dir',str(output)])==0
    with (output/'tables/baseline_timecourse.csv').open(encoding='utf-8') as handle:
        first=next(csv.DictReader(handle))
    assert float(first['A_central_venous'])==2
    assert (output/'registry_used.csv').read_bytes()==path.read_bytes()
    record=json.loads((output/'registry_source.json').read_text(encoding='utf-8'))
    assert record['selection']=='explicit --registry'
    assert record['sha256']==hashlib.sha256(path.read_bytes()).hexdigest()
    assert v13_default_parameters()['D0']==baseline


def test_external_registry_manifest_hash_matches_actual_input(tmp_path):
    path=custom_registry(tmp_path)
    output=tmp_path/'comparison'
    with parameter_registry(path):
        run_v13_design_comparison(output)
    record=json.loads((output/'run_manifest.json').read_text(encoding='utf-8'))
    assert record['parameters']['D0']==2
    assert record['source_sha256']['parameter_registry.csv']==hashlib.sha256(path.read_bytes()).hexdigest()


def test_selection_restored_after_exception(tmp_path):
    original=default_registry_path()
    path=custom_registry(tmp_path)
    with pytest.raises(RuntimeError):
        with parameter_registry(path):
            assert default_registry_path()==path.resolve()
            raise RuntimeError('test interruption')
    assert default_registry_path()==original


def test_nonexistent_registry_is_an_error_not_a_fallback(tmp_path):
    with pytest.raises(FileNotFoundError):
        main(['validate-registry','--registry',str(tmp_path/'missing.csv')])
