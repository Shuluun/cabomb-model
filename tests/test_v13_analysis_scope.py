"""Guard against inactive sweep dimensions and data-independent reporting."""
import csv
from pathlib import Path

import numpy as np

from cabomb_model.v13_analysis import (
    run_v13_feasibility_domain,
    run_v13_local_sensitivity,
    run_v13_monte_carlo,
)
from cabomb_model.v13_circulation import (
    v13_default_parameters, v13_detailed_metrics, v13_fluxes, v13_simulate,
)
from cabomb_model.v13_states import v13_initial_state, V13_STATE_INDEX
from cabomb_model.scenarios import make_scenario


def test_binding_coefficient_is_inactive_only_for_sequential_dual():
    params=v13_default_parameters()
    state=np.array(v13_initial_state(params['D0']))
    state[V13_STATE_INDEX['A_tumor_free']]=.1
    scenario=make_scenario(params)
    changed=dict(params,k_ysa_on=params['k_ysa_on']*2)
    for design in ['ysa','ocn-ysa']:
        a=v13_fluxes(state,params,design,scenario,'100nm','sequential_v1_2')
        b=v13_fluxes(state,changed,design,scenario,'100nm','sequential_v1_2')
        if design=='ysa':
            assert b['ysa_from_free']==2*a['ysa_from_free']>0
        else:
            assert a==b


def test_default_sensitivity_and_sampling_exclude_inactive_routes(tmp_path):
    rows=run_v13_local_sensitivity(tmp_path/'local',target_metrics=('tumor_auc',))
    symbols={r['parameter_symbol'] for r in rows}
    assert 'alpha_ysa_direct_access' not in symbols
    assert 'k_ysa_on' not in symbols
    assert 'k_cross' in symbols and 'R_EphA2_eff' in symbols
    samples=run_v13_monte_carlo(tmp_path/'mc',n_samples=2)
    assert all('k_ysa_on' not in r and 'alpha_ysa_direct_access' not in r for r in samples)
    assert all(r['targeting_model']=='sequential_v1_2' for r in rows+samples)


def test_dual_vs_ysa_grid_has_five_real_dimensions(tmp_path):
    rows=run_v13_feasibility_domain(tmp_path,levels=(.5,1.))
    symbols=['k_ocn_on','k_cross','k_ysa_on','k_int','R_EphA2_eff']
    assert len(rows)==2**5
    assert len({tuple(r[s] for s in symbols) for r in rows})==2**5
    assert all('alpha_ysa_direct_access' not in r for r in rows)
    with (tmp_path/'tables/feasibility_domain_summary.csv').open(encoding='utf-8') as handle:
        summary=next(csv.DictReader(handle))
    assert int(summary['grid_points'])==32
    assert int(summary['grid_dimension'])==5
    assert summary['targeting_model']=='sequential_v1_2'


def test_report_uses_actual_grid_metadata_and_win_fractions(tmp_path,monkeypatch):
    import cabomb_model.v13_reporting as reporting
    params=v13_default_parameters()
    result=v13_simulate(params,design='unmodified')
    first=v13_simulate(params,design='unmodified',recirculate=False)
    detailed=v13_detailed_metrics(result,params,first_pass=first)
    payload={name:[] for name in ['design_comparison','size_comparison','scenario_comparison',
        'local_sensitivity','monte_carlo','mass_balance_audit']}
    payload.update(baseline={'result':result,'params':params},baseline_detailed=[detailed],feasibility_domain=[{}])
    table=tmp_path/'tables/feasibility_domain_summary.csv'
    table.parent.mkdir()
    record={'grid_dimension':5,'levels_per_parameter':2,'grid_points':32,
        'targeting_model':'sequential_v1_2','grid_symbols':'k_ocn_on;k_cross;k_ysa_on;k_int;R_EphA2_eff',
        'dual_wins_tumor_auc_fraction':0.,'dual_wins_internalized_fraction':.25,
        'dual_wins_released_cargo_fraction':.5}
    with table.open('w',encoding='utf-8',newline='') as handle:
        writer=csv.DictWriter(handle,fieldnames=list(record))
        writer.writeheader(); writer.writerow(record)
    monkeypatch.setattr(reporting,'run_v13_all_analyses',lambda output_dir:payload)
    path=reporting.build_v13_markdown_report(tmp_path,tmp_path/'report.md')
    text=path.read_text(encoding='utf-8')
    assert '32 points' in text and '729 points' not in text
    assert '0.0% for tumor AUC' in text
    assert '25.0% for final internalization' in text
    assert 'wins tumor retention (AUC) across the whole' not in text
