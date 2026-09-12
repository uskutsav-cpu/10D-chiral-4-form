from copy import deepcopy
from pathlib import Path
import json
import pytest
from chiral4form.fitting import fit_degree,vector_fields_from_fits,fields_to_json
from chiral4form.pipeline import validate_config,execution_plan
from chiral4form.polynomial import Poly
from chiral4form.completion import search
from chiral4form.lifting import lift_models
from chiral4form.cli import main

ROOT=Path(__file__).resolve().parents[1]

def samples():
    return [{'seed':i,'prime':101,'basis_values':{'4':[x]},'targets':[
        {'generator':'tr2','degree':4,'monomial':[],'value':2*x%101}]} for i,x in enumerate([2,3,7,11])]

def test_holdout_fit_retains_coupling_monomials():
    s=samples();fit=fit_degree(s,['I4_1'],4,101,2)
    assert fit['coordinates'][0]['row']==[2]
    assert fit['holdout_seed_indices']==[2,3]
    ids,fields=vector_fields_from_fits([fit])
    assert ids==('I4_1',) and fields['tr2'][0]==2

@pytest.mark.parametrize("kind", ['seed','holdout','rank','noprime'])
def test_malformed_fitting_inputs(kind):
    s=samples()
    if kind=='seed':s[-1]['seed']=s[0]['seed']
    elif kind=='holdout':s[-1]['targets'][0]['value']+=1
    elif kind=='rank':
        for x in s:x['basis_values']['4']=[0]
    else:s[-1]['prime']=103
    with pytest.raises(ValueError):fit_degree(s,['I4_1'],4,101,2)

def test_missing_generalized_rows_not_invented():
    f=Poly.constant(1,1,101)
    record=fields_to_json(['I4_1'],{'tr2':(f,)},'test')
    record['generator_catalogue']=[{'id':'tr2','factors':['tr2']}]
    with pytest.raises(ValueError):search(record,['I6_2'],1)

def test_bounded_lift_still_not_symbolic_identity():
    models=[]
    for p in (30011,30013,30029):
        models.append(fields_to_json(['x'],{'X':(Poly.constant(1,'32/3',p),)},'test'))
    r=lift_models(models,500)
    assert 'candidate' in str(r).lower()
    assert '32/3' in str(r)

@pytest.mark.parametrize("change", [{'max_degree':14},{'primes':[101,101]},{'holdouts':0},{'lie_depth':0},{'lie_depth':9},{'typo':5}])
def test_bad_run_config(change):
    c=json.loads((ROOT/'configs/smoke.json').read_text());c.update(change)
    with pytest.raises((ValueError,TypeError)):validate_config(c)

def test_fixture_does_not_pretend_degree12_is_available():
    c=json.loads((ROOT/'configs/degree12.json').read_text())
    with pytest.raises(ValueError):execution_plan(c,ROOT)

def test_cli_plan_and_error_codes(capsys):
    assert main(['--root',str(ROOT),'plan'])==0
    assert 'execution_plan_only' in capsys.readouterr().out
    assert main(['--root',str(ROOT),'run','--config','missing.json'])==2

def test_cli_symbolic_outputs(tmp_path):
    for command in ('sextic','degree8','modmax'):
        out=tmp_path/(command+'.json')
        assert main([command,'--output',str(out)])==0
        assert out.exists()

@pytest.mark.parametrize('change',[{'fit_margin':1.5},{'holdouts':True},{'exact_blas':'false'},{'primes':[3]},{'max_bytes':12.5}])
def test_config_types_are_strict(change):
    c=json.loads((ROOT/'configs/smoke.json').read_text());c.update(change)
    with pytest.raises((TypeError,ValueError)):validate_config(c)

def test_worker_index_is_not_guessed(tmp_path):
    from chiral4form.worker import run_prime
    c=json.loads((ROOT/'configs/smoke.json').read_text())
    with pytest.raises(ValueError):run_prime(c,3,ROOT,tmp_path)
