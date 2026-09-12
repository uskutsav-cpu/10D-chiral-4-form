import json
from pathlib import Path
from copy import deepcopy
import pytest
from chiral4form.pipeline import run

ROOT=Path(__file__).resolve().parents[1]

@pytest.mark.integration
def test_actual_degree4_pipeline_and_checkpoint_replay(tmp_path):
    c=json.loads((ROOT/'configs/smoke.json').read_text())
    c.update(max_degree=4,primes=[30011],fit_margin=1,holdouts=2)
    logs=[]
    report=run(c,ROOT,tmp_path,progress=lambda msg,**kwargs:logs.append(msg))
    assert report['state']=='completed_declared_computations'
    assert report['prime_reports'][0]['basis_ranks']=={'4':1}
    assert report['prime_reports'][0]['finite_depth_lie_rank']['rank']==1
    logs.clear()
    repeated=run(c,ROOT,tmp_path,progress=lambda msg,**kwargs:logs.append(msg))
    assert repeated['run_fingerprint']==report['run_fingerprint']
    assert any('cached' in x for x in logs)
    assert all('computed' not in x for x in logs if 'sample' in x)


def test_resource_preflight_refuses_before_tensor_work(tmp_path):
    c=json.loads((ROOT/'configs/smoke.json').read_text());c['max_bytes']=1
    with pytest.raises(MemoryError):run(c,ROOT,tmp_path,progress=lambda *a,**k:None)
    status=json.loads((tmp_path/'status.json').read_text())
    assert status['state']=='preflight_refused'
    assert not (tmp_path/'checkpoints').exists()
