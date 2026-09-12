import json,subprocess
from pathlib import Path
import pytest
from chiral4form.provenance import atomic_json,semantic_hash,safe_child,freeze_checkout
from chiral4form.checkpoints import CheckpointStore

@pytest.mark.parametrize("name", ['../escape','/tmp/escape','a/../../escape'])
def test_path_traversal_rejected(tmp_path,name):
    with pytest.raises(ValueError):safe_child(tmp_path,name)

def test_symlink_escape_rejected(tmp_path):
    (tmp_path/'out').symlink_to('/tmp',target_is_directory=True)
    with pytest.raises(ValueError):safe_child(tmp_path,'out/escaped.json')

def test_atomic_semantic_json(tmp_path):
    p=tmp_path/'dir'/'file.json'
    atomic_json(p,{'b':2,'a':1})
    assert json.loads(p.read_text())=={'a':1,'b':2}
    assert semantic_hash({'a':1,'b':2})==semantic_hash({'b':2,'a':1})

def test_cache_resume_and_corruption(tmp_path):
    c=CheckpointStore(tmp_path,{'version':1});calls=[]
    fn=lambda:(calls.append(1) or {'answer':42})
    assert c.run('step',{'x':1},fn)==({'answer':42},False)
    assert c.run('step',{'x':1},fn)==({'answer':42},True)
    assert calls==[1]
    path=tmp_path/(c.key('step',{'x':1})+'.json')
    data=json.loads(path.read_text());data['payload']['answer']=43;path.write_text(json.dumps(data))
    with pytest.raises(ValueError):c.run('step',{'x':1},fn)

def test_changed_inputs_invalidate_cache(tmp_path):
    c=CheckpointStore(tmp_path,{'version':1})
    assert c.key('x',{'seed':1})!=c.key('x',{'seed':2})
    other=CheckpointStore(tmp_path,{'version':2})
    assert c.key('x',{})!=other.key('x',{})

def test_failure_and_lock_cleanup(tmp_path):
    c=CheckpointStore(tmp_path,{})
    def fail():raise ValueError('intentional negative control')
    with pytest.raises(ValueError):c.run('x',{},fail)
    assert list(tmp_path.glob('*.failed.json'))
    assert not list(tmp_path.glob('*.lock'))
    assert c.run('x',{},lambda:5)==(5,False)
    assert not list(tmp_path.glob('*.failed.json'))

def test_lock_not_silently_stolen(tmp_path):
    c=CheckpointStore(tmp_path,{})
    with c.lock('key'):
        with pytest.raises(RuntimeError):
            with c.lock('key'):pass

def test_freeze_git_checkout_requires_exact_clean_state(tmp_path):
    def git(*args):
        return subprocess.check_output(['git','-C',str(tmp_path),*args],text=True).strip()
    git('init','-q');git('config','user.name','Test');git('config','user.email','test@example.invalid')
    (tmp_path/'data.json').write_text('{"value":1}\n')
    git('add','data.json');git('commit','-qm','fixture')
    sha=git('rev-parse','HEAD')
    result=freeze_checkout(tmp_path,['data.json'],sha)
    assert sha in str(result)
    with pytest.raises(FileNotFoundError):freeze_checkout(tmp_path,['absent.json'],sha)
    with pytest.raises((ValueError,RuntimeError)):freeze_checkout(tmp_path,['data.json'],'0'*40)
    (tmp_path/'data.json').write_text('{"value":2}\n')
    with pytest.raises((ValueError,RuntimeError)):freeze_checkout(tmp_path,['data.json'],sha)
