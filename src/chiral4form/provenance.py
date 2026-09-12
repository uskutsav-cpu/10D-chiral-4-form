"""Atomic artifacts, semantic hashes, and explicit Git provenance."""
from __future__ import annotations
import hashlib
import json
import os
import subprocess
import tempfile
from pathlib import Path


def sha256_file(path: str|Path) -> str:
    h=hashlib.sha256()
    with Path(path).open('rb') as f:
        for chunk in iter(lambda:f.read(1024*1024),b''):
            h.update(chunk)
    return h.hexdigest()


def semantic_hash(value) -> str:
    return hashlib.sha256(json.dumps(value,sort_keys=True,separators=(',',':'),allow_nan=False).encode()).hexdigest()


def atomic_json(path: str|Path,value) -> None:
    path=Path(path)
    path.parent.mkdir(parents=True,exist_ok=True)
    fd,tmp=tempfile.mkstemp(prefix='.'+path.name+'.',dir=path.parent)
    try:
        with os.fdopen(fd,'w') as f:
            json.dump(value,f,indent=2,sort_keys=True,allow_nan=False)
            f.write('\n');f.flush();os.fsync(f.fileno())
        os.replace(tmp,path)
    except BaseException:
        Path(tmp).unlink(missing_ok=True)
        raise


def safe_child(root: str|Path,relative: str) -> Path:
    base=Path(root).resolve()
    if Path(relative).is_absolute():
        raise ValueError('absolute artifact paths are forbidden')
    q=(base/relative).resolve()
    if not q.is_relative_to(base):
        raise ValueError('artifact path escapes root')
    return q


def git(root,*args):
    result=subprocess.run(['git','-C',str(root),*args],text=True,capture_output=True,check=True)
    return result.stdout.strip()


def freeze_checkout(root,relative_paths,expected_commit):
    root=Path(root).resolve()
    head=git(root,'rev-parse','HEAD')
    if head!=expected_commit:
        raise ValueError(f'expected pinned commit {expected_commit}, got {head}')
    if git(root,'status','--porcelain','--untracked-files=no'):
        raise ValueError('source checkout has tracked modifications')
    files=[]
    for rel in relative_paths:
        q=safe_child(root,rel)
        if not q.is_file():
            raise FileNotFoundError(f'required upstream file is missing: {rel}')
        expected=git(root,'rev-parse',f'{head}:{rel}')
        actual=git(root,'hash-object','--',rel)
        if actual!=expected:
            raise ValueError(f'Git object differs from pinned source: {rel}')
        files.append({'path':rel,'git_blob':expected,'sha256':sha256_file(q),'bytes':q.stat().st_size})
    manifest={'schema':1,'commit':head,'files':files,'status':'hash_verified_not_physics_verified'}
    manifest['fingerprint']=semantic_hash(manifest)
    return manifest
