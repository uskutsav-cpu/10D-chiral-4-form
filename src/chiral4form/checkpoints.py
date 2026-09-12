"""Content-addressed resumable stages. Corrupt/stale success records never pass."""
from __future__ import annotations
from contextlib import contextmanager
import json
import os
import time
from pathlib import Path
from .provenance import atomic_json,semantic_hash

class CheckpointStore:
    def __init__(self,root,namespace):
        self.root=Path(root);self.root.mkdir(parents=True,exist_ok=True)
        self.namespace=namespace

    def key(self,stage,inputs):
        return semantic_hash({'namespace':self.namespace,'stage':stage,'inputs':inputs})

    def load(self,stage,inputs):
        key=self.key(stage,inputs);path=self.root/(key+'.json')
        if not path.exists():
            return None
        data=json.loads(path.read_text())
        if data.get('key')!=key or data.get('payload_sha256')!=semantic_hash(data.get('payload')):
            raise ValueError(f'corrupt checkpoint: {path}')
        return data['payload']

    @contextmanager
    def lock(self,key):
        path=self.root/(key+'.lock')
        try:
            fd=os.open(path,os.O_CREAT|os.O_EXCL|os.O_WRONLY,0o600)
        except FileExistsError as exc:
            raise RuntimeError(f'active/stale lock {path}; inspect PID before removing it') from exc
        try:
            with os.fdopen(fd,'w') as f:
                json.dump({'pid':os.getpid(),'created_unix':time.time()},f)
            yield
        finally:
            path.unlink(missing_ok=True)

    def run(self,stage,inputs,fn):
        key=self.key(stage,inputs)
        with self.lock(key):
            cached=self.load(stage,inputs)
            if cached is not None:
                return cached,True
            try:
                value=fn()
                atomic_json(self.root/(key+'.json'),{'key':key,'stage':stage,
                    'payload_sha256':semantic_hash(value),'payload':value})
                (self.root/(key+'.failed.json')).unlink(missing_ok=True)
                return value,False
            except BaseException as exc:
                atomic_json(self.root/(key+'.failed.json'),{'key':key,'stage':stage,
                    'exception':type(exc).__name__,'message':str(exc),'status':'failed'})
                raise
