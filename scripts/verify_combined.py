#!/usr/bin/env python3
"""One-command combined verification with per-stage logs and honest scope.

quick: complete regression suite, foundation saved-arithmetic checks, exact
independent orbit replay, and packaged degree-eight map arithmetic.
fresh: quick plus fresh characters, rank81, engine bridge, degree-eight tensor
map and degree-twelve heldouts.
full: fresh plus complete bounded degree-ten and degree-twelve map proofs.
"""
import argparse
import hashlib
import json
import os
from pathlib import Path
import shutil
import subprocess
import sys
import time

ROOT=Path(__file__).resolve().parents[1]


def main():
    ap=argparse.ArgumentParser(description=__doc__)
    ap.add_argument('--level',choices=('quick','fresh','full'),default='quick')
    ap.add_argument('--workers',type=int,default=2)
    args=ap.parse_args();os.chdir(ROOT)
    out=Path('reports/combined-run');out.mkdir(parents=True,exist_ok=True)
    py=sys.executable
    stages=[('regression',[py,'-m','pytest','-q']),
            ('foundation-arithmetic',[py,'foundation/science/verify_science.py']),
            ('independent-orbits',[py,'independent/verify.py','--output','reports/independent-algebra.json']),
            ('generalized-completion',[py,'independent/completion.py']),
            ('map8-arithmetic',[py,'-m','chiral4form','verify-map8','verification/map8/map_certificate.json'])]
    if args.level in ('fresh','full'):
        compiler=shutil.which('clang++') or shutil.which('g++')
        if not compiler:raise SystemExit('A C++17 compiler is required for fresh character generation')
        Path('build/combined').mkdir(parents=True,exist_ok=True)
        stages += [('compile-character',[compiler,'-std=c++17','-O3','foundation/science/character/d5_symmetric_character.cpp','-o','build/combined/d5-character']),
                   ('fresh-character',['build/combined/d5-character','22','reports/fresh-d5']),
                   ('python-character',[py,'foundation/science/character/check_character.py','reports/fresh-d5','--full-through','12','--output','reports/fresh-character-python.json']),
                   ('rank81',[py,'foundation/science/rank81/verify.py','foundation/science/rank81/input.json','--output','reports/fresh-rank81.json']),
                   ('tensor-bridge',[py,'scripts/verify_combined_tensor_bridge.py']),
                   ('map8-fresh',[py,'-m','chiral4form','verify-map8','verification/map8/map_certificate.json','--reevaluate-tensors']),
                   ('map12-holdouts',[py,'scripts/check_degree12_physics_holdouts.py'])]
    if args.level=='full':
        stages += [('map10-full',[py,'-c',"from chiral4form.degree10_map_proof import run; run('reports/degree10-map', registry_path='data/combined/registry.json')"]),
                   ('map12-full',[py,'scripts/resume_degree12_proof.py','--workers',str(args.workers)]),
                   ('map12-strict-replay',[py,'scripts/verify_degree12_checkpoint_inputs.py'])]
    env=dict(os.environ,OPENBLAS_NUM_THREADS='1',OMP_NUM_THREADS='1',VECLIB_MAXIMUM_THREADS='1')
    report={'schema':1,'level':args.level,'status':'running','stages':[],
            'scope':{'quick':'regression and saved exact arithmetic; no fresh physics','fresh':'fresh tensor/character checks and physical holdouts; not full degree12 identities','full':'fresh verification plus bounded arithmetic tensor map proofs through degree12'}[args.level]}
    destination=out/f'{args.level}-manifest.json'
    for name,cmd in stages:
        print(f'Running {name}',flush=True);start=time.monotonic();path=out/f'{name}.log'
        with path.open('w') as log:
            result=subprocess.run(cmd,stdout=log,stderr=subprocess.STDOUT,env=env)
        stage={'name':name,'command':cmd,'exit_code':result.returncode,'seconds':time.monotonic()-start,
               'log':str(path),'log_sha256':hashlib.sha256(path.read_bytes()).hexdigest()}
        report['stages'].append(stage)
        report['status']='failed' if result.returncode else 'running'
        destination.write_text(json.dumps(report,indent=2)+'\n')
        if result.returncode:
            print(path.read_text()[-5000:]);raise SystemExit(result.returncode)
    report['status']='passed';destination.write_text(json.dumps(report,indent=2)+'\n')
    print(f'PASS: {args.level}; evidence at {destination}',flush=True)

if __name__=='__main__':main()
