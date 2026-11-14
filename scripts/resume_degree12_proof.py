#!/usr/bin/env python3
"""Resume hash-bound degree-12 prime jobs and report failures immediately."""
import argparse
import json
import os
from pathlib import Path
import subprocess
import sys
import time
from prove_combined_degree12_map import parameters, namespace, verify


def main():
    ap=argparse.ArgumentParser(description=__doc__)
    ap.add_argument('--output',type=Path,default=Path('reports/degree12-map'))
    ap.add_argument('--workers',type=int,default=2)
    args=ap.parse_args()
    if not 1<=args.workers<=2:raise ValueError('one or two workers required')
    reg,rows,primes=parameters();digest,config=namespace(reg,rows,primes)
    args.output.mkdir(parents=True,exist_ok=True)
    config_path=args.output/'configuration.json'
    if not config_path.exists():
        config_path.write_text(json.dumps({"namespace_sha256":digest,"configuration":config},indent=2)+"\n")
    configuration=json.loads(config_path.read_text())
    if configuration['namespace_sha256']!=digest:raise ValueError('checkpoint namespace changed')
    pending=[]
    for p in primes:
        path=args.output/f'prime{p}.json'
        if path.exists():
            record=json.loads(path.read_text())
            if record['namespace_sha256']!=digest:raise ValueError(f'checkpoint hash changed: {p}')
            if len(record['samples'])==72:continue
        pending.append(p)
    running={};env=dict(os.environ,OPENBLAS_NUM_THREADS='1',OMP_NUM_THREADS='1',VECLIB_MAXIMUM_THREADS='1')
    code='from scripts.prove_combined_degree12_map import worker; import sys; worker(sys.argv[1],int(sys.argv[2]),sys.argv[3])'
    status_path=args.output/'supervisor.json'
    try:
        while pending or running:
            while pending and len(running)<args.workers:
                p=pending.pop(0);log=(args.output/f'worker{p}.log').open('a')
                process=subprocess.Popen([sys.executable,'-u','-c',code,str(args.output),str(p),digest],stdout=log,stderr=subprocess.STDOUT,env=env)
                running[p]=(process,log);print(f'started prime={p} pid={process.pid}',flush=True)
            status_path.write_text(json.dumps({'status':'running','active':{str(p):proc.pid for p,(proc,log) in running.items()},'pending':pending},indent=2)+'\n')
            for p,(process,log) in list(running.items()):
                code_out=process.poll()
                if code_out is None:continue
                log.close();del running[p]
                if code_out:
                    print((args.output/f'worker{p}.log').read_text()[-5000:],flush=True)
                    raise RuntimeError(f'prime {p} failed with exit {code_out}')
                print(f'completed prime={p}',flush=True)
            time.sleep(1)
        verify(args.output)
        status_path.write_text(json.dumps({'status':'complete','all_prime_jobs_complete':True},indent=2)+'\n')
    except BaseException as exc:
        for process,log in running.values():
            process.terminate();process.wait();log.close()
        status_path.write_text(json.dumps({'status':'stopped','reason':str(exc)},indent=2)+'\n')
        raise

if __name__=='__main__':main()
