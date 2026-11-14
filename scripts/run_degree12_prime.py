#!/usr/bin/env python3
"""Run one declared prime shard, or plan missing compatible proof shards."""
import argparse
import json
from pathlib import Path
from prove_combined_degree12_map import parameters, namespace, worker
from verify_degree12_checkpoint_inputs import validate_checkpoint


def main():
    ap=argparse.ArgumentParser(description=__doc__)
    mode=ap.add_mutually_exclusive_group(required=True)
    mode.add_argument('--plan',action='store_true')
    mode.add_argument('--prime',type=int)
    ap.add_argument('--output',type=Path,default=Path('reports/degree12-map'))
    args=ap.parse_args();reg,rows,primes=parameters();digest,config=namespace(reg,rows,primes)
    if args.plan:
        pending=[]
        for p in primes:
            path=args.output/f'prime{p}.json'
            if path.exists():
                record=json.loads(path.read_text())
                if record.get('status')=='complete':
                    validate_checkpoint(record,p,digest)
                    continue
                if record.get('namespace_sha256')!=digest:
                    raise ValueError(f'incompatible partial shard: {p}')
            pending.append(p)
        print(json.dumps(pending));return
    if args.prime not in primes:
        raise ValueError('prime is outside the declared bounded-proof set')
    args.output.mkdir(parents=True,exist_ok=True)
    path=args.output/'configuration.json'
    expected={'namespace_sha256':digest,'configuration':config}
    if path.exists() and json.loads(path.read_text())!=json.loads(json.dumps(expected)):
        raise ValueError('incompatible proof configuration')
    if not path.exists():
        path.write_text(json.dumps(expected,indent=2,sort_keys=True)+'\n')
    worker(args.output,args.prime,digest)
    validate_checkpoint(json.loads((args.output/f'prime{args.prime}.json').read_text()),args.prime,digest)

if __name__=='__main__':main()
