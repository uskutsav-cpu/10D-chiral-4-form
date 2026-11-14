#!/usr/bin/env python3
"""Strictly validate raw degree-twelve checkpoint inputs before proof replay."""
import argparse
import hashlib
import json
from pathlib import Path
from prove_combined_degree12_map import SEEDS, parameters, namespace, verify


def validate_checkpoint(record, prime, digest):
    if record.get('schema') != 1 or record.get('prime') != prime or type(record.get('prime')) is not int:
        raise ValueError('checkpoint prime/schema mismatch')
    if record.get('namespace_sha256') != digest or record.get('seeds') != list(SEEDS):
        raise ValueError('checkpoint namespace/seed mismatch')
    samples=record.get('samples')
    if not isinstance(samples,list) or len(samples)!=72 or record.get('status')!='complete':
        raise ValueError('72 completed physical samples required')
    for i,sample in enumerate(samples):
        if sample.get('seed')!=SEEDS[i]:raise ValueError('sample seed order mismatch')
        for name,size in [('basis',72),('targets',36)]:
            values=sample.get(name)
            if not isinstance(values,list) or len(values)!=size:
                raise ValueError(f'{name} vector shape mismatch')
            if any(type(value) is not int or not 0<=value<prime for value in values):
                raise ValueError(f'{name} entries are not canonical finite-field integers')


def main():
    ap=argparse.ArgumentParser(description=__doc__);ap.add_argument('--output',type=Path,default=Path('reports/degree12-map'));args=ap.parse_args()
    reg,rows,primes=parameters();digest,_=namespace(reg,rows,primes)
    for p in primes:
        validate_checkpoint(json.loads((args.output/f'prime{p}.json').read_text()),p,digest)
    certificate=verify(args.output,write=False)
    report={'schema':1,'passed':True,'raw_checkpoint_schema_verified':True,'all_bounded_arithmetic_replayed':True,
            'map_namespace':digest,'source_map_status':certificate['status'],
            'validator_sha256':hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
            'scope':'strict canonical prime/point/vector inputs plus complete bounded-CRT identity replay'}
    (args.output/'strict-input-verification.json').write_text(json.dumps(report,indent=2)+'\n')
    print(json.dumps(report))

if __name__=='__main__':main()
