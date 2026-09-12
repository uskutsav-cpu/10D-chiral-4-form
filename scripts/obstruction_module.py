#!/usr/bin/env python3
"""Exact row-span annihilators. This does NOT certify nonlinear reachability."""
import argparse,json
from dataclasses import asdict
from pathlib import Path
from chiral4form.certificates import build_quotient_certificate
from chiral4form.provenance import atomic_json


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('input',type=Path);parser.add_argument('--out',type=Path)
    args=parser.parse_args();data=json.loads(args.input.read_text())
    cert=build_quotient_certificate(data['reachable_rows'],data['ambient_dim'],data['prime'])
    out=asdict(cert)
    if args.out:atomic_json(args.out,out)
    else:print(json.dumps(out,indent=2))

if __name__=='__main__':main()
