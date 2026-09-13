#!/usr/bin/env python3
from __future__ import annotations
import argparse, json
from pathlib import Path

def main():
    p=argparse.ArgumentParser()
    p.add_argument("--discovery",type=Path,default=Path("verification/degree12/generalized_discovery.json"))
    p.add_argument("--output",type=Path,default=Path("configs/degree12_generalized_selected_verify.json"))
    p.add_argument("--primes",nargs="+",type=int,default=[33037,33049,33053])
    p.add_argument("--seed-start",type=int,default=2026092501)
    args=p.parse_args()
    d=json.loads(args.discovery.read_text())
    extras=list(d["selected_extras"])
    if len(extras)!=int(d["selected_count"]): raise SystemExit("selected_count mismatch")
    cfg={"schema":1,"max_degree":12,"primes":list(args.primes),"seed_start":args.seed_start,
         "fit_margin":4,"holdouts":4,"registry":"data/fixtures/low_degree_registry.json",
         "extras":extras,"max_bytes":5583457484,"max_multiply_adds":200000000000,
         "exact_blas":True,"lie_depth":3,"flow_analysis":False}
    args.output.parent.mkdir(parents=True,exist_ok=True)
    args.output.write_text(json.dumps(cfg,indent=2)+"\n")
    print("PASS: STEP 3C FRESH CONFIG")
    print("selected extras:",len(extras))
    print("fresh primes:",args.primes)
    print("WROTE:",args.output)
if __name__=="__main__": main()
