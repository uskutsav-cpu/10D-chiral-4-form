#!/usr/bin/env python3
from __future__ import annotations
import argparse,json
from pathlib import Path
from chiral4form.obstruction_ideal import degree14_probe_readiness
from chiral4form.registry import Registry

def main():
    p=argparse.ArgumentParser()
    p.add_argument("--registry",type=Path,default=Path("data/fixtures/low_degree_registry.json"))
    p.add_argument("--output",type=Path,default=Path("verification/all_orders/degree14_probe.json"))
    a=p.parse_args(); r=degree14_probe_readiness(Registry.load(a.registry))
    a.output.parent.mkdir(parents=True,exist_ok=True);a.output.write_text(json.dumps(r,indent=2)+"\n")
    print("DEGREE-14 FALSIFICATION PROBE");print("status:",r["status"]);print("WROTE:",a.output)
if __name__=="__main__":main()
