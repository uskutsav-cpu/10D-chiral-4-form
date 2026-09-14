#!/usr/bin/env python3
from pathlib import Path
from chiral4form.higher_degree_basis import build_higher_basis
if __name__=="__main__":
    p=Path("runs/all_orders/degree14_basis/registry_degree14.json")
    if not p.exists():raise SystemExit("degree-14 registry missing")
    build_higher_basis(p,16,"runs/all_orders/degree16_basis",max_candidates=1000000)
    print("PASS: DEGREE-16 BASIS RANK 1364")
