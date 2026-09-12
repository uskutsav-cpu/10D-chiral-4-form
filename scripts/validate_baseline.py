#!/usr/bin/env python3
from pathlib import Path

from chiral4form import load_baseline

ROOT = Path(__file__).resolve().parents[1]
BASELINE = ROOT / "data" / "baseline" / "paper2_baseline.json"

baseline = load_baseline(BASELINE)
print(f"status: {baseline.status}")
print(f"source: {baseline.source_repository}@{baseline.source_ref}")
print("degree  full  static  reachable  quotient  forcing")
for r in baseline.records:
    print(
        f"{r.degree:>6} {r.full_dim:>5} {r.static_stress_dim:>7} "
        f"{r.dynamic_reachable_dim:>10} {r.quotient_dim:>9} "
        f"{str(r.new_forcing_dim):>8}"
    )
print("PASS: baseline is internally consistent (not an independent physics recomputation)")
