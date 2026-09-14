#!/usr/bin/env python3
"""Run the maximal program through theorem gate #11.

Heavy stages are resumable. A mathematically blocked theorem writes a blocker
certificate instead of being promoted.
"""
from __future__ import annotations
import json,subprocess,sys
from pathlib import Path

def run(name,required=True):
    here=Path(__file__).resolve().parent
    print(f"\n===== {name} =====",flush=True)
    cp=subprocess.run([sys.executable,str(here/name)])
    if cp.returncode and required: raise SystemExit(cp.returncode)
    return cp.returncode

def main():
    # 1-5: low-order all-orders promotion.
    run("step01_prove_degree10_map.py")
    run("step02_promote_degree10.py")
    run("step03_build_degree12_ideal.py")
    run("step04_prove_degree12_tangency.py")
    run("step05_promote_degree12.py")

    # 6-8: first genuinely new finite jet.
    run("step06_build_degree14_basis.py")
    run("step07_run_degree14_map.py")
    run("step08_compute_degree14_quotient.py")

    q14p=Path("verification/all_orders/degree14_quotient.json")
    q14=json.loads(q14p.read_text()) if q14p.exists() else {}
    if q14.get("stable_new_dimension")==0:
        run("step09_build_degree16_basis.py")
        run("step09_run_degree16_map.py")
        run("step09_compute_degree16_quotient.py")
    else:
        print("\nDEGREE 16 SKIPPED: degree-14 quotient is nonzero or unresolved.",flush=True)
        Path("verification/all_orders/degree16_probe.json").write_text(
            json.dumps({"schema":1,"status":"degree16_probe_not_triggered",
                        "reason":"degree14 quotient nonzero or unresolved"},indent=2)+"\n")

    # 10-11 are theorem gates, never promoted from finite evidence alone.
    run("step10_stabilization_theorem.py")
    run("step11_formal_orbit_equality.py")

    print("\n==============================================")
    print("MAXIMAL ALL-ORDERS PIPELINE REACHED GATE #11")
    print("==============================================")
    print("Read verification/all_orders/*.json for exact theorem/blocker status.")

if __name__=="__main__":main()
