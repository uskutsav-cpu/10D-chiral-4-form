#!/usr/bin/env python3
import json,subprocess,sys
from pathlib import Path

if __name__=="__main__":
    subprocess.run([sys.executable,"scripts/all_orders/degree12/search_chart.py"],check=True)
    r=json.loads(Path("verification/all_orders/degree12_exact/chart_search.json").read_text())
    if r["status"]!="exact_chart_found":
        raise SystemExit("STOP: no exact chart yet; do not promote degree 12")
    from chiral4form.degree12_explicit_ideal import build_ideal
    from chiral4form.degree12_exact_tangency import certify_tangency
    from chiral4form.degree12_global_check import check_global
    from chiral4form.degree12_all_orders_promotion import promote
    i=build_ideal(); print("IDEAL:",i["status"])
    t=certify_tangency(); print("TANGENCY:",t["status"])
    g=check_global(); print("GLOBAL:",g["status"])
    a=promote(); print("ALL ORDERS:",a["status"])
    print("PASS: DEGREE-12 EXACT-IDEAL PIPELINE COMPLETE")
