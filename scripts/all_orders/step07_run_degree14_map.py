#!/usr/bin/env python3
from chiral4form.higher_degree_pipeline import run_higher
if __name__=="__main__":
    run_higher("runs/all_orders/degree14_basis/registry_degree14.json",14,
        "runs/all_orders/degree14_map")
    print("PASS: DEGREE-14 MODULAR STRESS MAP FIT + HOLDOUTS")
