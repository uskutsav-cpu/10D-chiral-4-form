#!/usr/bin/env python3
from chiral4form.higher_degree_pipeline import run_higher
if __name__=="__main__":
    run_higher("runs/all_orders/degree16_basis/registry_degree16.json",16,
        "runs/all_orders/degree16_map")
    print("PASS: DEGREE-16 MODULAR STRESS MAP FIT + HOLDOUTS")
