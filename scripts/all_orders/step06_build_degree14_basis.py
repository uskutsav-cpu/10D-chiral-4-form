#!/usr/bin/env python3
from chiral4form.higher_degree_basis import build_higher_basis
if __name__=="__main__":
    build_higher_basis("data/fixtures/low_degree_registry.json",14,
        "runs/all_orders/degree14_basis")
    print("PASS: DEGREE-14 BASIS RANK 247")
