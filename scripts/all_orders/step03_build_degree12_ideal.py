#!/usr/bin/env python3
from chiral4form.degree12_ideal_builder import build_degree12_ideal
if __name__=="__main__":
    r=build_degree12_ideal()
    print("DEGREE-12 IDEAL STATUS:",r["status"])
    if not r["status"].startswith("conditional_exact"):
        print("This is a real blocker; add independent reduced degree-12 prime models and rerun.")
