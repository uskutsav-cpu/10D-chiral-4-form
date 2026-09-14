#!/usr/bin/env python3
from pathlib import Path
from chiral4form.degree10_map_proof import run
if __name__=="__main__":
    r=run(Path("verification/all_orders/degree10_map"))
    print("PASS: UNCONDITIONAL DEGREE-10 PHYSICS/BASIS MAP")
    print("certificate:",r["certificate_hash"])
