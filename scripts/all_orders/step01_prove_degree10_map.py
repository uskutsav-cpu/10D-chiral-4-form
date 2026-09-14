#!/usr/bin/env python3
from pathlib import Path
from chiral4form.degree10_map_proof import run

if __name__=="__main__":
    r=run(
        Path("verification/all_orders/degree10_map"),
        source_root=".cache/upstream-degree12-fresh-20260911-223942",
    )
    print("PASS: UNCONDITIONAL DEGREE-10 PHYSICS/BASIS MAP")
    print("registry:",r["registry_provenance"])
    print("certificate:",r["certificate_hash"])
