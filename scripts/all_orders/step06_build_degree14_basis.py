#!/usr/bin/env python3
import subprocess
import sys

from chiral4form.higher_degree_basis import build_higher_basis
from chiral4form.registry_resolution import freeze_resolved_registry

if __name__=="__main__":
    subprocess.run(
        [sys.executable,"scripts/all_orders/diagnose_genbg.py"],
        check=True,
    )

    base=freeze_resolved_registry(
        12,
        "runs/all_orders/base_registry_degree12.json",
        source_root=".cache/upstream-degree12-fresh-20260911-223942",
    )

    build_higher_basis(
        base,
        14,
        "runs/all_orders/degree14_basis",
    )
    print("PASS: DEGREE-14 BASIS RANK 247")
