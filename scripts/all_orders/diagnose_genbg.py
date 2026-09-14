#!/usr/bin/env python3
import json
from pathlib import Path

from chiral4form.nauty_multigraphs import (
    enumerator_smoke_test,
    genbg_command,
)
from chiral4form.provenance import atomic_json

if __name__=="__main__":
    report={}
    for degree in (4, 14, 16):
        print()
        print("="*60)
        print("ENUMERATOR SMOKE TEST", degree)
        print("="*60)
        if degree == 4:
            # Degree 4 doesn't satisfy the explicit circulant smoke constructor,
            # so only test streaming.
            from itertools import islice
            from chiral4form.nauty_multigraphs import stream_multigraphs, graph_label
            gs=list(islice(stream_multigraphs(4),3))
            r={
                "command":genbg_command(4),
                "streamed_count":len(gs),
                "first_labels":[graph_label(g) for g in gs],
                "passed":len(gs)==3,
            }
        else:
            r=enumerator_smoke_test(degree,1)
        report[str(degree)]=r
        print(json.dumps(r,indent=2))
        if not r["passed"]:
            raise SystemExit(f"degree-{degree} enumerator smoke test failed")

    out=Path("verification/all_orders/genbg_diagnostic.json")
    out.parent.mkdir(parents=True,exist_ok=True)
    atomic_json(out,{
        "schema":1,
        "status":"genbg_large_graph_enumerator_verified",
        "degrees":report,
    })
    print()
    print("PASS: genbg/genbgL diagnostics")
    print("WROTE:",out)
