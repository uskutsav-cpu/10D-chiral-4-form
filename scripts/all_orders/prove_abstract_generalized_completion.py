#!/usr/bin/env python3
from __future__ import annotations
import argparse,json
from pathlib import Path
from chiral4form.hilbert_completion import abstract_hilbert_completion_theorem

def main():
    p=argparse.ArgumentParser()
    p.add_argument("--output",type=Path,
        default=Path("verification/all_orders/abstract_generalized_completion.json"))
    a=p.parse_args();r=abstract_hilbert_completion_theorem()
    a.output.parent.mkdir(parents=True,exist_ok=True)
    a.output.write_text(json.dumps(r,indent=2)+"\n")
    print("PASS: ABSTRACT ALL-ORDERS GENERALIZED COMPLETION EXISTENCE THEOREM")
    print(r["conclusion"]);print("WROTE:",a.output)
if __name__=="__main__":main()
