#!/usr/bin/env python3
from __future__ import annotations
import argparse,json
from pathlib import Path
from chiral4form.formal_filtration import all_orders_filtration_theorem

def main():
    p=argparse.ArgumentParser()
    p.add_argument("--output",type=Path,default=Path("verification/all_orders/filtration_theorem.json"))
    a=p.parse_args(); r=all_orders_filtration_theorem()
    a.output.parent.mkdir(parents=True,exist_ok=True)
    a.output.write_text(json.dumps(r,indent=2)+"\n")
    print("PASS: FORMAL FILTRATION THEOREM")
    print(r["theorem"]);print("WROTE:",a.output)
if __name__=="__main__":main()
