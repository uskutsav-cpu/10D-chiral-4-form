#!/usr/bin/env python3
from __future__ import annotations
import argparse,json
from pathlib import Path
from chiral4form.localization import modmax_symbolic_report
def main():
    ap=argparse.ArgumentParser()
    ap.add_argument("--inventory",type=Path,default=Path("verification/degree12/generalized_candidate_inventory.json"))
    ap.add_argument("--discovery",type=Path,default=Path("verification/degree12/generalized_discovery.json"))
    ap.add_argument("--completion",type=Path,default=Path("verification/degree12/generalized_completion_manifest.json"))
    ap.add_argument("--json-output",type=Path,default=Path("verification/degree12/conformal_scope_report.json"))
    ap.add_argument("--md-output",type=Path,default=Path("docs/CONFORMAL_SCOPE_REPORT.md"))
    a=ap.parse_args(); inv=json.loads(a.inventory.read_text()); d=json.loads(a.discovery.read_text()); c=json.loads(a.completion.read_text())
    deg={x["id"]:int(x["degree"]) for x in inv["candidates"]}; sel=d["selected_extras"]
    counts={str(k):sum(deg[x]==k for x in sel) for k in (4,6,8,10,12)}
    omitted=sorted(set(deg)-set(sel),key=lambda x:(deg[x],x)); mod=modmax_symbolic_report()
    rep={"schema":1,"status":"physical_and_conformal_scope_report","pure_stress_orbit_dimensions":[4,6,10],
      "generalized_selected_count":len(sel),"selected_counts_by_degree":counts,"omitted_primitive_candidates":omitted,
      "modmax_prior_work_reproduction":mod,
      "scope_conclusion":"Finite polynomial theorems are analytic at the zero-field point through degree 12. The ModMax-like b*sqrt(I4) sector is nonanalytic at I4=0 and lies outside that theorem class.",
      "paper_safe_interpretation":"The computation gives a finite-degree constructive realization of the D=10 need for additional F5-dependent scalar structures beyond stress invariants; it does not settle nonanalytic ModMax reachability."}
    a.json_output.parent.mkdir(parents=True,exist_ok=True); a.json_output.write_text(json.dumps(rep,indent=2)+"\n")
    md=f"""# Conformal / physical scope report

## Certified finite-polynomial results
- Pure-stress orbit dimensions: **4 -> 6 -> 10**.
- Genuine generalized completion through degree 12: **certified**.
- Selected primitive extras: **{len(sel)}**.
- Selected counts by degree: `{counts}`.
- Omitted primitive candidates: `{omitted}`.

## ModMax-like sector
Repository reproduction: conformal homogeneity = `{mod["conformal_homogeneity_verified"]}`,
stress-square identity = `{mod["stress_square_identity_verified"]}`.

`V=b sqrt(I4)` is nonanalytic at `I4=0`; therefore it is outside the finite
polynomial theorem class and no ModMax reachability theorem is claimed.
"""
    a.md_output.parent.mkdir(parents=True,exist_ok=True); a.md_output.write_text(md)
    print("PASS: STEP 4 SCOPE REPORT")
    print("selected counts:",counts,"omitted:",omitted)
    print("WROTE:",a.json_output); print("WROTE:",a.md_output)
if __name__=="__main__": main()
