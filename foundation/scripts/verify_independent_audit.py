#!/usr/bin/env python3
"""Verify saved independent audit evidence without production imports.

This is a low-cost algebra/evidence check. It does not recompute tensor
contractions or rederive the stated Hilbert upper bounds.
"""
import argparse
import hashlib
import importlib.util
import json
from pathlib import Path
import subprocess
import sys


def require(condition,message):
    if not condition:raise ValueError(message)


def load(path):return json.loads(path.read_text())


def main():
    require(__debug__, "verification requires Python without -O or PYTHONOPTIMIZE")
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--source-only",action="store_true",help="verify exact source certificate and canonical maps only")
    parser.add_argument("--repository",type=Path,default=Path(__file__).resolve().parents[1])
    parser.add_argument("--audit-dir",type=Path,help="override results/audit evidence directory")
    args=parser.parse_args()
    repository=args.repository.resolve()
    audit=args.audit_dir.resolve() if args.audit_dir else repository/"results/audit"
    math=audit/"math-review";root=audit/"root-review";source=audit/"source-review"
    certificate=load(source/"exact-source-evaluation-certificate.json")
    for degree in [8,10]:
        canonical=load(repository/f"results/order{degree}_change_of_basis.json")
        require(canonical==certificate["maps"][f"degree{degree}"],f"canonical degree{degree} map differs from exact certificate")
    checked={"canonical_source_maps":True}
    if not args.source_only:
        spec=importlib.util.spec_from_file_location("independent_audit_bareiss",math/"bareiss_minors.py")
        module=importlib.util.module_from_spec(spec);spec.loader.exec_module(module)
        bareiss=module.bareiss
        stored=load(repository/"results/rank81_certificate.json")
        stored_results=load(math/"bareiss_results.json")
        require(hashlib.sha256((repository/"results/rank81_certificate.json").read_bytes()).hexdigest()==stored_results["certificate_sha256"],"frozen certificate hash differs")
        expected={(r["prime"],r["seed"]):r for r in stored_results["results"]}
        require(len(stored["witnesses"])==4,"four stored graph cells required")
        for w in stored["witnesses"]:
            p=w["prime"];columns=w["pivot_columns"]
            require(len(w["jacobian"])==81 and all(len(row)==126 for row in w["jacobian"]),"stored graph matrix shape")
            require(len(columns)==len(set(columns))==81,"stored graph pivot columns")
            determinant,_=bareiss([[row[c] for c in columns] for row in w["jacobian"]])
            prior=expected[p,w["seed"]]
            require(str(determinant)==prior["integer_determinant"],"stored integer graph determinant differs")
            require(determinant%p==w["determinant_mod_p"]==prior["determinant_mod_p"]!=0,"stored graph minor differs")
        checked["stored_graph_minors"]=4
        reverse=load(math/"reverse_gradient_results.json")
        fresh_minors=load(math/"fresh_minor_bareiss_results.json")
        dual=load(math/"independent_fresh_results.json")
        new_witnesses={}
        for index,cell in enumerate(reverse["cells"]):
            p=cell["prime"];w=load(root/f"fresh-witness-p{p}.json")
            new_witnesses[p]=w
            expected_minor=next(r for r in fresh_minors["cells"] if (r["prime"],r["seed"])==(p,cell["seed"]))
            require(len(cell["records"])==83,"83 independent graph rows required")
            for j,record in enumerate(cell["records"][:81]):
                require(record["jacobian_row"]==w["jacobian"][j],"fresh full Jacobian row differs")
                require(record["value"]==w["values"][j]==dual["cells"][index]["values"][j],"fresh graph value differs")
                tangent=sum(a*b for a,b in zip(record["jacobian_row"],cell["direction"]))%p
                require(tangent==dual["cells"][index]["directional_derivatives"][j],"fresh directional derivative differs")
            columns=expected_minor["pivot_columns"]
            determinant,_=bareiss([[r["jacobian_row"][c] for c in columns] for r in cell["records"][:81]])
            require(str(determinant)==expected_minor["integer_determinant"],"fresh integer graph determinant differs")
            require(determinant%p==expected_minor["determinant_mod_prime"]==w["determinant_mod_p"]!=0,"fresh graph minor differs")
        checked["fresh_independent_graph_minors"]=len(reverse["cells"])
        orbit=load(root/"orbit_tangent_certificate.json")
        for cell in orbit["cells"]:
            p=cell["prime"];rows=cell["generator_rows_45x126"];columns=cell["minor_coordinate_indices"]
            require(len(rows)==45 and all(len(r)==126 for r in rows),"orbit matrix shape")
            minor=[[r[c] for c in columns] for r in rows]
            require(minor==cell["minor_45x45"],"orbit minor extraction differs")
            determinant,_=bareiss(minor)
            require(determinant%p==cell["minor_determinant_mod_p"]!=0,"orbit minor differs")
            require(all(sum(a*b for a,b in zip(row,tangent))%p==0 for row in new_witnesses[p]["jacobian"] for tangent in rows),"Jacobian does not annihilate saved orbit directions")
        checked["orbit_minors"]=len(orbit["cells"])
        polynomial=load(math/"degree12_polynomial_gradient_certificate.json")
        for field in polynomial["fields"]:
            p=field["prime"];stacked=[[] for _ in range(72)]
            for pi,point in enumerate(field["points"]):
                if pi==0:
                    records=next(c["records"] for c in reverse["cells"] if c["prime"]==p)
                    records=records[:9]+records[21:]
                else:records=point["records"]
                require(len(records)==71,"degree12 point graph inventory differs")
                low=records[:9];high=records[9:]
                q,a,b=[r["value"] for r in low[:3]]
                gq,ga,gb=[r["jacobian_row"] for r in low[:3]]
                combine=lambda ts:[sum(c*r[j] for c,r in ts)%p for j in range(126)]
                products=[combine([(3*q*q,gq)]),combine([(2*a,ga)]),combine([(b,ga),(a,gb)]),combine([(2*b,gb)])]
                products += [combine([(r["value"],gq),(q,r["jacobian_row"])]) for r in low[3:]]
                rows=products+[r["jacobian_row"] for r in high]
                stacked=[before+after for before,after in zip(stacked,rows)]
            require(stacked==field["stacked_gradient"],"degree12 stacked product/graph rows differ")
            require(len(stacked)==72 and all(len(row)==126*len(field["points"]) for row in stacked),"degree12 stacked matrix shape")
            columns=field["pivot_columns"]
            require(len(columns)==len(set(columns))==72,"degree12 pivot columns")
            determinant,_=bareiss([[row[c] for c in columns] for row in stacked])
            require(str(determinant)==field["integer_minor_determinant"],"degree12 integer determinant differs")
            require(determinant%p==field["determinant_mod_prime"]!=0,"degree12 minor differs")
        checked["degree12_polynomial_minors"]=len(polynomial["fields"])
    command=[sys.executable,str(math/"verify_source_crt.py"),"--math-dir",str(math),"--source-dir",str(source),
             "--repository",str(repository),"--frozen-basis",str(repository/"results/rank81_basis.json"),"--no-write"]
    completed=subprocess.run(command,text=True,capture_output=True)
    if completed.returncode:
        sys.stderr.write(completed.stderr)
        raise RuntimeError("independent exact source CRT verifier rejected evidence")
    source_result=json.loads(completed.stdout)
    checked["source_certificate"]={k:source_result[k] for k in ["point_files_verified","reconstructed_integer_entries","ranks"]}
    print(json.dumps({"status":"PASS","source_only":args.source_only,"checked":checked,
                      "scope":"Saved evidence algebra and hashes; no tensor recomputation. Polynomial identities retain stated invariant-space upper-bound and invariance premises."},indent=2))


if __name__=="__main__":main()
