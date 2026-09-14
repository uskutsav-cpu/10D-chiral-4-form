"""Search low-height affine charts for the degree-12 obstruction kernel.

The frozen global theorem certifies a four-dimensional degree-12 fiber and
selects columns [0,1,2,4] as valid vertical coordinates.  Rational RREF height
depends strongly on which four leading coordinates are chosen free, so search
several valid charts before spending compute on new primes.
"""
from __future__ import annotations

import itertools,json
from pathlib import Path

from .degree12_certificate import obstruction_equations
from .degree12_kernel_lift import discover_model_directory,theorem_record
from .finite_field import nullspace,rref
from .normalization_audit import load_records
from .provenance import atomic_json
from .subspace_lift import lift_stable_row_space

N=81
LEADING=72
CERTIFIED_FREE=(0,1,2,4)

def load_modular_kernels():
    theorem=theorem_record()
    model_dir,primes=discover_model_directory()
    records=load_records([model_dir/f"fields_prime{p}.json" for p in primes])
    kernels={}
    labels_ref=None
    for p in primes:
        _,_,labels,equations,_=obstruction_equations(records[p],p)
        if labels_ref is None:
            labels_ref=list(labels)
        elif list(labels)!=labels_ref:
            raise ValueError(f"ansatz labels changed at prime {p}")
        K=nullspace(equations,p,ncols=N)
        if len(K)!=68:
            raise ValueError(f"kernel dimension at prime {p} is {len(K)}, expected 68")
        kernels[p]=K
    return theorem,model_dir,primes,labels_ref,kernels

def chart_spaces(kernels,free):
    free=tuple(sorted(int(x) for x in free))
    if len(free)!=4 or len(set(free))!=4 or any(not 0<=x<LEADING for x in free):
        raise ValueError("chart needs four distinct leading columns")
    piv=[c for c in range(LEADING) if c not in free]
    perm=piv+list(free)+list(range(LEADING,N))
    spaces={}
    for p,K in kernels.items():
        kp=[[row[c] for c in perm] for row in K]
        rr,pv=rref(kp,p,ncols=N)
        if pv[:68] != list(range(68)) or len(pv)!=68:
            return None,None
        rows=[]
        for r in rr[:68]:
            out=[0]*N
            for newc,oldc in enumerate(perm):
                out[oldc]=r[newc]
            rows.append(out)
        spaces[p]=rows
    return spaces,tuple(piv)

def try_chart(kernels,free):
    spaces,piv=chart_spaces(kernels,free)
    if spaces is None:
        return {"free":list(free),"valid":False,"exact":False,"unresolved":10**9}
    lift=lift_stable_row_space(spaces,piv,ncols=N,holdout_primes=(),max_bad_primes=0)
    unresolved=len(lift.get("unresolved",[]))
    return {
        "free":list(free),
        "valid":True,
        "exact":bool(lift.get("qq_verified")),
        "unresolved":unresolved,
        "status":lift.get("status"),
        "method_counts":lift.get("method_counts",{}),
        "lift":lift if lift.get("qq_verified") else None,
    }

def candidate_charts():
    seen=set()
    def emit(x):
        t=tuple(sorted(x))
        if t not in seen:
            seen.add(t);return t
    for c in (CERTIFIED_FREE,(0,1,4,67)):
        t=emit(c)
        if t is not None: yield t
    # Search simple decomposable leading coordinates first.
    for c in itertools.combinations(range(10),4):
        t=emit(c)
        if t is not None: yield t
    # Then one primitive coordinate plus three simple coordinates.
    for z in range(10,72):
        for base in itertools.combinations(range(10),3):
            t=emit((*base,z))
            if t is not None: yield t

def search_charts(
    output="verification/all_orders/degree12_exact/chart_search.json",
    kernel_output="verification/all_orders/degree12_exact/kernel_lift.json",
    max_charts=1200,
    progress=print,
):
    theorem,model_dir,primes,labels,kernels=load_modular_kernels()
    best=None;attempts=[]
    for i,free in enumerate(candidate_charts(),1):
        if i>max_charts: break
        rec=try_chart(kernels,free)
        attempts.append({k:v for k,v in rec.items() if k!="lift"})
        if rec["valid"] and (best is None or rec["unresolved"]<best["unresolved"]):
            best=rec
            progress(f"chart {i}: free={free}, unresolved={rec['unresolved']}, exact={rec['exact']}",flush=True)
        if rec["exact"]:
            result={
                "schema":1,
                "status":"exact_degree12_kernel_lift",
                "model_directory":str(model_dir),
                "primes":list(primes),
                "source_theorem":"verification/degree12/qq/degree12_local_theorem_certificate.json",
                "source_qq_theorem":True,
                "ansatz_labels":labels,
                "ansatz_dimension":81,
                "kernel_dimension":68,
                "leading_rank":68,
                "pivots":rec["lift"]["pivots"],
                "free_degree12_columns":list(free),
                "free_degree12_count":4,
                "lift":rec["lift"],
                "exact":True,
                "chart_search_selected":True,
            }
            kpath=Path(kernel_output);kpath.parent.mkdir(parents=True,exist_ok=True)
            atomic_json(kpath,result)
            summary={
                "schema":1,"status":"exact_chart_found","selected_free":list(free),
                "attempt_count":i,"best_unresolved":0,"attempts":attempts,
            }
            opath=Path(output);opath.parent.mkdir(parents=True,exist_ok=True)
            atomic_json(opath,summary)
            return summary
    summary={
        "schema":1,
        "status":"no_exact_chart_with_current_primes",
        "attempt_count":len(attempts),
        "best_free":None if best is None else best["free"],
        "best_unresolved":None if best is None else best["unresolved"],
        "attempts":attempts,
        "next_action":"add fresh degree-12 model primes adaptively, then rerun chart search",
    }
    opath=Path(output);opath.parent.mkdir(parents=True,exist_ok=True)
    atomic_json(opath,summary)
    return summary
