"""Lift the degree-12 obstruction kernel directly from frozen modular models."""
from __future__ import annotations
import json,re
from pathlib import Path
from .degree12_certificate import obstruction_equations
from .finite_field import nullspace
from .normalization_audit import load_records
from .provenance import atomic_json
from .subspace_lift import canonical_row_spaces,lift_stable_row_space

THEOREM_PATH=Path("verification/degree12/qq/degree12_local_theorem_certificate.json")

def theorem_record():
    rec=json.loads(THEOREM_PATH.read_text())
    for k,v in {
        "qq_theorem":True,
        "obstruction_nullity":68,
        "degree12_leading_obstruction_rank":68,
        "reduced_ambient_dimension":78,
        "local_orbit_dimension":10,
    }.items():
        if rec.get(k)!=v:
            raise ValueError(f"frozen theorem gate failed: {k}={rec.get(k)!r}, expected {v!r}")
    eq=rec.get("equation_rank",{})
    if eq.get("lower_bound")!=13 or eq.get("upper_bound")!=13:
        raise ValueError("frozen degree-12 equation rank is not exactly 13")
    return rec

def _prime_from_name(path):
    m=re.search(r"fields_prime(\d+)\.json$",path.name)
    return int(m.group(1)) if m else None

def discover_model_directory():
    theorem=theorem_record()
    required=tuple(int(x) for x in theorem["primes"])
    need=set(required);groups={}
    for path in Path("runs").rglob("fields_prime*.json"):
        p=_prime_from_name(path)
        if p is not None:
            groups.setdefault(path.parent,set()).add(p)
    exact=[p for p,ps in groups.items() if need<=ps]
    if not exact:
        partial=sorted(
            ((len(need&ps),str(p),sorted(need-ps)) for p,ps in groups.items()),
            reverse=True
        )[:10]
        raise FileNotFoundError("No run directory has all theorem primes. Best matches: "+json.dumps(partial))
    return sorted(exact,key=lambda p:(len(str(p)),str(p)))[0],required

def lift_kernel(output="verification/all_orders/degree12_exact/kernel_lift.json"):
    theorem=theorem_record()
    model_dir,primes=discover_model_directory()
    paths=[model_dir/f"fields_prime{p}.json" for p in primes]
    records=load_records(paths)
    kernels={};labels_ref=None
    for p in primes:
        _,_,labels,equations,_=obstruction_equations(records[p],p)
        if labels_ref is None:
            labels_ref=list(labels)
        elif list(labels)!=labels_ref:
            raise ValueError(f"ansatz labels changed at prime {p}")
        K=nullspace(equations,p,ncols=81)
        if len(K)!=68:
            raise ValueError(f"expected kernel dimension 68 at prime {p}, got {len(K)}")
        kernels[p]=K
    spaces,pivots=canonical_row_spaces(kernels,ncols=81,expected_rank=68)
    if any(c>=72 for c in pivots):
        raise AssertionError("leading kernel rank is not 68")
    free=[c for c in range(72) if c not in pivots]
    if len(free)!=4:
        raise AssertionError("expected exactly 4 free degree-12 coordinates")
    lift=lift_stable_row_space(spaces,pivots,ncols=81,holdout_primes=(),max_bad_primes=0)
    result={
        "schema":1,
        "status":"exact_degree12_kernel_lift" if lift.get("qq_verified") else "blocked_degree12_kernel_lift",
        "model_directory":str(model_dir),
        "primes":list(primes),
        "source_theorem":str(THEOREM_PATH),
        "source_qq_theorem":bool(theorem["qq_theorem"]),
        "ansatz_labels":labels_ref,
        "ansatz_dimension":81,
        "kernel_dimension":68,
        "leading_rank":68,
        "pivots":list(pivots),
        "free_degree12_columns":free,
        "free_degree12_count":len(free),
        "lift":lift,
        "exact":bool(lift.get("qq_verified")),
    }
    path=Path(output);path.parent.mkdir(parents=True,exist_ok=True)
    atomic_json(path,result)
    return result
