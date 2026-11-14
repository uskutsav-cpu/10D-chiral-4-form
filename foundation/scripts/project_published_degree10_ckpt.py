#!/usr/bin/env python3
"""Project the equation-(4.24) candidates into Q10, with durable checkpoints.

Why this exists alongside `project_published_degree10.py`
---------------------------------------------------------
The dominant cost of a degree-10 projection is rebuilding the 14 atlas columns
at 22 samples for each of 6 primes. That work is *identical* every time a new
published candidate is implemented, and re-paying it turned a one-formula
addition into a full multi-minute rerun. Here every evaluation is an immutable
checkpoint unit, so adding a formula only evaluates that formula.

Checkpoint location
-------------------
Active checkpoints must NOT live on an eventually-consistent filesystem. The
canonical tree is under ~/Documents, which is iCloud-synced, so the default
root is a local temp path and never inside the repository. Override with
SDINV_CKPT_ROOT.

Unit layout
-----------
Atlas columns occupy column ids 0..len(names)-1. Published candidates occupy
ids 1000+ so the two namespaces can never collide, and so a change in atlas
size cannot silently realias a stored formula value.
"""
import hashlib
import json
import os
import subprocess
import sys
import time
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))
sys.path.insert(0, str(ROOT / "scripts"))
sys.path.insert(0, str(ROOT / "tests"))

from sdinv.forms import selfdual_projector, to_dense, random_form
from sdinv.exactmap import rank_mod
from sdinv.projection_checkpoint import (
    ProjectionCheckpoint, block_fingerprint, evaluator_fingerprint, peak_rss_mb)
from sdinv.published_degree10_invariants import (
    AMBIGUITY_VARIANTS, NOT_IMPLEMENTED, PUBLISHED_DEGREE10)
from solve_intrinsic_quotients import rref, project
from stress_flow_closure import closure_span
from test_M_only_quotients import (
    registry_items, evaluate_atlas_element, solve_exact)

CERT = ROOT / "results" / "stress_flow" / "certificates"
FIT_PRIMES = (32749, 32719, 32693, 32771)
HOLDOUT_PRIMES = (32713, 32717)
PRIMES = FIT_PRIMES + HOLDOUT_PRIMES
SEED8 = {4: ["I4_1"], 6: ["I6_2"], 8: ["I8_3", "I8_4", "I8_5", "I8_6"]}
FORMULA_ID_BASE = 1000
EVALUATOR_VERSION = "p10-2026-07-30-boostfix"
MODULAR_BACKEND = "mod_einsum-int64-pairwise-v1"

DEFAULT_CKPT = Path(
    os.environ.get("SDINV_CKPT_ROOT")
    or Path(os.environ.get("TMPDIR", "/tmp")) / "sdinv_ckpt") / "p10"


def sample(prime, seed):
    projector = selfdual_projector(10, 5, True, prime)
    raw = random_form(10, 5, np.random.default_rng(seed), prime)
    return to_dense((projector @ raw) % prime, 10, 5, prime)


def source_commit():
    try:
        return subprocess.run(
            ["git", "-C", str(ROOT), "rev-parse", "HEAD"],
            capture_output=True, text=True, check=True).stdout.strip()
    except (subprocess.CalledProcessError, OSError):
        return "unknown"


def main(seed_base=41000, extra_samples=8):
    items = registry_items()
    names_by_prime = {}
    out = {}
    started = time.time()

    # A prime subset may be requested so a first pass can establish the rank
    # before the full six-prime certificate is paid for. Checkpoints make the
    # later primes incremental rather than a rerun.
    only = os.environ.get("SDINV_PRIMES")
    primes = tuple(int(p) for p in only.split(",")) if only else PRIMES

    for prime in primes:
        path = CERT / f"interacting_degree12_{prime}.json"
        if not path.exists():
            continue
        with path.open() as stream:
            cert = json.load(stream)
        _, bmap, span = closure_span(cert, SEED8, prime)
        names = bmap[10]
        names_by_prime[prime] = names
        ech, piv = rref(span[10], prime)
        free = [j for j in range(len(names)) if j not in set(piv)]
        n_samples = len(names) + extra_samples

        # Atlas hash is ORDER-SENSITIVE on purpose: atlas coordinates are
        # reported against this exact column order, so a reordering changes
        # the meaning of every stored vector even though the set is identical.
        atlas_sha = hashlib.sha256(
            json.dumps(names).encode()).hexdigest()[:32]
        # Quotient hash covers the projector itself -- the echelon form and the
        # free-column choice. A different D10 closure yields different quotient
        # vectors from identical atlas coordinates, so cached PROJECTIONS must
        # not survive a change here.
        quotient_sha = hashlib.sha256(
            json.dumps({"ech": [[int(v) for v in row] for row in ech],
                        "piv": list(piv), "free": free},
                       sort_keys=True).encode()).hexdigest()[:32]
        store = ProjectionCheckpoint(
            DEFAULT_CKPT / f"prime_{prime}",
            {"source_commit": source_commit(), "atlas_sha256": atlas_sha,
             "quotient_sha256": quotient_sha,
             "block_fingerprint": block_fingerprint(),
             "modular_backend": MODULAR_BACKEND,
             "evaluator_version": EVALUATOR_VERSION, "degree": 10,
             "prime": prime, "seed_base": seed_base})

        forms = [sample(prime, seed_base + 11 * i) for i in range(n_samples)]

        # --- atlas columns, checkpointed ------------------------------------
        atlas_fp = block_fingerprint()
        cols = []
        for j, nm in enumerate(names):
            column = []
            for i, form in enumerate(forms):
                unit = store.load_unit(prime, i, j, fingerprint=atlas_fp)
                if unit is None:
                    t0 = time.time()
                    val = evaluate_atlas_element(items[nm], items, form,
                                                 prime, {})
                    store.save_unit(prime, i, j, nm, val, time.time() - t0,
                                    fingerprint=atlas_fp)
                    column.append(int(val) % prime)
                else:
                    column.append(int(unit["value"]) % prime)
            cols.append(column)
            store.flush_manifest()

        A = [[cols[j][i] for j in range(len(names))] for i in range(n_samples)]

        # --- published candidates, checkpointed -----------------------------
        # The AMB-01/AMB-02 alternative readings are projected alongside the
        # primary ones. They occupy a separate column-id band so a reading can
        # never be confused with the candidate it is a variant of.
        targets = [(cid, spec["evaluator"],
                    FORMULA_ID_BASE + int(cid.split("_")[1]))
                   for cid, spec in sorted(PUBLISHED_DEGREE10.items())]
        targets += [(f"{cid}[{name}]", fn,
                     2 * FORMULA_ID_BASE + int(cid.split("_")[1]))
                    for cid, (amb, name, fn) in sorted(AMBIGUITY_VARIANTS.items())]

        proj = {}
        for cid, evaluator, col_id in targets:
            fp = evaluator_fingerprint(evaluator)
            b = []
            for i, form in enumerate(forms):
                unit = store.load_unit(prime, i, col_id, fingerprint=fp)
                if unit is None:
                    t0 = time.time()
                    val = evaluator(form, prime)
                    store.save_unit(prime, i, col_id, cid, val,
                                    time.time() - t0, fingerprint=fp)
                    b.append(int(val) % prime)
                else:
                    b.append(int(unit["value"]) % prime)
            store.flush_manifest()

            x, ok = solve_exact(A, b, prime)
            if not ok:
                proj[cid] = {"status": "not_in_atlas_span"}
                print(f"      {cid}: NOT IN ATLAS SPAN -- not a Lorentz "
                      f"scalar, or a wrong index placement", flush=True)
                continue
            q = project(x, ech, piv, free, prime)
            proj[cid] = {"status": "solved", "quotient_vector": q,
                         "atlas_coordinates": [int(v) % prime for v in x],
                         "nonzero": any(v % prime for v in q)}

        def _rank(keys, field):
            rows = [proj[k][field] for k in keys
                    if proj[k].get("status") == "solved"]
            return (rank_mod(np.asarray(rows, dtype=np.int64) % prime, prime)
                    if rows else 0)

        primary = [k for k in proj if "[" not in k]
        variants = [k for k in proj if "[" in k]
        # The headline rank is the TWELVE published candidates only. Variant
        # readings are ranked separately -- folding them in would inflate the
        # count with alternative transcriptions of the same equation.
        qr = _rank(primary, "quotient_vector")
        out[prime] = {
            "projections": proj,
            "Q10_rank_from_published": qr,
            "Q10_rank_including_variants": _rank(primary + variants,
                                                 "quotient_vector"),
            "published_atlas_rank": _rank(primary, "atlas_coordinates"),
            "dim_Q10": len(free), "atlas_dim": len(names),
            "role": "fit" if prime in FIT_PRIMES else "holdout"}
        atlas_rank = out[prime]["published_atlas_rank"]
        print(f"  prime {prime} ({out[prime]['role']}): Q10 rank = {qr} / "
              f"{len(free)}, atlas rank = {atlas_rank} / {len(names)}",
              flush=True)
        for cid, p in sorted(proj.items()):
            print(f"      {cid}: {p.get('status')} q={p.get('quotient_vector')}",
                  flush=True)

        _write(out, started)

    _write(out, started)


def _write(out, started):
    """Persist the artifact after every prime, so a kill loses at most one."""
    ranks = {v["Q10_rank_from_published"] for v in out.values()}
    fit_ranks = {v["Q10_rank_from_published"] for p, v in out.items()
                 if p in FIT_PRIMES}
    hold_ranks = {v["Q10_rank_from_published"] for p, v in out.items()
                  if p in HOLDOUT_PRIMES}
    payload = {
        "schema": 2, "degree": 10,
        "implemented": sorted(PUBLISHED_DEGREE10),
        "not_implemented": sorted(NOT_IMPLEMENTED),
        "n_not_implemented": len(NOT_IMPLEMENTED),
        "evaluator_version": EVALUATOR_VERSION,
        "source_commit": source_commit(),
        "fit_primes": list(FIT_PRIMES), "holdout_primes": list(HOLDOUT_PRIMES),
        "per_prime": {str(k): v for k, v in out.items()},
        "consistent": len(ranks) == 1,
        "fit_holdout_agree": (fit_ranks == hold_ranks and len(fit_ranks) == 1),
        "Q10_rank_from_implemented_published":
            ranks.pop() if len(ranks) == 1 else None,
        "runtime_seconds": round(time.time() - started, 1),
        "peak_rss_mb": round(peak_rss_mb(), 1),
    }
    target = ROOT / "results" / "intrinsic_candidates" / "published_degree10_map.json"
    target.parent.mkdir(parents=True, exist_ok=True)
    tmp = target.with_suffix(".json.tmp")
    tmp.write_text(json.dumps(payload, indent=1, sort_keys=True) + "\n")
    os.replace(tmp, target)
    print(f"wrote {target}", flush=True)
    print(f"runtime {payload['runtime_seconds']}s peak_rss "
          f"{payload['peak_rss_mb']} MB", flush=True)


if __name__ == "__main__":
    main()
