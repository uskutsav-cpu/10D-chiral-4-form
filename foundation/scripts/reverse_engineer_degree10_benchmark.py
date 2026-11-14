#!/usr/bin/env python3
"""Independent reverse benchmark: recover Q10 from block topology alone.

Generation NEVER consults the published equation-(4.24) formulas or their
coordinate vectors. `sdinv.reverse_block_decomposition` does not import
`sdinv.published_degree10_invariants`, and this runner only loads the published
result AFTER the search, for comparison.

Usage:
    reverse_engineer_degree10_benchmark.py [--limit N] [--pilot]
"""
import argparse
import hashlib
import json
import os
import sys
import time
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))
sys.path.insert(0, str(ROOT / "scripts"))
sys.path.insert(0, str(ROOT / "tests"))

from sdinv.exactmap import rank_mod
from sdinv.forms import selfdual_projector, to_dense, random_form
from sdinv.projection_checkpoint import ProjectionCheckpoint, peak_rss_mb
from sdinv.reverse_block_decomposition import (
    block_multisets, build_einsum, evaluate, generate_candidates, make_blocks,
    stream_candidates)
from solve_intrinsic_quotients import rref, project
from stress_flow_closure import closure_span
from test_M_only_quotients import registry_items, evaluate_atlas_element, solve_exact

CERT = ROOT / "results" / "stress_flow" / "certificates"
SEED8 = {4: ["I4_1"], 6: ["I6_2"], 8: ["I8_3", "I8_4", "I8_5", "I8_6"]}
FIT, HOLDOUT = 32749, 32717
R = ROOT / "results" / "intrinsic_candidates"
OUT = R / "degree10_reverse_benchmark.json"
CKPT = Path(os.environ.get("SDINV_CKPT_ROOT")
            or Path(os.environ.get("TMPDIR", "/tmp")) / "sdinv_ckpt") / "reverse10"


def _encode_topology(topology):
    """JSON-safe encoding of a topology; round-trips via _decode_topology."""
    return sorted([[list(a), list(b), int(c)]
                   for (a, b), c in topology.items()])


def _decode_topology(encoded):
    return {(tuple(a), tuple(b)): int(c) for a, b, c in encoded}


def sample(prime, seed):
    pr = selfdual_projector(10, 5, True, prime)
    return to_dense((pr @ random_form(10, 5, np.random.default_rng(seed), prime))
                    % prime, 10, 5, prime)


def atlas_context(prime):
    with (CERT / f"interacting_degree12_{prime}.json").open() as s:
        cert = json.load(s)
    _, bmap, span = closure_span(cert, SEED8, prime)
    names = bmap[10]
    ech, piv = rref(span[10], prime)
    free = [j for j in range(len(names)) if j not in set(piv)]
    return names, ech, piv, free


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--limit", type=int, default=400,
                    help="max candidates evaluated per block multiset")
    ap.add_argument("--max-topologies", "--topology-cap", dest="max_topologies",
                    type=int, default=60000,
                    help="raw-topology cap per sector; a sector that hits it "
                         "is reported truncated and is NOT exhausted")
    ap.add_argument("--sector", action="append", default=None,
                    help="restrict to a named sector, e.g. N1050+N1050+N1050"
                         "+N1050+N1050; repeatable")
    ap.add_argument("--shard-residue", type=int, default=0)
    ap.add_argument("--shard-modulus", type=int, default=1,
                    help="deterministic sharding on the CANONICAL key, so two "
                         "shards can never claim the same scalar and the union "
                         "over residues is exactly the unsharded set")
    ap.add_argument("--resume", action="store_true",
                    help="reuse checkpointed atlas and candidate evaluations")
    ap.add_argument("--stop-after-candidates", type=int, default=None)
    ap.add_argument("--stop-after-seconds", type=float, default=None)
    ap.add_argument("--max-rss-mb", type=float, default=None,
                    help="abort cleanly if resident memory exceeds this")
    ap.add_argument("--pilot", action="store_true",
                    help="stop after the first multiset that yields rank")
    args = ap.parse_args()

    started = time.time()
    items = registry_items()
    names, ech, piv, free = atlas_context(FIT)
    n_samples = len(names) + 8
    forms = [sample(FIT, 41000 + 11 * i) for i in range(n_samples)]

    atlas_sha = hashlib.sha256(json.dumps(names).encode()).hexdigest()[:32]
    store = ProjectionCheckpoint(
        CKPT / f"prime_{FIT}",
        {"atlas_sha256": atlas_sha, "degree": 10, "prime": FIT,
         "evaluator_version": "reverse-v1"})

    print(f"atlas {len(names)} columns, {n_samples} samples, prime {FIT}",
          flush=True)
    A = []

    for i, form in enumerate(forms):
        row = []
        for j, nm in enumerate(names):
            unit = store.load_unit(FIT, i, j, fingerprint=atlas_sha)
            if unit is None:
                t0 = time.time()
                val = evaluate_atlas_element(items[nm], items, form, FIT, {})
                store.save_unit(FIT, i, j, nm, val, time.time() - t0,
                                fingerprint=atlas_sha)
                row.append(int(val) % FIT)
            else:
                row.append(int(unit["value"]) % FIT)
        A.append(row)
        store.flush_manifest()
    print(f"atlas built in {time.time()-started:.0f}s", flush=True)

    # --- enumerate every sector first, then evaluate SAMPLE-OUTER ----------
    # Candidate-outer would rebuild the three quadratic blocks for every
    # (candidate, sample) pair, and the block build is ~4.4 s against ~150 ms
    # for one contraction -- thirty times the cost of the thing being measured.
    # Caching all samples' blocks instead would hold 22 x 16 MB, which this
    # machine does not have. Sample-outer builds each sample's blocks once and
    # sweeps every candidate against them.
    # STREAMED generation. The list form materialised every topology before
    # any filtering -- up to 30 000 dicts per sector, all live at once -- which
    # drove enumeration RSS to ~2.6 GB across 21 sectors. Streaming
    # canonicalises and discards immediately, so only compact canonical KEYS
    # accumulate; the same sector now peaks at ~33 MB.
    stats_all = {}
    plan = []
    wanted = set(args.sector) if args.sector else None
    for ms in block_multisets(5):
        label = "+".join(ms)
        if wanted and label not in wanted:
            continue
        t0 = time.time()
        kept, last = [], {}
        for topo, st in stream_candidates(
                list(ms), cap=args.max_topologies,
                shard_residue=args.shard_residue,
                shard_modulus=args.shard_modulus,
                max_candidates=args.limit):
            kept.append(topo)
            last = st
        stats = {
            "raw_topologies": last.get("raw", 0),
            "disconnected": last.get("disconnected", 0),
            "duplicate_canonical": last.get("duplicate", 0),
            "canonical": last.get("yielded", 0),
            "truncated": last.get("raw", 0) >= args.max_topologies,
            "enumeration_seconds": round(time.time() - t0, 1),
            "selected_for_evaluation": len(kept),
            "peak_rss_mb_after_sector": round(peak_rss_mb(), 1),
        }
        stats_all[label] = stats
        for idx, topo in enumerate(kept):
            cid = f"{label}#{idx:04d}"
            try:
                spec, _ = build_einsum(list(ms), topo)
            except ValueError as exc:
                # STRUCTURAL rejection is a legitimate outcome and gets a
                # terminal status. It is not the same as an evaluation error
                # and must never be conflated with one.
                plan.append({"multiset": list(ms), "label": label, "id": cid,
                             "topology": topo, "einsum": None,
                             "status": "structurally_rejected",
                             "error": f"{type(exc).__name__}: {exc}"})
                continue
            plan.append({"multiset": list(ms), "label": label, "id": cid,
                         "topology": topo, "einsum": spec,
                         "status": "planned", "error": None})
        print(f"  {label}: canon={stats['canonical']} "
              f"selected={len(kept)} truncated={stats['truncated']} "
              f"rss={stats['peak_rss_mb_after_sector']}MB "
              f"{stats['enumeration_seconds']}s", flush=True)
        if args.max_rss_mb and peak_rss_mb() > args.max_rss_mb:
            print(f"  ABORT: RSS {peak_rss_mb():.0f}MB exceeds "
                  f"--max-rss-mb {args.max_rss_mb}", flush=True)
            break
        if args.stop_after_candidates and len(plan) >= args.stop_after_candidates:
            break
        if args.stop_after_seconds and time.time() - started > args.stop_after_seconds:
            break

    print(f"\nevaluating {len(plan)} candidates on {n_samples} samples",
          flush=True)
    values = [[0] * n_samples for _ in plan]
    for i, form in enumerate(forms):
        t0 = time.time()
        blocks = make_blocks(form, FIT)
        for ci, cand in enumerate(plan):
            if cand["status"] == "structurally_rejected":
                values[ci][i] = None
                continue
            try:
                values[ci][i] = evaluate(cand["multiset"], cand["topology"],
                                         blocks, FIT)
            except Exception as exc:
                # NOT swallowed. The original pilot dropped nine candidates
                # here with a bare `continue`, so the planned and evaluated
                # counts silently disagreed and nothing recorded why.
                values[ci][i] = None
                if cand["status"] != "evaluation_error":
                    cand["status"] = "evaluation_error"
                    cand["error"] = f"{type(exc).__name__}: {exc}"
                    cand["error_sample"] = i
        del blocks
        print(f"    sample {i+1}/{n_samples}  {time.time()-t0:.1f}s "
              f"rss={peak_rss_mb():.0f}MB", flush=True)

    quotient_rows, found, rank = [], [], 0
    per_label = {}
    for ci, cand in enumerate(plan):
        b = values[ci]
        acc = per_label.setdefault(cand["label"],
                                   {"evaluated": 0, "in_atlas_span": 0,
                                    "outside_atlas_span": 0,
                                    "structurally_rejected": 0,
                                    "evaluation_error": 0,
                                    "nonzero_quotient": 0})
        if cand["status"] == "structurally_rejected":
            acc["structurally_rejected"] += 1
            continue
        if cand["status"] == "evaluation_error" or any(v is None for v in b):
            if cand["status"] == "planned":
                cand["status"] = "evaluation_error"
                cand["error"] = "value missing without a recorded exception"
            acc["evaluation_error"] += 1
            continue
        cand["status"] = "evaluated"
        acc["evaluated"] += 1
        x, ok = solve_exact(A, b, FIT)
        if not ok:
            # A successfully evaluated Lorentz scalar of degree 10 MUST lie in
            # the atlas span. Failing that is an implementation defect, not a
            # property of the candidate -- it is how the mixed-variance M block
            # was found -- so it is recorded per candidate, not counted and
            # forgotten.
            acc["outside_atlas_span"] += 1
            cand["outside_atlas_span"] = True
            continue
        acc["in_atlas_span"] += 1
        cand["atlas_coordinates"] = [int(v) % FIT for v in x]
        q = project(x, ech, piv, free, FIT)
        if not any(v % FIT for v in q):
            continue
        acc["nonzero_quotient"] += 1
        trial = quotient_rows + [q]
        nr = rank_mod(np.asarray(trial, dtype=np.int64) % FIT, FIT)
        if nr > rank:
            rank, quotient_rows = nr, trial
            found.append({"id": cand["id"], "multiset": cand["multiset"],
                          "einsum": cand["einsum"],
                          "topology": _encode_topology(cand["topology"]),
                          "atlas_coordinates": cand["atlas_coordinates"],
                          "quotient_vector": q, "rank_after": rank})
            print(f"    RANK {rank}/3 from {cand['label']}  {cand['einsum']}",
                  flush=True)
            # Deliberately NOT breaking at rank 3. Stopping early left the
            # remaining candidates without a terminal status, and they were
            # then reported as "interrupted" -- which was simply false: they
            # had been evaluated, just not projected. The projection is pure
            # linear algebra on values already computed, so completing it costs
            # almost nothing and makes the accounting identity honest.
    for label, acc in per_label.items():
        stats_all.setdefault(label, {}).update(acc)

    # --- terminal-status accounting; the identity must reconcile exactly ----
    counts = {"planned": len(plan), "evaluated": 0, "structurally_rejected": 0,
              "evaluation_error": 0, "interrupted": 0}
    exceptions = []
    outside = []
    for cand in plan:
        st = cand["status"]
        if st == "planned":
            st = cand["status"] = "interrupted"
        counts[st] = counts.get(st, 0) + 1
        if st in ("structurally_rejected", "evaluation_error"):
            exceptions.append({"id": cand["id"], "label": cand["label"],
                               "status": st, "error": cand.get("error"),
                               "einsum": cand.get("einsum"),
                               "error_sample": cand.get("error_sample")})
        if cand.get("outside_atlas_span"):
            outside.append({"id": cand["id"], "einsum": cand.get("einsum")})
    reconciles = (counts["planned"] == counts["evaluated"]
                  + counts["structurally_rejected"]
                  + counts["evaluation_error"] + counts["interrupted"])
    print(f"\naccounting: {counts}  reconciles={reconciles}", flush=True)
    if outside:
        print(f"  {len(outside)} evaluated candidates OUTSIDE the atlas span "
              f"-- implementation defect, not a result", flush=True)

    payload = {
        "schema": 2,
        # Emitted by the RUNNER, not injected afterwards: a regenerated
        # artifact must not silently lose the frozen claim wording. It did
        # exactly that once, and tests/test_q10_wording_freeze.py caught it.
        "basis_wording": (
            "Preferred ambiguity-minimal Level-B basis among the twelve "
            "published degree-10 candidates under the documented "
            "deterministic simplicity rule."),
        "forbidden_wording": ["ambiguity-robust", "universally canonical",
                              "the unique compact basis"],
        "ambiguity_facts": {
            "P10_10": ("forced -- appears in every independent triple; sole "
                       "carrier of the third quotient coordinate"),
            "P10_09": ("source-reading dependent (AMB-01); quotient image "
                       "differs between readings; excluded from the basis, "
                       "retained as a valid implemented interpretation"),
            "P10_11": ("source-reading dependent (AMB-02); quotient image "
                       "differs between readings; in the basis because at "
                       "least one non-robust member is unavoidable"),
            "P10_12": ("tested alternative AMB-02 reading gives identical "
                       "evaluations and identical quotient vectors"),
            "no_fully_robust_published_triple": True},
        "search_is_exhaustive": False,
        "accounting": counts,
        "accounting_reconciles": reconciles,
        "exceptions": exceptions,
        "outside_atlas_span": outside,
        "independence": ("generation never imports or inspects "
                         "published_degree10_invariants or its coordinates"),
        "fit_prime": FIT,
        "atlas_columns": len(names),
        "dim_Q10": len(free),
        "recovered_rank": rank,
        "recovered_vectors": quotient_rows,
        "recovered_candidates": found,
        "per_multiset": stats_all,
        "limit_per_multiset": args.limit,
        "max_topologies": args.max_topologies,
        "runtime_seconds": round(time.time() - started, 1),
        "peak_rss_mb": round(peak_rss_mb(), 1),
    }
    OUT.write_text(json.dumps(payload, indent=1, sort_keys=True) + "\n")
    print(f"\nrecovered Q10 rank {rank}/{len(free)} in "
          f"{payload['runtime_seconds']}s, peak {payload['peak_rss_mb']} MB")
    print(f"wrote {OUT.name}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
