#!/usr/bin/env python3
"""Degree-12 reverse pilot. PREPARED, NOT LAUNCHED.

Refuses to run without --i-mean-it, because an unbounded degree-12 sweep is a
multi-day job on this machine and the Goal that produced this file explicitly
forbids launching it automatically.

Reuses the degree-10 machinery with six blocks instead of five. Independence is
the same hard constraint: `sdinv.reverse_block_decomposition` does not import
`published_degree12_invariants`, and nothing here loads the published
structures during generation.

Known input, which is a NEGATIVE result and therefore safe to state up front:
P12_01/02/03 have combined Q12 rank 0 of 4, so all four compact Q12 directions
are unknown and there is no published span for the search to accidentally copy.

Defaults are deliberately small. Measure before forecasting: degree-12
contractions run over six 6-index operands rather than five, so the degree-10
figure of ~150 ms per candidate is an extrapolation, not a measurement. This
project has already recorded one wrong runtime forecast built from too few
points.
"""
import argparse
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

from sdinv.projection_checkpoint import peak_rss_mb
from sdinv.reverse_block_decomposition import (
    build_einsum, evaluate, make_blocks, multiset_slot_count, stream_candidates)

R = ROOT / "results" / "intrinsic_candidates"
PLAN = R / "degree12_reverse_pilot_plan.json"
FIT, HOLDOUT = 32749, 32717


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--i-mean-it", action="store_true",
                    help="required; this pilot does not run by accident")
    ap.add_argument("--sector", default=None,
                    help="default: the first shard named in the plan")
    ap.add_argument("--shard-residue", type=int, default=0)
    ap.add_argument("--shard-modulus", type=int, default=4)
    ap.add_argument("--topology-cap", type=int, default=4000)
    ap.add_argument("--stop-after-candidates", type=int, default=400)
    ap.add_argument("--stop-after-seconds", type=float, default=5400)
    ap.add_argument("--max-rss-mb", type=float, default=1500)
    ap.add_argument("--measure-only", action="store_true",
                    help="time a handful of contractions and exit, so a cost "
                         "model is measured rather than extrapolated")
    args = ap.parse_args()

    plan = json.loads(PLAN.read_text())
    sector = args.sector or plan["first_shard"]["sector"]
    kinds = sector.split("+")
    print(f"degree-12 pilot: sector {sector} ({multiset_slot_count(kinds)} slots)")
    print(f"  shard {args.shard_residue}/{args.shard_modulus}, "
          f"cap {args.topology_cap}, limits: "
          f"{args.stop_after_candidates} candidates / "
          f"{args.stop_after_seconds}s / {args.max_rss_mb}MB")

    if not (args.i_mean_it or args.measure_only):
        print("\nREFUSING TO RUN. This is a prepared plan, not a launch.")
        print("Re-invoke with --measure-only for a cost measurement, or")
        print("--i-mean-it to actually search. See "
              "docs/DEGREE12_REVERSE_PILOT_PLAN.md.")
        return 2

    from sdinv.forms import selfdual_projector, to_dense, random_form
    pr = selfdual_projector(10, 5, True, FIT)
    form = to_dense((pr @ random_form(10, 5, np.random.default_rng(5), FIT))
                    % FIT, 10, 5, FIT)
    t0 = time.time()
    blocks = make_blocks(form, FIT)
    build_s = time.time() - t0

    times, n = [], 0
    for topo, _ in stream_candidates(kinds, cap=args.topology_cap,
                                     shard_residue=args.shard_residue,
                                     shard_modulus=args.shard_modulus,
                                     max_candidates=8 if args.measure_only
                                     else args.stop_after_candidates):
        try:
            build_einsum(kinds, topo)
        except ValueError:
            continue
        t1 = time.time()
        try:
            evaluate(kinds, topo, blocks, FIT)
        except Exception as exc:
            print(f"    evaluation failed: {type(exc).__name__}: {exc}")
            continue
        times.append(time.time() - t1)
        n += 1
        if args.max_rss_mb and peak_rss_mb() > args.max_rss_mb:
            print(f"  ABORT: RSS {peak_rss_mb():.0f}MB > {args.max_rss_mb}")
            break
        if time.time() - t0 > args.stop_after_seconds:
            print("  stopping: --stop-after-seconds reached")
            break

    if times:
        mean = sum(times) / len(times)
        print(f"\n  block build   {build_s:.1f}s")
        print(f"  contraction   {mean*1000:.0f} ms mean over {len(times)} "
              f"(min {min(times)*1000:.0f}, max {max(times)*1000:.0f})")
        print(f"  peak RSS      {peak_rss_mb():.0f} MB")
        print(f"\n  MEASURED, not extrapolated. A full-sector forecast needs "
              f"points spread across the whole range, not this prefix.")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
