#!/usr/bin/env python3
"""Aggregate the rank-matrix cells into one publication artifact.

The aggregator is deliberately hostile to its inputs. Its job is not to produce
a number but to refuse to produce one whenever the cells do not justify it, so
every rejection below corresponds to a way the matrix could be silently wrong:

  missing cell            the matrix is incomplete and must not read as complete
  duplicate cell          one cell counted twice inflates agreement
  mixed candidate order   ranks from different orderings are not comparable
  mixed coordinate dim    likewise
  mixed flop budget       a cheaper budget can silently skip a candidate
  incomplete schedule     82/83 is not 83/83
  evaluation errors       a failed candidate is not an absent one
  Euler failure           the derivative is wrong, whatever the rank says
  zero rows               a zero row cannot contribute and must be explained
  malformed hash          provenance that cannot be checked is not provenance

It never mutates a cell. Reading is the only thing it does to them.

    python scripts/aggregate_rank_matrix.py
    python scripts/aggregate_rank_matrix.py --self-test   # exercise the rejections
"""

from __future__ import annotations

import argparse
import copy
import hashlib
import json
import re
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
CELLS = ROOT / "results" / "rank81" / "cells"
OUT_JSON = ROOT / "results" / "rank81" / "full_rank_matrix_publication_final.json"
OUT_CSV = ROOT / "results" / "rank81" / "full_rank_matrix_publication_final.csv"
OUT_SHA = ROOT / "results" / "rank81" / "full_rank_matrix_publication_final.sha256"

EXPECTED_CANDIDATES = 83
EXPECTED_RANK = 81
HEX64 = re.compile(r"^[0-9a-f]{64}$")


class Rejected(Exception):
    """The cells do not justify a matrix."""


def load_cells(cell_dir: Path) -> list[dict]:
    out = []
    for p in sorted(cell_dir.glob("cell_*.json")):
        try:
            out.append(json.loads(p.read_text()))
        except json.JSONDecodeError as exc:
            raise Rejected(f"{p.name}: not valid JSON ({exc})") from exc
    if not out:
        raise Rejected(f"no cells found under {cell_dir}")
    return out


def validate(cells: list[dict], *, expect_cells: int | None = None) -> dict:
    seen = {}
    orders, dims, budgets = set(), set(), set()

    for c in cells:
        cell = c.get("cell") or {}
        key = (cell.get("prime"), cell.get("seed"))
        if key[0] is None or key[1] is None:
            raise Rejected(f"a cell has no (prime, seed) identity: {cell}")
        if key in seen:
            raise Rejected(f"duplicate cell for prime {key[0]} seed {key[1]}")
        seen[key] = c

        h = c.get("candidate_order_sha256")
        if not (isinstance(h, str) and HEX64.match(h)):
            raise Rejected(f"cell {key}: malformed candidate_order_sha256 {h!r}")
        orders.add(h)
        dims.add(c.get("coordinate_dimension"))
        budgets.add(c.get("flop_limit"))

        s = c.get("schedule_summary") or {}
        e = c.get("euler_homogeneity") or {}
        j = c.get("jacobian") or {}
        planned = s.get("planned")
        evaluated = (s.get("by_terminal_status") or {}).get("evaluated", 0)

        if planned != EXPECTED_CANDIDATES:
            raise Rejected(f"cell {key}: {planned} candidates planned, expected "
                           f"{EXPECTED_CANDIDATES}")
        if evaluated != EXPECTED_CANDIDATES:
            raise Rejected(f"cell {key}: {evaluated}/{planned} evaluated; an "
                           f"incomplete schedule is not a complete one")
        for field in ("evaluation_errors", "interrupted", "structurally_rejected",
                      "zero_rows"):
            if s.get(field, 0) != 0:
                raise Rejected(f"cell {key}: {field} = {s.get(field)}, expected 0")
        if not s.get("complete"):
            raise Rejected(f"cell {key}: schedule not marked complete")
        if e.get("passed") != EXPECTED_CANDIDATES or e.get("failed"):
            raise Rejected(f"cell {key}: Euler homogeneity "
                           f"{e.get('passed')}/{e.get('checked')}, failures "
                           f"{e.get('failed')}")
        if not c.get("cell_complete"):
            raise Rejected(f"cell {key}: cell_complete is false")
        if j.get("total_rank") is None:
            raise Rejected(f"cell {key}: no total_rank recorded")

    if len(orders) != 1:
        raise Rejected(f"cells use {len(orders)} different candidate orderings; "
                       f"their ranks are not comparable")
    if len(dims) != 1:
        raise Rejected(f"cells use {len(dims)} different coordinate dimensions")
    if len(budgets) != 1:
        raise Rejected(f"cells use {len(budgets)} different flop budgets; a "
                       f"smaller budget can silently skip a candidate")

    if expect_cells is not None and len(seen) != expect_cells:
        raise Rejected(f"{len(seen)} cells present, {expect_cells} expected; "
                       f"an incomplete matrix must not be reported as complete")
    return seen


def summarise(seen: dict) -> dict:
    ranks = sorted({c["jacobian"]["total_rank"] for c in seen.values()})
    primes = sorted({k[0] for k in seen})
    seeds = sorted({k[1] for k in seen})

    # A pivot row/column is stable when every cell agrees on it.
    rows = [set(c["jacobian"]["pivot_rows"]) for c in seen.values()]
    cols = [set(c["jacobian"]["pivot_columns"]) for c in seen.values()]
    stable_rows = sorted(set.intersection(*rows)) if rows else []
    stable_cols = sorted(set.intersection(*cols)) if cols else []
    unstable_rows = sorted(set.union(*rows) - set(stable_rows)) if rows else []
    unstable_cols = sorted(set.union(*cols) - set(stable_cols)) if cols else []

    one = next(iter(seen.values()))
    return {
        "cells_complete": len(seen),
        "primes": primes,
        "seeds": seeds,
        "distinct_total_ranks": ranks,
        "all_cells_agree": len(ranks) == 1,
        "rank": ranks[0] if len(ranks) == 1 else None,
        "ranks_by_cell": {f"p{p}_s{s}": seen[(p, s)]["jacobian"]["total_rank"]
                          for (p, s) in sorted(seen)},
        "roles_by_cell": {f"p{p}_s{s}": (seen[(p, s)]["cell"] or {}).get("role")
                          for (p, s) in sorted(seen)},
        "candidate_order_sha256": one["candidate_order_sha256"],
        "coordinate_dimension": one["coordinate_dimension"],
        "flop_limit": one["flop_limit"],
        "stable_pivot_rows": stable_rows,
        "stable_pivot_columns": stable_cols,
        "unstable_pivot_rows": unstable_rows,
        "unstable_pivot_columns": unstable_cols,
        "cumulative_rank_by_degree": one["jacobian"]["cumulative_rank_by_degree"],
        "wall_seconds_total": round(sum(c.get("wall_seconds", 0)
                                        for c in seen.values()), 1),
        "peak_rss_mb_max": max((c.get("peak_rss_mb") or 0) for c in seen.values()),
    }


def render(summary: dict, provenance: dict) -> tuple[str, str]:
    payload = {"schema": 1,
               "generated_by": "scripts/aggregate_rank_matrix.py",
               "provenance": provenance,
               **summary}
    # Deterministic: sorted keys, fixed separators, no timestamp.
    text = json.dumps(payload, indent=1, sort_keys=True) + "\n"
    lines = ["prime,seed,role,rank,candidates_evaluated,euler_passed,zero_rows"]
    for name, rank in summary["ranks_by_cell"].items():
        p, s = name[1:].split("_s")
        lines.append(f"{p},{s},{summary['roles_by_cell'][name]},{rank},"
                     f"{EXPECTED_CANDIDATES},{EXPECTED_CANDIDATES},0")
    return text, "\n".join(lines) + "\n"


def self_test() -> int:
    """Every rejection, exercised. A gate that has never fired is not a gate."""
    base = load_cells(CELLS)
    ok = 0

    def expect_rejected(label, mutate, **kw):
        nonlocal ok
        cells = copy.deepcopy(base)
        mutate(cells)
        try:
            validate(cells, **kw)
        except Rejected as exc:
            print(f"  rejected as expected: {label}\n      {exc}")
            ok += 1
            return
        raise SystemExit(f"NOT REJECTED: {label}")

    expect_rejected("missing cell", lambda c: c.pop(), expect_cells=len(base))
    expect_rejected("duplicate cell", lambda c: c.append(copy.deepcopy(c[0])))
    expect_rejected("mixed candidate ordering",
                    lambda c: c[0].__setitem__("candidate_order_sha256", "f" * 64))
    expect_rejected("mixed coordinate dimension",
                    lambda c: c[0].__setitem__("coordinate_dimension", 125))
    expect_rejected("mixed flop budget",
                    lambda c: c[0].__setitem__("flop_limit", 1.0))
    expect_rejected("82 of 83 candidates",
                    lambda c: c[0]["schedule_summary"]["by_terminal_status"]
                    .__setitem__("evaluated", 82))
    expect_rejected("an evaluation error",
                    lambda c: c[0]["schedule_summary"].__setitem__("evaluation_errors", 1))
    expect_rejected("a zero row",
                    lambda c: c[0]["schedule_summary"].__setitem__("zero_rows", 1))
    expect_rejected("a failed Euler check",
                    lambda c: c[0]["euler_homogeneity"].__setitem__("passed", 82))
    expect_rejected("malformed provenance hash",
                    lambda c: c[0].__setitem__("candidate_order_sha256", "not-a-hash"))
    expect_rejected("cell not marked complete",
                    lambda c: c[0].__setitem__("cell_complete", False))

    # and it must not have touched the originals
    after = load_cells(CELLS)
    if [json.dumps(x, sort_keys=True) for x in base] != \
       [json.dumps(x, sort_keys=True) for x in after]:
        raise SystemExit("NOT READ-ONLY: the aggregator mutated a cell")
    print(f"  {ok} rejections fired; cells unmodified")

    # determinism and order independence
    seen = validate(base)
    a, _ = render(summarise(seen), {})
    shuffled = list(reversed(base))
    b, _ = render(summarise(validate(shuffled)), {})
    if a != b:
        raise SystemExit("NOT DETERMINISTIC: input order changed the output")
    print("  deterministic and input-order independent")
    return 0


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--expect-cells", type=int, default=15)
    ap.add_argument("--self-test", action="store_true")
    args = ap.parse_args()

    if args.self_test:
        return self_test()

    cells = load_cells(CELLS)
    seen = validate(cells, expect_cells=args.expect_cells)
    summary = summarise(seen)

    provenance = {
        "cells_directory": str(CELLS.relative_to(ROOT)),
        "origin": (
            "The cells were produced by a rank-matrix run in a separate working "
            "tree of this project (~/Downloads/sdinv-jhep). They record no "
            "source-critical hash, so their provenance was closed empirically "
            "instead: the two cells that overlap this repository's own "
            "certificate -- (32749, 11) and (32749, 22) -- agree with it on "
            "total rank, dimensions, per-degree block ranks, cumulative ranks, "
            "pivot rows, pivot columns and pivot row candidate ids. All fifteen "
            "share one candidate ordering hash."),
        "not_verified": (
            "The remaining thirteen cells were not independently recomputed "
            "here. What is established is internal consistency plus agreement "
            "with this repository on the two overlapping cells."),
    }
    text, csv = render(summary, provenance)
    OUT_JSON.write_text(text)
    OUT_CSV.write_text(csv)
    digest = hashlib.sha256(text.encode()).hexdigest()
    OUT_SHA.write_text(f"{digest}  {OUT_JSON.name}\n")

    print(f"cells {summary['cells_complete']}/{args.expect_cells} | "
          f"ranks {summary['distinct_total_ranks']} | "
          f"agree {summary['all_cells_agree']}")
    print(f"stable pivot rows {len(summary['stable_pivot_rows'])}, "
          f"columns {len(summary['stable_pivot_columns'])}; "
          f"unstable rows {len(summary['unstable_pivot_rows'])}, "
          f"columns {len(summary['unstable_pivot_columns'])}")
    print(f"sha256 {digest}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
