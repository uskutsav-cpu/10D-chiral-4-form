#!/usr/bin/env python3
"""Checkpointed degree runner with separate polynomial and functional sieves.

By default, revalidate the discovered graph inventory at the chosen degree.
For new searches, --candidates accepts edge-list JSONL records, optionally
wrapped as {"id": ..., "graph": ...}. Nothing is exhaustively enumerated by
default. Exact nauty generation remains in the existing degree pipelines.
"""

import argparse
from concurrent.futures import ProcessPoolExecutor
import hashlib
import json
import os
from pathlib import Path
import sys

for name in ("OPENBLAS_NUM_THREADS", "OMP_NUM_THREADS", "VECLIB_MAXIMUM_THREADS"):
    os.environ.setdefault(name, "1")

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))

import numpy as np
from sdinv.catalog import atomic_write_json
from sdinv.certificate import digest, evaluate_item, require_prime
from sdinv.checkpoint import load_checkpoint, write_checkpoint
from sdinv.graphs import graph_from_label, validate_graph
from sdinv.modp import RankSieve
from sdinv.primitive import HILBERT_DIMENSIONS, HILBERT_FACTOR_EXPONENTS, PolynomialSieve, product_monomials
from sdinv.serialize import loads, to_record


def inventory():
    result = []
    for degree in (8, 10, 12):
        for item in json.loads((ROOT / f"results/10d_order{degree}.json").read_text())["generators"]:
            result.append({"id": item["id"], "degree": item["order"],
                           "graph": to_record(graph_from_label(item["graph"]))})
    return result


def candidate_stream(path, known, degree):
    if path is None:
        yield from (item for item in known if item["degree"] == degree)
        return
    with Path(path).open() as stream:
        for line in stream:
            if not line.strip():
                continue
            record = json.loads(line)
            matrix = loads(record.get("graph", record))
            validate_graph(matrix, 5, 4)
            if len(matrix) != degree:
                raise ValueError("candidate degree mismatch")
            graph = to_record(matrix)
            yield {"id": record.get("id", digest(graph)), "degree": degree, "graph": graph}


def _evaluate(job):
    item, prime, points = job
    return {"item": item, "samples": [evaluate_item(item, prime, point) for point in points]}


def _product_gradient(factors, cache, sample, prime):
    scalar, row = 1, np.zeros(126, dtype=np.int64)
    for factor in factors:
        value = cache[factor][sample]["value"]
        derivative = np.array(cache[factor][sample]["jacobian_row"])
        row = (row * value + scalar * derivative) % prime
        scalar = scalar * value % prime
    return row


def run(degree, prime, seeds, candidates=None, checkpoint=None, workers=1, max_candidates=None):
    require_prime(prime)
    if degree <= 0 or degree % 2 or degree * 5 // 2 > 52:
        raise ValueError("positive even degree with at most 52 contracted pairs required")
    if len(seeds) < 2 or len(set(seeds)) != len(seeds) or workers < 1:
        raise ValueError("at least two distinct sample seeds and positive workers required")
    known = inventory()
    lower = [x for x in known if x["degree"] < degree]
    # Above degree 14, the committed lower-degree inventory is no longer
    # polynomially complete. Rank evidence is still useful, but do not call
    # a candidate primitive relative to uncomputed lower generators.
    lower_polynomial_inventory_complete = degree <= 14
    if degree > 12 and candidates is None:
        raise ValueError("supply --candidates JSONL for degrees above the discovered inventory")
    points = [np.random.default_rng(s).integers(0, prime, 126).tolist() for s in seeds]
    sources = [ROOT / "src/sdinv" / f"{name}.py" for name in ("certificate", "contract", "forms", "modp", "primitive", "serialize", "graphs")]
    sources.append(Path(__file__))
    identity = {"degree": degree, "prime": prime, "seeds": seeds, "points": points,
                "inventory_sha256": digest(known), "engine_sha256": digest([hashlib.sha256(p.read_bytes()).hexdigest() for p in sources]),
                "candidates_sha256": hashlib.sha256(Path(candidates).read_bytes()).hexdigest() if candidates else None}
    state = load_checkpoint(checkpoint, identity) if checkpoint and Path(checkpoint).exists() else {"lower": [], "accepted": [], "cursor": 0}
    if [r["item"] for r in state["lower"]] != lower[:len(state["lower"])]:
        raise ValueError("checkpoint lower inventory mismatch")

    def save():
        if checkpoint:
            write_checkpoint(checkpoint, identity, state)

    # Lower computations are also resumable; never lose expensive warmup rows.
    for item in lower[len(state["lower"]):]:
        state["lower"].append(_evaluate((item, prime, points)))
        save()
    functional = RankSieve(126, prime)
    cache = {r["item"]["id"]: r["samples"] for r in state["lower"]}
    for r in state["lower"]:
        functional.add(r["samples"][0]["jacobian_row"])
    polynomial = PolynomialSieve(degree, 126 * len(points), prime, "stacked_gradients")
    products = product_monomials({x["id"]: x["degree"] for x in lower}, degree)
    for factors in products:
        row = np.concatenate([_product_gradient(factors, cache, k, prime) for k in range(len(points))])
        polynomial.add("*".join(factors), row, product=True)
    for r in state["accepted"]:
        polynomial.add(r["item"]["id"], np.concatenate([s["jacobian_row"] for s in r["samples"]]))
        functional.add(r["samples"][0]["jacobian_row"])
    stream = iter(candidate_stream(candidates, known, degree))
    for _ in range(state["cursor"]):
        next(stream)
    pool = ProcessPoolExecutor(max_workers=workers) if workers > 1 else None
    try:
        while not polynomial.complete and (max_candidates is None or state["cursor"] < max_candidates):
            batch = []
            for _ in range(workers):
                if max_candidates is not None and state["cursor"] + len(batch) >= max_candidates:
                    break
                item = next(stream, None)
                if item is None:
                    break
                batch.append((item, prime, points))
            if not batch:
                break
            for evaluated in (pool.map(_evaluate, batch) if pool else map(_evaluate, batch)):
                state["cursor"] += 1
                item, samples = evaluated["item"], evaluated["samples"]
                primitive = polynomial.add(item["id"], np.concatenate([s["jacobian_row"] for s in samples]))
                functional_increment = functional.add(samples[0]["jacobian_row"])
                if functional_increment and not primitive:
                    raise ValueError("functional pivot without a polynomial pivot")
                if primitive:
                    evaluated.update({"functional_increment": functional_increment, "polynomial_rank": polynomial.rank,
                                      "functional_rank": functional.rank, "functional_pivots": list(functional.pivots),
                                      "polynomial_pivots": list(polynomial.sieve.pivots)})
                    state["accepted"].append(evaluated)
                    print(f"degree {degree} p={prime}: polynomial {polynomial.rank}, functional {functional.rank}", flush=True)
                save()
                if polynomial.complete:
                    break
    finally:
        if pool:
            pool.shutdown(wait=True, cancel_futures=True)
    return {"schema": 1, "degree": degree, "prime": prime, "seeds": seeds, "coordinates": points,
            "polynomial_method": "gradients stacked across independently sampled points",
            "polynomial_rank": polynomial.rank, "product_rank": len(polynomial.product_ids),
            "primitive_directions_found": len(polynomial.primitive_ids), "functional_rank": functional.rank,
            "primitive_count_certified_relative_to_all_lower_degrees": lower_polynomial_inventory_complete,
            "hilbert_dimension": HILBERT_DIMENSIONS.get(degree),
            "hilbert_factor_exponent": HILBERT_FACTOR_EXPONENTS.get(degree),
            "polynomial_complete": polynomial.complete, "candidates_evaluated": state["cursor"],
            "accepted": state["accepted"], "identity": identity,
            "scope": "Rank lower bounds from fixed exact samples. Polynomial completeness uses the supplied Hilbert upper bound; a rejected row alone is not a symbolic relation. Above degree 14, candidate pivots are only relative to the supplied, incomplete lower inventory."}


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--degree", type=int, required=True)
    parser.add_argument("--primes", type=int, nargs="+", default=[32749, 32719])
    parser.add_argument("--seeds", type=int, nargs="+", default=[20260907, 20260908])
    parser.add_argument("--candidates", type=Path)
    parser.add_argument("--checkpoint-dir", type=Path, default=ROOT / "work/degrees")
    parser.add_argument("--workers", type=int, default=1)
    parser.add_argument("--max-candidates", type=int)
    parser.add_argument("--out", type=Path)
    args = parser.parse_args(argv)
    results = []
    for prime in args.primes:
        results.append(run(args.degree, prime, args.seeds, args.candidates,
                           args.checkpoint_dir / f"degree{args.degree}_p{prime}.json",
                           args.workers, args.max_candidates))
    selection = [[r["item"]["graph"] for r in result["accepted"]] for result in results]
    same_basis = all(basis == selection[0] for basis in selection)
    payload = {"runs": results, "same_selected_basis": same_basis,
               "complete": same_basis and all(r["polynomial_complete"] for r in results)}
    atomic_write_json(args.out or ROOT / f"results/degree{args.degree}_verification.json", payload)
    if not payload["complete"]:
        print("INCOMPLETE: search evidence saved; no degree-completeness claim")
        return 2
    print("PASS: same selected basis, independently checked polynomial and functional ranks")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
