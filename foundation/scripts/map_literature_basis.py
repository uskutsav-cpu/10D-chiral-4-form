#!/usr/bin/env python3
"""Verify certified source maps, or explicitly reproduce historical fits.

The default reads the corrected source certificate. ``--legacy`` retains the
old fitting experiment with an output path distinct from the canonical maps.
"""

import argparse
from fractions import Fraction
import json
from pathlib import Path
import sys
import subprocess

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))

import numpy as np
from sdinv.catalog import atomic_write_json
from sdinv.certificate import context, digest
from sdinv.contract import planned_value
from sdinv.exactmap import reconstruct_vector, solve_full_column_rank, rank_mod
from sdinv.graphs import graph_from_label
from sdinv.published_degree8_invariants import SELECTED_BASIS as BASIS, SOURCE, evaluate_selected_basis


def rational_rank(matrix):
    a = [list(map(Fraction, row)) for row in matrix]
    rank = 0
    for j in range(len(a[0])):
        pivot = next((i for i in range(rank, len(a)) if a[i][j]), None)
        if pivot is None:
            continue
        a[rank], a[pivot] = a[pivot], a[rank]
        scale = a[rank][j]
        a[rank] = [x / scale for x in a[rank]]
        for i in range(len(a)):
            if i != rank:
                scale = a[i][j]
                a[i] = [x - scale * y for x, y in zip(a[i], a[rank])]
        rank += 1
        if rank == len(a):
            break
    return rank


def inverse(matrix):
    n = len(matrix)
    a = [[Fraction(x) for x in row] + [Fraction(i == j) for j in range(n)]
         for i, row in enumerate(matrix)]
    for j in range(n):
        pivot = next((i for i in range(j, n) if a[i][j]), None)
        if pivot is None:
            raise ValueError("rational matrix is singular")
        a[j], a[pivot] = a[pivot], a[j]
        scale = a[j][j]
        a[j] = [x / scale for x in a[j]]
        for i in range(n):
            if i != j:
                scale = a[i][j]
                a[i] = [x - scale * y for x, y in zip(a[i], a[j])]
    return [row[n:] for row in a]


def reduce_matrix(matrix, prime):
    return np.array([[int(x.numerator * pow(x.denominator, -1, prime) % prime)
                      for x in row] for row in matrix], dtype=np.int64)


def _sample8(prime, seed, graph_items):
    a = np.random.default_rng(seed).integers(0, prime, 126).tolist()
    form, _, _ = context(prime, a)
    graph_values = [planned_value(graph_from_label(g["graph"]), form, 10, 5, True, prime)
                    for g in graph_items]
    quartic_squared = graph_values[0] ** 2 % prime
    graphs = graph_values[1:] + [quartic_squared]
    tensors = evaluate_selected_basis(form, prime) + [quartic_squared]
    return {"seed": seed, "coordinates": a, "graphs": graphs, "tensors": tensors}


def order8(out, cache_dir):
    low = json.loads((ROOT / "results/10d_order8.json").read_text())
    graph_items = [g for g in low["generators"] if g["order"] in (4, 8)]
    primes, holdout = [32749, 32719, 32693, 32713, 32717], 32771
    fit_seeds, fresh_seeds = list(range(2026090700, 2026090709)), [2026090801, 2026090802, 2026090803]
    evidence = {}
    sources = ["src/sdinv/published_degree8_invariants.py", "src/sdinv/stress.py",
               "src/sdinv/contract.py", "src/sdinv/modp.py", "scripts/map_literature_basis.py"]
    identity = digest({"graphs": graph_items, "sources": {s: (ROOT / s).read_bytes().hex() for s in sources}})
    for p in primes + [holdout]:
        cache = cache_dir / f"order8_p{p}.json"
        data = json.loads(cache.read_text()) if cache.exists() else {"identity": identity, "prime": p, "samples": []}
        if data["identity"] != identity:
            raise ValueError("order-8 cache engine mismatch; use a new cache directory")
        seeds = fit_seeds + fresh_seeds if p in primes else fresh_seeds
        if [s["seed"] for s in data["samples"]] != seeds[:len(data["samples"])]:
            raise ValueError("order-8 cache sample mismatch")
        for seed in seeds[len(data["samples"]):]:
            data["samples"].append(_sample8(p, seed, graph_items))
            atomic_write_json(cache, data)
            print(f"order8 p={p}: sample {len(data['samples'])}/{len(seeds)}", flush=True)
        if p in primes:
            design = np.array([s["graphs"] for s in data["samples"][:len(fit_seeds)]])
            target = np.array([s["tensors"] for s in data["samples"][:len(fit_seeds)]])
            data["tensor_to_graph_matrix_mod_p"] = [solve_full_column_rank(design, target[:, j], p).tolist() for j in range(7)]
            data["graph_rank"] = rank_mod(design, p)
            data["tensor_rank"] = rank_mod(target, p)
        evidence[str(p)] = data
    tensor_to_graph = [reconstruct_vector([evidence[str(p)]["tensor_to_graph_matrix_mod_p"][i] for p in primes], primes) for i in range(7)]
    graph_to_tensor = inverse(tensor_to_graph)
    for p in primes + [holdout]:
        forward, backward = reduce_matrix(tensor_to_graph, p), reduce_matrix(graph_to_tensor, p)
        for sample in evidence[str(p)]["samples"]:
            g, t = np.array(sample["graphs"]), np.array(sample["tensors"])
            if not np.array_equal(forward @ g % p, t) or not np.array_equal(backward @ t % p, g):
                raise ValueError(f"order-8 holdout mismatch at prime {p}")
    result = {"schema": 1, "degree": 8, "literature_source": SOURCE,
              "literature_basis": BASIS, "graph_basis": graph_items[1:],
              "basis_selection": "T8_1 through T8_5 from (4.12)--(4.15), and H8_1 from (4.18); this selection spans all seven degree-8 directions after adjoining the quartic square.",
              "product": "I4_1^2 (repository normalization)",
              "matrix_convention": "[G1,...,G6,P]^T = graph_to_literature * [T1,...,T6,P]^T",
              "graph_to_literature": [[str(x) for x in row] for row in graph_to_tensor],
              "literature_to_graph": [[str(x) for x in row] for row in tensor_to_graph],
              "primitive_quotient_matrix_6x6": [[str(x) for x in row[:6]] for row in graph_to_tensor[:6]],
              "product_correction": [str(row[6]) for row in graph_to_tensor[:6]],
              "exact_6x6_without_products": all(row[6] == 0 for row in graph_to_tensor[:6]),
              "fitting_primes": primes, "holdout_prime": holdout,
              "verification_seeds": fresh_seeds, "prime_witnesses": evidence,
              "verified": True,
              "proof_scope": "Exact modular solves, unique bounded rational reconstruction and fresh value holdouts; polynomial identity validity over Q is not proved by finite samples alone."}
    atomic_write_json(out, result)
    print("PASS: order-8 rational map and fresh holdouts")


def order10(out):
    source = ROOT / "results/degree10/B10_coordinates_per_prime.json"
    store = json.loads(source.read_text())["per_prime"]
    primes = sorted(map(int, store))
    fit, holdout = primes[:-1], primes[-1]
    names = store[str(primes[0])]["basis"]
    ids = sorted(store[str(primes[0])]["coordinates"])
    matrices = [reconstruct_vector([store[str(p)]["coordinates"][cid] for p in fit], fit) for cid in ids]
    reduced = reduce_matrix(matrices, holdout)
    if reduced.tolist() != [store[str(holdout)]["coordinates"][cid] for cid in ids]:
        raise ValueError("degree-10 held-out prime does not match rational reconstruction")
    products = [[int(i == j) for i in range(len(names))] for j, name in enumerate(names) if "*" in name]
    b_rank, p_rank, union_rank = rational_rank(matrices), rational_rank(products), rational_rank(matrices + products)
    result = {"schema": 1, "degree": 10,
              "status": "requested invertible 12x12 map does not exist for the implemented source readings",
              "literature_basis": ids, "graph_and_product_basis": names,
              "matrix_convention": "T_i = sum_j matrix[i][j] * atlas_j",
              "matrix_12x14": [[str(x) for x in row] for row in matrices],
              "published_span_rank": b_rank, "product_span_rank": p_rank,
              "union_rank": union_rank, "intersection_rank": b_rank + p_rank - union_rank,
              "primitive_quotient_rank": union_rank - p_rank,
              "fitting_primes": fit, "holdout_prime": holdout,
              "prime_witnesses": store, "verified_reconstruction": True,
              "source_artifact": str(source.relative_to(ROOT)), "source_sha256": digest(json.loads(source.read_text())),
              "source_readings": "See src/sdinv/published_degree10_invariants.py and docs/PUBLISHED_DEGREE10_INDEX_AUDIT.md; AMB-01/AMB-02 readings remain explicit.",
              "proof_scope": "Rational matrix ranks are exact. Inherited polynomial maps were fitted and checked on finite samples; this regeneration checks arithmetic and a held-out prime, not new tensor evaluations."}
    atomic_write_json(out, result)
    print(f"PASS: degree-10 rational map has ranks published={b_rank}, products={p_rank}, union={union_rank}")


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--degree", type=int, choices=[8, 10], required=True)
    parser.add_argument("--out")
    parser.add_argument("--legacy", action="store_true",
                        help="run historical modular fitting, never write a canonical map")
    parser.add_argument("--cache-dir", type=Path, default=ROOT / "work/literature-maps")
    args = parser.parse_args()
    canonical = ROOT / f"results/order{args.degree}_change_of_basis.json"
    protected = {(ROOT / f"results/order{degree}_change_of_basis.json").resolve()
                 for degree in (8, 10)}
    if args.out and Path(args.out).resolve() in protected:
        if args.legacy or Path(args.out).resolve() != canonical.resolve():
            parser.error("output cannot overwrite a canonical map with legacy or different-degree data")
    if not args.legacy:
        subprocess.run([sys.executable, str(ROOT / "scripts/verify_independent_audit.py"),
                        "--source-only"], check=True)
        result = json.loads(canonical.read_text())
        if args.out and Path(args.out).resolve() != canonical.resolve():
            atomic_write_json(args.out, result)
        print(f"PASS: certified degree-{args.degree} source map: {canonical}")
        return
    out = Path(args.out) if args.out else ROOT / f"work/legacy-order{args.degree}-change-of-basis.json"
    if args.degree == 8:
        order8(out, args.cache_dir)
    else:
        order10(out)


if __name__ == "__main__":
    main()
