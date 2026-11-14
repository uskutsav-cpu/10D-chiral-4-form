"""Self-contained exact Jacobian witnesses for the ordered graph basis.

All directions and polynomials are integral. A nonzero minor modulo one
prime therefore proves a characteristic-zero lower bound. Recomputing graph
contractions is a separate verification step from checking the saved matrix.
"""

import hashlib
import json
from math import isqrt

import numpy as np

from .contract import CompactDerivativeBasis, value_and_jacobian_row
from .forms import basis_tuples, hodge_matrix, to_dense
from .graphs import validate_graph
from .modp import RankSieve
from .serialize import from_record


def digest(data):
    return hashlib.sha256(json.dumps(data, sort_keys=True, separators=(",", ":")).encode()).hexdigest()


def require_prime(prime):
    # These modest primes also keep all compact linear algebra below int64.
    if type(prime) is not int or not 13 < prime <= 65521:
        raise ValueError("certificate prime must lie between 13 and 65521")
    if any(prime % d == 0 for d in range(2, isqrt(prime) + 1)):
        raise ValueError("certificate modulus is not prime")


def integral_directions():
    """Columns e_I + *e_I for lexicographic I containing the time index 0.

    The first 126 sorted components form an identity matrix. Thus the
    coordinates are exactly F_{0abcd}, with no 1/2 projector denominators.
    """
    hodge = hodge_matrix(10, 5, True, 32749)
    hodge = np.where(hodge > 32749 // 2, hodge - 32749, hodge)
    indices = [k for k, t in enumerate(basis_tuples(10, 5)) if t[0] == 0]
    directions = (np.eye(252, dtype=np.int64) + hodge)[:, indices]
    if not np.array_equal(directions[indices], np.eye(126, dtype=np.int64)):
        raise ValueError("unexpected self-dual coordinate convention")
    if not np.array_equal(hodge @ directions, directions):
        raise ValueError("integral directions are not self-dual")
    return directions, indices


def context(prime, coordinates):
    require_prime(prime)
    directions, indices = integral_directions()
    a = np.asarray(coordinates)
    if a.dtype.kind not in "iu":
        raise ValueError("point coordinates must be integers")
    if a.shape != (126,) or np.any(a < 0) or np.any(a >= prime):
        raise ValueError("point must contain 126 reduced coordinates")
    a = a.astype(np.int64)
    compact = directions @ a % prime
    basis = CompactDerivativeBasis(10, 5, directions, indices, prime)
    return to_dense(compact, 10, 5, prime), basis, compact


def determinant_mod(matrix, prime):
    """Python-integer Gaussian determinant; independent of RankSieve."""
    require_prime(prime)
    rows = [[int(x) % prime for x in row] for row in matrix]
    n = len(rows)
    if any(len(row) != n for row in rows):
        raise ValueError("determinant requires a square matrix")
    result = 1
    for k in range(n):
        pivot = next((i for i in range(k, n) if rows[i][k]), None)
        if pivot is None:
            return 0
        if pivot != k:
            rows[k], rows[pivot] = rows[pivot], rows[k]
            result = -result
        result = result * rows[k][k] % prime
        inverse = pow(rows[k][k], -1, prime)
        for i in range(k + 1, n):
            scale = rows[i][k] * inverse % prime
            for j in range(k + 1, n):
                rows[i][j] = (rows[i][j] - scale * rows[k][j]) % prime
            rows[i][k] = 0
    return result % prime


def validate_basis(manifest):
    items = manifest["invariants"]
    if len(items) != 81 or len({x["id"] for x in items}) != 81:
        raise ValueError("certificate requires exactly 81 distinct invariant IDs")
    if [x["degree"] for x in items] != [4, 6, 6] + [8] * 6 + [10] * 12 + [12] * 60:
        raise ValueError("unexpected ordered degree inventory")
    for item in items:
        graph = from_record(item["graph"])
        validate_graph(graph, valence=5, max_mult=4)
        if len(graph) != item["degree"] or item["formula_sha256"] != digest(item["graph"]):
            raise ValueError("graph degree or formula hash mismatch")
        from .graphs import graph_from_label
        from .latex import graph_to_latex
        if item["adjacency_matrix"] != graph.tolist() or not np.array_equal(graph_from_label(item["label"]), graph):
            raise ValueError("graph representations disagree")
        if item["einstein_formula"] != graph_to_latex(graph):
            raise ValueError("Einstein formula does not match graph")
    return items


def evaluate_item(item, prime, coordinates, max_memory_bytes=2 * 1024**3):
    form, basis, _ = context(prime, coordinates)
    scalar, row = value_and_jacobian_row(from_record(item["graph"]), form, basis,
                                        10, 5, True, prime,
                                        max_memory_bytes=max_memory_bytes)
    if int(row @ np.asarray(coordinates, dtype=np.int64) % prime) != item["degree"] * scalar % prime:
        raise ValueError(f"Euler homogeneity failed for {item['id']}")
    return {"id": item["id"], "degree": item["degree"], "value": int(scalar),
            "jacobian_row": row.tolist()}


def make_witness(manifest, prime, seed, coordinates, evaluated):
    items = validate_basis(manifest)
    if [r["id"] for r in evaluated] != [x["id"] for x in items]:
        raise ValueError("evaluated rows do not match the fixed basis")
    form, _, compact = context(prime, coordinates)
    rows = np.asarray([r["jacobian_row"] for r in evaluated], dtype=np.int64)
    if rows.shape != (81, 126):
        raise ValueError("Jacobian has wrong shape")
    sieve, ranks = RankSieve(126, prime), {}
    for item, row in zip(items, rows):
        if not sieve.add(row):
            raise ValueError(f"fixed basis loses rank at {item['id']}")
        ranks[str(item["degree"])] = sieve.rank
    columns = sieve.pivots
    determinant = determinant_mod(rows[:, columns], prime)
    if not determinant:
        raise ValueError("minor is zero")
    return {"prime": prime, "seed": seed, "basis_sha256": digest(manifest),
            "coordinates": list(map(int, coordinates)), "selfdual_components": compact.tolist(),
            "invariant_ids": [x["id"] for x in items], "values": [r["value"] for r in evaluated],
            "jacobian": rows.tolist(), "pivot_columns": columns,
            "determinant_mod_p": determinant, "rank": sieve.rank,
            "cumulative_rank_by_degree": ranks, "euler_homogeneity": True}


def verify_witness(manifest, witness, recompute=False):
    items = validate_basis(manifest)
    prime = witness["prime"]
    require_prime(prime)
    if witness["basis_sha256"] != digest(manifest):
        raise ValueError("basis hash mismatch")
    if witness["invariant_ids"] != [x["id"] for x in items]:
        raise ValueError("invariant ordering mismatch")
    coordinates = witness["coordinates"]
    _, _, compact = context(prime, coordinates)
    if witness["selfdual_components"] != compact.tolist():
        raise ValueError("self-dual point mismatch")
    rows = np.asarray(witness["jacobian"])
    values = np.asarray(witness["values"])
    if rows.shape != (81, 126) or rows.dtype.kind not in "iu":
        raise ValueError("Jacobian must contain 81 by 126 integers")
    if values.shape != (81,) or values.dtype.kind not in "iu":
        raise ValueError("values must contain 81 integers")
    if np.any(rows < 0) or np.any(rows >= prime) or np.any(values < 0) or np.any(values >= prime):
        raise ValueError("values and rows must be reduced modulo the prime")
    degrees = np.asarray([x["degree"] for x in items], dtype=np.int64)
    if not np.array_equal(rows @ np.asarray(coordinates) % prime, degrees * values % prime):
        raise ValueError("Euler homogeneity mismatch")
    if witness.get("euler_homogeneity") is not True:
        raise ValueError("Euler validation metadata mismatch")
    columns = witness["pivot_columns"]
    if len(columns) != 81 or len(set(columns)) != 81 or any(type(c) is not int or not 0 <= c < 126 for c in columns):
        raise ValueError("invalid pivot columns")
    determinant = determinant_mod(rows[:, columns], prime)
    if determinant == 0 or determinant != witness["determinant_mod_p"]:
        raise ValueError("minor determinant mismatch")
    sieve, ranks = RankSieve(126, prime), {}
    for item, row in zip(items, rows):
        sieve.add(row)
        ranks[str(item["degree"])] = sieve.rank
    if witness["rank"] != 81 or ranks != witness["cumulative_rank_by_degree"]:
        raise ValueError("rank metadata mismatch")
    if recompute:
        for item, saved_row, saved_value in zip(items, rows, values):
            fresh = evaluate_item(item, prime, coordinates)
            if fresh["value"] != saved_value or fresh["jacobian_row"] != saved_row.tolist():
                raise ValueError(f"graph recomputation mismatch: {item['id']}")
    return {"prime": prime, "seed": witness["seed"], "rank": 81,
            "determinant_mod_p": determinant, "graph_rows_recomputed": recompute}
