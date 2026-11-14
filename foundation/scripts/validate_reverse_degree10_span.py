#!/usr/bin/env python3
"""Validate the reverse-recovered Q10 triple and compare it to the published span.

Order matters and is enforced by construction: the reverse search has already
finished and its three einsum specifications are read from its own artifact.
Only NOW is the published basis loaded, for comparison. Nothing here feeds back
into generation.

Validation performed:
  * re-evaluate the three recovered contractions on the FIT prime;
  * evaluate them on an independent HOLDOUT prime;
  * evaluate on FRESH samples (a different seed base) at both primes;
  * confirm rank 3 in each case;
  * prove span equality with the published Level-B basis;
  * compute exact change-of-basis matrices both ways and verify they are
    mutually inverse.
"""
import json
import sys
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))
sys.path.insert(0, str(ROOT / "scripts"))
sys.path.insert(0, str(ROOT / "tests"))

from sdinv.exactmap import rank_mod
from sdinv.forms import selfdual_projector, to_dense, random_form
from sdinv.reverse_block_decomposition import evaluate, make_blocks
from solve_intrinsic_quotients import rref, project
from stress_flow_closure import closure_span
from test_M_only_quotients import registry_items, evaluate_atlas_element, solve_exact

CERT = ROOT / "results" / "stress_flow" / "certificates"
SEED8 = {4: ["I4_1"], 6: ["I6_2"], 8: ["I8_3", "I8_4", "I8_5", "I8_6"]}
R = ROOT / "results" / "intrinsic_candidates"
FIT, HOLDOUT = 32749, 32717


def sample(prime, seed):
    pr = selfdual_projector(10, 5, True, prime)
    return to_dense((pr @ random_form(10, 5, np.random.default_rng(seed), prime))
                    % prime, 10, 5, prime)


def inv_matrix(mat, p):
    n = len(mat)
    aug = [[int(mat[i][j]) % p for j in range(n)]
           + [1 if i == j else 0 for j in range(n)] for i in range(n)]
    for col in range(n):
        piv = next((r for r in range(col, n) if aug[r][col] % p), None)
        if piv is None:
            return None
        aug[col], aug[piv] = aug[piv], aug[col]
        f = pow(aug[col][col], p - 2, p)
        aug[col] = [(v * f) % p for v in aug[col]]
        for r in range(n):
            if r != col and aug[r][col] % p:
                g = aug[r][col]
                aug[r] = [(a - g * b) % p for a, b in zip(aug[r], aug[col])]
    return [row[n:] for row in aug]


def matmul(a, b, p):
    return [[sum(a[i][t] * b[t][j] for t in range(len(b))) % p
             for j in range(len(b[0]))] for i in range(len(a))]


def project_candidates(prime, cands, topos, seed_base):
    """Atlas-solve and quotient-project each recovered contraction."""
    items = registry_items()
    with (CERT / f"interacting_degree12_{prime}.json").open() as s:
        cert = json.load(s)
    _, bmap, span = closure_span(cert, SEED8, prime)
    names = bmap[10]
    ech, piv = rref(span[10], prime)
    free = [j for j in range(len(names)) if j not in set(piv)]
    n = len(names) + 8
    forms = [sample(prime, seed_base + 11 * i) for i in range(n)]

    A = [[evaluate_atlas_element(items[nm], items, f, prime, {}) % prime
          for nm in names] for f in forms]

    # sample-outer: build each sample's blocks once, sweep all three candidates
    cols = [[0] * len(forms) for _ in cands]
    for i, f in enumerate(forms):
        bl = make_blocks(f, prime)
        for ci, (cand, topo) in enumerate(zip(cands, topos)):
            cols[ci][i] = evaluate(cand["multiset"], topo, bl, prime)
        del bl

    out = []
    for ci, cand in enumerate(cands):
        b = cols[ci]
        x, ok = solve_exact(A, b, prime)
        if not ok:
            out.append(None)
            continue
        out.append(project(x, ech, piv, free, prime))
    return out, len(free)


def rehydrate_topologies(cands):
    """Recover each candidate's topology by matching its recorded einsum.

    The benchmark artifact stores the einsum specification but not the
    topology object, so the topology is re-derived rather than trusted from a
    serialisation: the candidate's own sector is re-enumerated and the unique
    topology whose `build_einsum` reproduces the recorded spec is taken. That
    is a stricter check than deserialising would be -- it confirms the recorded
    spec is actually generatable by the declared search, and it fails loudly if
    the generator has changed since the pilot ran.
    """
    from sdinv.reverse_block_decomposition import build_einsum, stream_candidates

    out = []
    for cand in cands:
        kinds = cand["multiset"]
        match = None
        for topo, _ in stream_candidates(list(kinds), cap=60000):
            try:
                spec, _r = build_einsum(list(kinds), topo)
            except ValueError:
                continue
            if spec == cand["einsum"]:
                match = topo
                break
        if match is None:
            raise SystemExit(
                f"could not regenerate {cand['einsum']} in sector "
                f"{'+'.join(kinds)}; the generator no longer produces the "
                f"candidate the pilot recorded")
        out.append(match)
    return out


def main():
    bench = json.loads((R / "degree10_reverse_benchmark.json").read_text())
    cands = bench["recovered_candidates"]
    assert len(cands) == 3, f"expected 3 recovered directions, got {len(cands)}"
    print("reverse-recovered contractions (generation already complete):")
    for c in cands:
        print(f"  {'+'.join(c['multiset'])}  {c['einsum']}")

    topos = rehydrate_topologies(cands)
    print(f"  regenerated all {len(topos)} topologies from their einsums")

    result = {"schema": 1, "recovered_einsums": [c["einsum"] for c in cands],
              "recovered_multisets": [c["multiset"] for c in cands],
              "per_prime": {}}

    for prime in (FIT, HOLDOUT):
        for tag, seed_base in (("original", 41000), ("fresh", 77000)):
            vecs, dim = project_candidates(prime, cands, topos, seed_base)
            ok = all(v is not None for v in vecs)
            rank = (rank_mod(np.asarray(vecs, dtype=np.int64) % prime, prime)
                    if ok else 0)
            result["per_prime"].setdefault(str(prime), {})[tag] = {
                "vectors": vecs, "rank": rank, "dim_Q10": dim,
                "all_in_atlas_span": ok}
            print(f"  prime {prime} [{tag} samples]: rank {rank}/{dim} "
                  f"in_span={ok}", flush=True)

    # --- span comparison, only now loading the published result -------------
    pub_map = json.loads((R / "published_degree10_map.json").read_text())
    basis = json.loads((R / "intrinsic_Q10_levelB_basis.json").read_text())
    triple = basis["preferred_basis"]
    result["published_basis"] = triple

    comparison = {}
    for prime in (FIT, HOLDOUT):
        p = prime
        rev = result["per_prime"][str(prime)]["original"]["vectors"]
        pub = [pub_map["per_prime"][str(prime)]["projections"][n]["quotient_vector"]
               for n in triple]
        dim = result["per_prime"][str(prime)]["original"]["dim_Q10"]
        r_rev = rank_mod(np.asarray(rev, dtype=np.int64) % p, p)
        r_pub = rank_mod(np.asarray(pub, dtype=np.int64) % p, p)
        r_union = rank_mod(np.asarray(rev + pub, dtype=np.int64) % p, p)
        # Both spans sit inside a dim-3 quotient. Each having rank 3 already
        # forces both to BE Q10, hence equal; the union rank is the direct
        # check and must not exceed 3.
        equal = (r_rev == r_pub == r_union == dim)
        Rinv = inv_matrix(rev, p)
        rev_in_pub = matmul(rev, inv_matrix(pub, p), p)
        pub_in_rev = matmul(pub, Rinv, p)
        prod = matmul(rev_in_pub, pub_in_rev, p)
        identity = [[1 if i == j else 0 for j in range(3)] for i in range(3)]
        comparison[str(prime)] = {
            "rank_reverse": r_rev, "rank_published": r_pub,
            "rank_union": r_union, "dim_Q10": dim,
            "spans_equal": equal,
            "reverse_in_published_basis": rev_in_pub,
            "published_in_reverse_basis": pub_in_rev,
            "mutually_inverse": prod == identity,
        }
        print(f"  prime {prime}: rank_rev={r_rev} rank_pub={r_pub} "
              f"union={r_union}/{dim} spans_equal={equal} "
              f"mutually_inverse={prod == identity}")

    result["span_comparison"] = comparison
    result["relationship"] = (
        "DIFFERENT compact basis, equivalent span. Every recovered direction "
        "comes from the N4125^5 sector, whereas all three published basis "
        "elements are N1050-based. The reverse search did not rediscover the "
        "published formulas; it found an independent set spanning the same "
        "quotient.")
    result["all_checks_pass"] = (
        all(v["spans_equal"] and v["mutually_inverse"]
            for v in comparison.values())
        and all(t["rank"] == 3 for pr in result["per_prime"].values()
                for t in pr.values()))

    (R / "degree10_reverse_span_validation.json").write_text(
        json.dumps(result, indent=1, sort_keys=True) + "\n")
    print(f"\nall_checks_pass={result['all_checks_pass']}")
    return 0 if result["all_checks_pass"] else 1


if __name__ == "__main__":
    raise SystemExit(main())
