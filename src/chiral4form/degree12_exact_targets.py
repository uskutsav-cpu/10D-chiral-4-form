"""Direct characteristic-zero recovery of all 36 reduced degree-12 source targets.

This module deliberately bypasses the superseded Step1C raw-target height
shortcut.  That shortcut bounded coefficient monomials *before* the exact
six-parameter degree-10 orbit substitution; consequently some nonzero reduced
rows (notably B4/B5/B7) carried a raw bound of zero.

The replacement proof is fail-closed:

1. Reuse the eight expensive 72x72 degree-12 basis checkpoints and reconstruct
   their exact integer sample matrix by bounded CRT.
2. Build the canonical tau=48*T traces from the repository stress engine.
3. Apply the exact certified degree-10 orbit substitution before extracting
   every reduced source monomial.
4. Derive, for each of the 36 source targets, an exact orbit denominator
   clearer and a rigorous post-restriction coefficient-L1/sample bound.
5. Compute all 36 physical targets together on the same 72 small self-dual
   forms at direct high primes, checkpointing every sample.
6. Add direct primes only until the largest cleared bound has a unique centered
   CRT lift.  Frozen theorem primes are never used in this reconstruction.
7. Solve 36 exact 72x72 systems over QQ and verify every recovered coordinate
   row against all 15 frozen reduced theorem-prime models.

Thus no finite-prime row-space guess is promoted to characteristic zero.
"""
from __future__ import annotations

import json
import random
from collections import defaultdict
from fractions import Fraction
from math import gcd
from pathlib import Path

import numpy as np
import sympy as sp

from .degree12_exact_common import load_theorem_equations, step1c_record, theorem_record
from .exact import crt, rref_q
from .finite_field import _check_prime_field_modulus, matrix_rank
from .forms import dense, layout
from .modular_reconstruction import reduce_fraction
from .polynomial import Poly
from .provenance import atomic_json
from .registry import Registry
from .stress import (
    coefficient_ids,
    matrix_series_product,
    scalar_series_product,
    tau_terms,
)
from .tensors import TensorBudget

REGISTRY_PATH = Path("runs/degree12-generalized-selected-verify/registry.json")
DEG10_OBS = Path("verification/degree10/obstruction_certificate.json")
OUTDIR = Path("verification/all_orders/degree12_exact/interpolation")
CERT_PATH = Path("verification/all_orders/degree12_exact/exact_source_targets.json")
# Compatibility path used by the previous hardened source-row loader.
LEGACY_CERT_PATH = Path("verification/all_orders/degree12_exact/exact_missing_targets.json")
STATUS_PATH = Path("verification/all_orders/degree12_exact/target_recovery_status.json")

FULL_PRIMES = (65521, 65519, 65497, 65479, 65449, 65447, 65437, 65423)
MAX_EXTRA_PHYSICAL_PRIMES = 64
ALGORITHM_VERSION = "degree12-hardened-v3-all-36-post-orbit-direct"
GRAPH12_SAMPLE_BOUND = 10**30
GRAPH_VALUE_BOUND = GRAPH12_SAMPLE_BOUND
B13_SAMPLE_BOUND = 10**30
TAU_TERM_ENTRY_BOUND = 2 * 10**27
TRACE_PAIR_BOUND = 100 * TAU_TERM_ENTRY_BOUND**2
VAR_INDEX = {"A": 0, "B": 1, "C": 2, "D": 3, "U": 4, "V": 5}

EXACT_D12_BUDGET = TensorBudget(
    max_bytes=2 * 1024**3,
    max_multiply_adds=100_000_000_000,
    exact_blas=True,
)


def _lcm(a: int, b: int) -> int:
    return abs(a // gcd(a, b) * b) if a and b else 0


def _ceil_fraction(q: Fraction) -> int:
    q = Fraction(q)
    return (q.numerator + q.denominator - 1) // q.denominator


def _small_selfdual(seed: int, p: int):
    rng = random.Random(int(seed))
    L = layout(10, 5)
    raw = np.zeros(len(L["basis"]), dtype=np.int64)
    for i in L["electric"]:
        raw[int(i)] = 1 if rng.getrandbits(1) else -1
    star = np.zeros_like(raw)
    for i in range(len(raw)):
        star[i] = int(L["star_sign"][i]) * int(raw[int(L["complements"][i])])
    values = raw + star
    if int(np.max(np.abs(values))) > 1:
        raise AssertionError("small self-dual sample is not {-1,0,1}-valued")
    return dense(values % p, p, 10, 5)


def _all_source_descriptors() -> tuple[dict, ...]:
    record = step1c_record()["equation"]
    rows = []
    for meta in record["basis"]:
        rows.append({
            "source_row": int(meta["source_row"]),
            "basis_id": str(meta["basis_id"]),
            "generator": str(meta["generator"]),
            "output_monomial": list(meta.get("output_monomial", ())),
            "output_monomial_text": str(meta["output_monomial_text"]),
            "kind": "basis",
        })
    for meta in record["dependencies"]:
        rows.append({
            "source_row": int(meta["source_row"]),
            "generator": str(meta["generator"]),
            "output_monomial": list(meta.get("output_monomial", ())),
            "output_monomial_text": str(meta["output_monomial_text"]),
            "kind": "dependent",
        })
    rows.sort(key=lambda x: x["source_row"])
    if [x["source_row"] for x in rows] != list(range(36)):
        raise ValueError("corrected Step1C metadata does not label exactly source rows 0..35")
    return tuple(rows)


def _exp_from_meta(meta: dict) -> tuple[int, ...]:
    exponent = [0] * 6
    for item in meta.get("output_monomial", ()):
        name = str(item["variable"])
        if name not in VAR_INDEX:
            raise ValueError(f"unknown reduced variable {name}")
        exponent[VAR_INDEX[name]] += int(item["power"])
    return tuple(exponent)


def _load_full_checkpoints() -> tuple[tuple[int, ...], dict[int, dict]]:
    data: dict[int, dict] = {}
    for p in FULL_PRIMES:
        path = OUTDIR / f"prime{p}.json"
        if not path.exists():
            raise FileNotFoundError(
                f"missing expensive checkpoint {path}; restore it rather than recomputing blindly"
            )
        rec = json.loads(path.read_text())
        if int(rec.get("prime", -1)) != p:
            raise ValueError(f"checkpoint prime mismatch in {path}")
        basis = rec.get("basis", ())
        if len(basis) != 72 or any(len(row) != 72 for row in basis):
            raise ValueError(f"checkpoint {path} does not contain a 72x72 basis matrix")
        data[p] = rec
    seeds = tuple(int(x) for x in data[FULL_PRIMES[0]].get("seeds", ()))
    if len(seeds) != 72 or len(set(seeds)) != 72:
        raise ValueError("full checkpoints do not carry 72 distinct seeds")
    for p in FULL_PRIMES:
        if tuple(int(x) for x in data[p].get("seeds", ())) != seeds:
            raise ValueError(f"interpolation seed ordering changed at full prime {p}")
    return seeds, data


def _centered_integer(residues, primes, bound: int) -> int:
    a, modulus = crt([int(x) for x in residues], list(primes))
    if modulus <= 2 * int(bound):
        raise ValueError(
            f"CRT uniqueness condition failed: modulus={modulus} <= 2*bound={2*int(bound)}"
        )
    z = int(a)
    if z > modulus // 2:
        z -= modulus
    if abs(z) > int(bound):
        raise AssertionError(f"reconstructed integer {z} exceeds certified bound {bound}")
    return z


def _reconstruct_exact_basis_matrix(data: dict[int, dict]) -> list[list[int]]:
    matrix = []
    for i in range(72):
        row = []
        for j in range(72):
            row.append(_centered_integer(
                [data[p]["basis"][i][j] for p in FULL_PRIMES],
                FULL_PRIMES,
                GRAPH12_SAMPLE_BOUND,
            ))
        matrix.append(row)
    for p in FULL_PRIMES:
        got = [[int(x) % p for x in row] for row in matrix]
        want = [[int(x) % p for x in row] for row in data[p]["basis"]]
        if got != want:
            raise AssertionError(f"exact basis matrix does not reduce to saved checkpoint at {p}")
    if matrix_rank(
        [[int(x) % FULL_PRIMES[0] for x in row] for row in matrix],
        FULL_PRIMES[0], ncols=72,
    ) != 72:
        raise AssertionError("saved degree-12 interpolation matrix is not full rank")
    return matrix


def _inverse_q(matrix):
    n = len(matrix)
    augmented = [
        [Fraction(x) for x in row] + [Fraction(int(i == j)) for j in range(n)]
        for i, row in enumerate(matrix)
    ]
    rr, pivots = rref_q(augmented, ncols=2 * n)
    if pivots[:n] != list(range(n)):
        raise ValueError("exact degree-10 orbit leading matrix is singular")
    return [row[n:] for row in rr[:n]]


def _degree10_orbit_substitution_qq(reg: Registry) -> dict[str, Poly]:
    obs = json.loads(DEG10_OBS.read_text())
    ids = tuple(coefficient_ids(reg, 10))
    if ids != tuple(obs.get("coordinate_ids", ())):
        raise ValueError("registry degree-10 coordinate order differs from frozen obstruction certificate")
    rows = [[Fraction(int(x)) for x in row] for row in obs["obstructions_integer_coordinates"]]
    if len(rows) != 12 or any(len(row) != 19 for row in rows):
        raise ValueError("malformed degree-10 obstruction matrix")
    leading = [row[:14] for row in rows]
    tails = [row[14:] for row in rows]

    n = 6
    A, B, C, D, U, V = [Poly.variable(n, i) for i in range(n)]
    zero = Poly(n)
    lower = [A**4, A**2 * B, A * C, A * D, B**2]
    free = {0: U, 12: V}
    dependent_positions = [j for j in range(14) if j not in free]
    inverse = _inverse_q([[row[j] for j in dependent_positions] for row in leading])

    rhs = []
    for i in range(12):
        value = zero - leading[i][0] * U - leading[i][12] * V
        for coefficient, monomial in zip(tails[i], lower):
            value = value - coefficient * monomial
        rhs.append(value)
    dependent = []
    for j in range(12):
        value = zero
        for i in range(12):
            value = value + inverse[j][i] * rhs[i]
        dependent.append(value)

    z10 = [None] * 14
    z10[0] = U
    z10[12] = V
    for pos, value in zip(dependent_positions, dependent):
        z10[pos] = value
    if any(value is None for value in z10):
        raise AssertionError("degree-10 orbit substitution is incomplete")

    substitution: dict[str, Poly] = {
        "I4_1": A,
        "I6_1": B,
        "I6_2": zero,
        "I8_1": C,
        "I8_2": -16 * A**3,
        "I8_3": zero,
        "I8_4": zero,
        "I8_5": zero,
        "I8_6": zero,
        "I4_1^2": D,
    }
    for name, value in zip(ids[10:24], z10):
        substitution[name] = value
    if set(substitution) != set(ids):
        raise AssertionError("degree-10 QQ orbit substitution does not cover all coefficients")
    return substitution


def _degree10_orbit_substitution_mod_p(reg: Registry, p: int) -> dict[str, Poly]:
    return {name: Poly(6, dict(q.terms), p) for name, q in _degree10_orbit_substitution_qq(reg).items()}


def _invariant_value_bound(reg: Registry, name: str, memo=None) -> int:
    memo = {} if memo is None else memo
    if name in memo:
        return memo[name]
    item = reg.items[name]
    if item.graph is not None:
        value = 10 ** (5 * int(item.degree) // 2)
    else:
        value = 1
        for factor in item.factors:
            value *= _invariant_value_bound(reg, factor, memo)
    memo[name] = int(value)
    return int(value)


def _invariant_gradient_bound(reg: Registry, name: str, value_memo=None, grad_memo=None) -> int:
    value_memo = {} if value_memo is None else value_memo
    grad_memo = {} if grad_memo is None else grad_memo
    if name in grad_memo:
        return grad_memo[name]
    item = reg.items[name]
    if item.graph is not None:
        value = int(item.degree) * 10 ** (5 * int(item.degree) // 2 - 5)
    else:
        factors = tuple(item.factors)
        value = 0
        for j, factor in enumerate(factors):
            term = _invariant_gradient_bound(reg, factor, value_memo, grad_memo)
            for k, other in enumerate(factors):
                if k != j:
                    term *= _invariant_value_bound(reg, other, value_memo)
            value += term
    grad_memo[name] = int(value)
    return int(value)


def _derived_tau_entry_bound(reg: Registry) -> dict:
    ids = tuple(coefficient_ids(reg, 10))
    value_memo, grad_memo = {}, {}
    free_bound = 10**4
    identity_bound = max(
        (int(reg.items[name].degree) - 2) * _invariant_value_bound(reg, name, value_memo)
        for name in ids
    )
    gradient_names = [name for name in ids if reg.items[name].degree + 2 <= 10]
    gradient_pair_bound = 0
    worst_pair = None
    for i, left in enumerate(gradient_names):
        for right in gradient_names[i:]:
            if reg.items[left].degree + reg.items[right].degree - 2 > 10:
                continue
            gl = _invariant_gradient_bound(reg, left, value_memo, grad_memo)
            gr = _invariant_gradient_bound(reg, right, value_memo, grad_memo)
            symmetry = 1 if left == right else 2
            bound = 25 * symmetry * 10**4 * gl * gr
            if bound > gradient_pair_bound:
                gradient_pair_bound, worst_pair = int(bound), (left, right)
    derived = max(int(free_bound), int(identity_bound), int(gradient_pair_bound))
    if derived > TAU_TERM_ENTRY_BOUND:
        raise AssertionError(
            f"derived tau-entry bound {derived} exceeds certified envelope {TAU_TERM_ENTRY_BOUND}"
        )
    return {
        "free_bound": int(free_bound),
        "identity_bound": int(identity_bound),
        "gradient_pair_bound": int(gradient_pair_bound),
        "worst_gradient_pair": list(worst_pair) if worst_pair else None,
        "derived_max": int(derived),
        "certified_envelope": int(TAU_TERM_ENTRY_BOUND),
    }


def _tau_bound_series(reg: Registry) -> dict[tuple[int, tuple[str, ...]], int]:
    out: dict[tuple[int, tuple[str, ...]], int] = defaultdict(int)
    ids = tuple(coefficient_ids(reg, 10))
    value_memo, grad_memo = {}, {}
    out[(2, ())] += 10**4
    for name in ids:
        d = int(reg.items[name].degree)
        out[(d, (name,))] += (d - 2) * _invariant_value_bound(reg, name, value_memo)
    gradient_names = [name for name in ids if reg.items[name].degree + 2 <= 10]
    for i, left in enumerate(gradient_names):
        for right in gradient_names[i:]:
            d = int(reg.items[left].degree) + int(reg.items[right].degree) - 2
            if d > 10:
                continue
            symmetry = 1 if left == right else 2
            bound = (
                25 * symmetry * 10**4
                * _invariant_gradient_bound(reg, left, value_memo, grad_memo)
                * _invariant_gradient_bound(reg, right, value_memo, grad_memo)
            )
            out[(d, tuple(sorted((left, right))))] += int(bound)
    return dict(out)


def _bound_convolve(left, right, *, max_degree=12, matrix=False):
    out = defaultdict(int)
    for (d1, m1), b1 in left.items():
        for (d2, m2), b2 in right.items():
            d = int(d1) + int(d2)
            if d > max_degree:
                continue
            factor = 10 if matrix else 1
            out[(d, tuple(sorted(m1 + m2)))] += factor * int(b1) * int(b2)
    return dict(out)


def _trace_bound_series(reg: Registry) -> dict[str, dict]:
    traces: dict[str, dict] = {}
    value_memo = {}
    t1 = defaultdict(int)
    for name in coefficient_ids(reg, 10):
        d = int(reg.items[name].degree)
        t1[(d, (name,))] += 10 * (d - 2) * _invariant_value_bound(reg, name, value_memo)
    traces["tr1"] = dict(t1)

    tau = _tau_bound_series(reg)
    current = {(0, ()): 1}
    # First matrix factor creates an entry without a summed matrix index.
    for k in range(1, 7):
        nxt = defaultdict(int)
        for (d1, m1), b1 in current.items():
            for (d2, m2), b2 in tau.items():
                d = d1 + d2
                if d > 12:
                    continue
                factor = 1 if k == 1 else 10
                nxt[(d, tuple(sorted(m1 + m2)))] += factor * int(b1) * int(b2)
        current = dict(nxt)
        if k >= 2:
            traces[f"tr{k}"] = {key: 10 * int(value) for key, value in current.items()}
    return traces


def _generator_bound_series(generator: str, traces: dict[str, dict]) -> dict:
    series = {(0, ()): 1}
    for token in generator.split("*"):
        if token not in traces:
            raise ValueError(f"unsupported stress factor in source target: {token}")
        series = _bound_convolve(series, traces[token], max_degree=12, matrix=False)
    return series


def _target_bound(reg: Registry, meta: dict, *, traces=None, substitution=None) -> tuple[int, int, dict]:
    theorem = theorem_record()
    if "integral" not in str(theorem.get("formal_integrality_basis", "")).lower():
        raise ValueError("frozen theorem lacks the full-map formal-integrality input")
    traces = _trace_bound_series(reg) if traces is None else traces
    substitution = _degree10_orbit_substitution_qq(reg) if substitution is None else substitution
    series = _generator_bound_series(str(meta["generator"]), traces)
    exponent = _exp_from_meta(meta)
    denominator = 1
    weighted = Fraction(0)
    contributing = 0
    for (degree, monomial), raw_bound in series.items():
        if int(degree) != 12:
            continue
        q = Poly.constant(6, 1)
        for name in monomial:
            if name not in substitution:
                raise ValueError(f"bound series references coefficient outside degree-10 orbit: {name}")
            q = q * substitution[name]
        coefficient = Fraction(q.terms.get(exponent, 0))
        if not coefficient:
            continue
        contributing += 1
        denominator = _lcm(denominator, coefficient.denominator)
        weighted += abs(coefficient) * int(raw_bound)
    if contributing == 0 or weighted <= 0:
        raise AssertionError(
            f"source row {meta['source_row']} {meta['generator']}|{meta['output_monomial_text']} "
            "has no nonzero post-orbit structural contribution"
        )
    cleared = _ceil_fraction(weighted * denominator)
    return int(denominator), int(cleared), {
        "source_row": int(meta["source_row"]),
        "generator": str(meta["generator"]),
        "output_monomial": str(meta["output_monomial_text"]),
        "orbit_denominator_clearer": str(denominator),
        "post_orbit_l1_bound": str(_ceil_fraction(weighted)),
        "cleared_integer_bound": str(cleared),
        "contributing_full_monomials": int(contributing),
        "full_map_integrality_input": theorem["formal_integrality_basis"],
    }


def _b11_characteristic_zero_bound(reg: Registry) -> tuple[int, int, dict]:
    """Compatibility wrapper; now uses the generic post-orbit bound engine."""
    meta = next(x for x in _all_source_descriptors() if x["source_row"] == 21)
    denominator, cleared, proof = _target_bound(reg, meta)
    proof = dict(proof)
    proof["B11_cleared_integer_bound"] = str(cleared)
    return denominator, cleared, proof


def _trace_value_series(reg: Registry, form, p: int) -> dict[str, dict]:
    cache: dict = {}
    traces: dict[str, dict] = {}
    t1 = {}
    for name in coefficient_ids(reg, 10):
        d = int(reg.items[name].degree)
        value = reg.evaluate(name, form, p, cache=cache, budget=EXACT_D12_BUDGET)
        coefficient = 10 * (d - 2) * int(value) % p
        if coefficient:
            t1[(d, (name,))] = coefficient
    traces["tr1"] = t1

    tau = tau_terms(reg, form, p, 10, cache=cache, budget=EXACT_D12_BUDGET)
    current = {(0, ()): np.eye(10, dtype=np.int64)}
    for k in range(1, 7):
        current = matrix_series_product(current, tau, p, 12, EXACT_D12_BUDGET)
        if k >= 2:
            traces[f"tr{k}"] = {
                key: int(np.trace(matrix)) % p
                for key, matrix in current.items()
                if int(np.trace(matrix)) % p
            }
    return traces


def _generator_value_series(generator: str, traces: dict[str, dict], p: int) -> dict:
    series = {(0, ()): 1}
    for token in generator.split("*"):
        if token not in traces:
            raise ValueError(f"unsupported stress factor in source target: {token}")
        series = scalar_series_product(series, traces[token], p, 12)
    return series


def _restricted_target_value(series: dict, substitution: dict[str, Poly], exponent: tuple[int, ...], p: int) -> int:
    reduced = Poly(6, p=p)
    for (degree, monomial), value in series.items():
        if int(degree) != 12:
            continue
        term = Poly.constant(6, int(value), p)
        for name in monomial:
            if name not in substitution:
                raise ValueError(f"physical source term references coefficient outside certified orbit: {name}")
            term = term * substitution[name]
        reduced = reduced + term
    return int(reduced.terms.get(exponent, 0)) % p


def _all_target_values_for_sample(reg: Registry, form, p: int, substitution: dict[str, Poly], descriptors) -> dict[int, int]:
    traces = _trace_value_series(reg, form, p)
    generator_cache: dict[str, dict] = {}
    out = {}
    for meta in descriptors:
        generator = str(meta["generator"])
        if generator not in generator_cache:
            generator_cache[generator] = _generator_value_series(generator, traces, p)
        out[int(meta["source_row"])] = _restricted_target_value(
            generator_cache[generator], substitution, _exp_from_meta(meta), p
        )
    return out


def _all_target_checkpoint(reg: Registry, seeds: tuple[int, ...], p: int, descriptors, full_data: dict[int, dict]) -> dict:
    path = OUTDIR / f"all_source_targets_prime{p}.json"
    rows = {str(i): [] for i in range(36)}
    if path.exists():
        rec = json.loads(path.read_text())
        compatible = (
            int(rec.get("prime", -1)) == p
            and tuple(int(x) for x in rec.get("seeds", ())) == seeds
            and rec.get("algorithm_version") == ALGORITHM_VERSION
        )
        if compatible:
            saved = rec.get("values_by_source_row", {})
            lengths = {len(saved.get(str(i), ())) for i in range(36)}
            if len(lengths) != 1:
                raise ValueError(f"source-target checkpoint {path} has uneven row lengths")
            done = next(iter(lengths), 0)
            if done > 72:
                raise ValueError(f"source-target checkpoint {path} has too many samples")
            rows = {str(i): [int(x) % p for x in saved.get(str(i), ())] for i in range(36)}
            if done == 72 and rec.get("status") == "complete":
                return rec
            print(f"resume all-source physical prime {p}: {done}/72 cached", flush=True)
        else:
            archive = path.with_name(path.stem + ".superseded.json")
            if not archive.exists():
                atomic_json(archive, rec)
    done = len(rows["0"])
    substitution = _degree10_orbit_substitution_mod_p(reg, p)
    base = {
        "schema": 1,
        "prime": p,
        "seeds": list(seeds),
        "algorithm_version": ALGORITHM_VERSION,
        "definition": "all_36_degree12_scalar_sources_after_exact_degree10_orbit_restriction",
    }
    for index in range(done, 72):
        values = _all_target_values_for_sample(
            reg, _small_selfdual(seeds[index], p), p, substitution, descriptors
        )
        for row in range(36):
            rows[str(row)].append(int(values[row]) % p)
        atomic_json(path, {
            **base,
            "status": "running",
            "samples_finished": index + 1,
            "values_by_source_row": rows,
        })
        print(f"all-source physical prime {p}: {index + 1}/72", flush=True)
    rec = {
        **base,
        "status": "complete",
        "samples_finished": 72,
        "values_by_source_row": rows,
    }
    # Cross-check B11 against the previously certified specialized checkpoint when available.
    old_b11 = OUTDIR / f"restricted_B11_prime{p}.json"
    if old_b11.exists():
        previous = json.loads(old_b11.read_text())
        if previous.get("status") == "complete" and len(previous.get("values", ())) == 72:
            if [int(x) % p for x in previous["values"]] != rows["21"]:
                raise AssertionError(f"generic all-source B11 disagrees with specialized B11 checkpoint at {p}")
    # Cross-check B13 against the original full interpolation target when available.
    if p in full_data and len(full_data[p].get("B13", ())) == 72:
        if [int(x) % p for x in full_data[p]["B13"]] != rows["35"]:
            raise AssertionError(f"generic all-source B13 disagrees with original B13 checkpoint at {p}")
    atomic_json(path, rec)
    return rec


def _extra_prime_stream(excluded: set[int]):
    candidate = 65419
    while candidate > 30000:
        if candidate not in excluded:
            try:
                _check_prime_field_modulus(candidate)
            except ValueError:
                pass
            else:
                yield candidate
        candidate -= 2


def _exact_solve_q(matrix: list[list[int]], values: list[Fraction]) -> list[Fraction]:
    M = sp.Matrix(matrix)
    y = sp.Matrix([sp.Rational(q.numerator, q.denominator) for q in values])
    x = M.inv(method="DM") * y
    result = [Fraction(int(sp.Rational(q).p), int(sp.Rational(q).q)) for q in x]
    sx = sp.Matrix([sp.Rational(q.numerator, q.denominator) for q in result])
    if M * sx != y:
        raise AssertionError("exact 72x72 QQ solve does not reproduce source target samples")
    return result


def _verify_row(name: str, row: list[Fraction], source_row: int, checkpoints: dict[int, dict], exact_basis, theorem_primes, equations):
    if len(row) != 72:
        raise ValueError(f"{name} exact target row is not 72-dimensional")
    for p, checkpoint in checkpoints.items():
        coords = [reduce_fraction(q, p) for q in row]
        matrix = [[int(x) % p for x in sample] for sample in exact_basis]
        values = checkpoint["values_by_source_row"][str(source_row)]
        for i, sample in enumerate(matrix):
            got = sum(int(a) * int(x) for a, x in zip(sample, coords)) % p
            if got != int(values[i]) % p:
                raise AssertionError(f"{name} exact row fails direct physical sample {p}/{i}")
    for p in theorem_primes:
        got = [reduce_fraction(q, p) for q in row]
        want = [int(x) % p for x in equations[p][source_row][:72]]
        if got != want:
            raise AssertionError(f"{name} exact row disagrees with frozen theorem model at {p}")


def recover_exact_targets(output: str | Path = CERT_PATH, status_output: str | Path = STATUS_PATH) -> dict:
    descriptors = _all_source_descriptors()
    model_dir, theorem_primes, equations, _ = load_theorem_equations()
    seeds, full_data = _load_full_checkpoints()
    exact_basis = _reconstruct_exact_basis_matrix(full_data)
    reg = Registry.load(REGISTRY_PATH)
    if len(reg.degree_bases.get(12, ())) != 72:
        raise ValueError("frozen registry lacks the complete 72-dimensional degree-12 basis")

    print("PASS: eight expensive full checkpoints loaded and exact basis matrix reconstructed", flush=True)
    print(f"theorem model directory: {model_dir}", flush=True)

    tau_audit = _derived_tau_entry_bound(reg)
    trace_bounds = _trace_bound_series(reg)
    substitution_q = _degree10_orbit_substitution_qq(reg)
    bound_records: dict[int, dict] = {}
    denominators: dict[int, int] = {}
    cleared_bounds: dict[int, int] = {}
    for meta in descriptors:
        d, b, proof = _target_bound(reg, meta, traces=trace_bounds, substitution=substitution_q)
        row = int(meta["source_row"])
        denominators[row], cleared_bounds[row], bound_records[row] = d, b, proof
    max_row = max(cleared_bounds, key=cleared_bounds.get)
    max_bound = cleared_bounds[max_row]
    print(
        f"post-orbit bounds certified for all 36 source targets; worst row={max_row} "
        f"cleared_bound_digits={len(str(max_bound))}",
        flush=True,
    )
    for row in sorted(cleared_bounds, key=cleared_bounds.get, reverse=True)[:6]:
        print(
            f"  row {row}: D={denominators[row]} bound_digits={len(str(cleared_bounds[row]))} "
            f"{bound_records[row]['generator']}|{bound_records[row]['output_monomial']}",
            flush=True,
        )

    # Use only direct physical primes for the characteristic-zero sample lift.
    excluded = set(theorem_primes)
    candidate_primes = iter(FULL_PRIMES)
    extra_stream = _extra_prime_stream(set(FULL_PRIMES) | excluded)
    common_denominator = 1
    for d in denominators.values():
        common_denominator = _lcm(common_denominator, d)
    used_primes: list[int] = []
    extra_used: list[int] = []
    checkpoints: dict[int, dict] = {}
    modulus = 1
    extras_attempted = 0
    while modulus <= 2 * max_bound or len(used_primes) < 2:
        try:
            p = next(candidate_primes)
            is_extra = False
        except StopIteration:
            if extras_attempted >= MAX_EXTRA_PHYSICAL_PRIMES:
                atomic_json(Path(status_output), {
                    "schema": 4,
                    "status": "blocked_all_source_target_height_modulus",
                    "worst_source_row": max_row,
                    "worst_cleared_bound": str(max_bound),
                    "direct_physical_primes": used_primes,
                    "modulus": str(modulus),
                    "instruction": "add direct physical target primes; do not weaken the bound gate",
                })
                raise RuntimeError("direct source-target prime product did not close the certified height bound")
            p = next(extra_stream)
            extras_attempted += 1
            is_extra = True
        if common_denominator % p == 0:
            print(f"SKIP all-source prime {p}: divides an orbit denominator", flush=True)
            continue
        try:
            checkpoint = _all_target_checkpoint(reg, seeds, p, descriptors, full_data)
        except ValueError as exc:
            if "prime divides coefficient denominator" in str(exc):
                print(f"SKIP all-source prime {p}: {exc}", flush=True)
                continue
            raise
        checkpoints[p] = checkpoint
        used_primes.append(p)
        modulus *= p
        if is_extra:
            extra_used.append(p)
        print(
            f"all-source bounded-CRT prime {p} accepted; usable_primes={len(used_primes)} "
            f"modulus_digits={len(str(modulus))}",
            flush=True,
        )

    recovered = []
    by_source_row = {}
    for meta in descriptors:
        row_index = int(meta["source_row"])
        denominator = denominators[row_index]
        bound = cleared_bounds[row_index]
        target_values: list[Fraction] = []
        for i in range(72):
            residues = [
                int(checkpoints[p]["values_by_source_row"][str(row_index)][i]) * (denominator % p) % p
                for p in used_primes
            ]
            numerator = _centered_integer(residues, used_primes, bound)
            target_values.append(Fraction(numerator, denominator))
        coords = _exact_solve_q(exact_basis, target_values)
        label = str(meta.get("basis_id") or f"R{row_index}")
        _verify_row(label, coords, row_index, checkpoints, exact_basis, theorem_primes, equations)
        record = {
            **meta,
            "coordinate_domain": "QQ",
            "coordinates": [str(x) for x in coords],
            "orbit_denominator_clearer": str(denominator),
            "cleared_sample_bound": str(bound),
            "bound_proof": bound_records[row_index],
        }
        recovered.append(record)
        by_source_row[str(row_index)] = record
        print(f"PASS: exact physical source target row {row_index}/35", flush=True)

    payload = {
        "schema": 5,
        "status": "exact_characteristic_zero_all_degree12_source_targets",
        "proof_method": (
            "all_36_direct_physical_post_orbit_targets + rigorous denominator-cleared bounded CRT "
            "+ exact 72x72 QQ interpolation + disjoint 15-prime reduced-model verification"
        ),
        "algorithm_version": ALGORITHM_VERSION,
        "full_checkpoint_primes": list(FULL_PRIMES),
        "direct_physical_primes": used_primes,
        "extra_direct_physical_primes": extra_used,
        "theorem_primes_used_in_reconstruction": [],
        "theorem_primes_verification_only": list(theorem_primes),
        "exact_basis_sample_bound": str(GRAPH12_SAMPLE_BOUND),
        "tau_entry_bound_audit": tau_audit,
        "worst_source_row": int(max_row),
        "worst_cleared_sample_bound": str(max_bound),
        "rows": recovered,
        "rows_by_source_row": by_source_row,
        "all_36_direct_physical_crosschecks": True,
        "all_36_theorem_prime_crosschecks": True,
    }
    # Compatibility named entries for existing tools/tests.
    for record in recovered:
        bid = record.get("basis_id")
        if bid:
            payload[str(bid)] = record
    out = Path(output)
    out.parent.mkdir(parents=True, exist_ok=True)
    atomic_json(out, payload)
    atomic_json(LEGACY_CERT_PATH, payload)
    atomic_json(Path(status_output), {
        "schema": 4,
        "status": "all_36_exact_source_targets_verified",
        "certificate": str(out),
    })
    return payload
