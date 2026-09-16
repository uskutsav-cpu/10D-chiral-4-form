"""Shared fail-closed helpers for the exact degree-12 obstruction theorem.

This module deliberately centralizes every structural assumption that used to
be duplicated across experimental scripts: theorem primes, reduced-coordinate
ordering, the 13 independent source-row labels, and model discovery.
"""
from __future__ import annotations

import json
import re
from pathlib import Path

from .degree12_certificate import obstruction_equations
from .normalization_audit import load_records

THEOREM_PATH = Path("verification/degree12/qq/degree12_local_theorem_certificate.json")
STEP1C_PATH = Path("verification/degree12/qq/step1c_source_identity_targets.json")

EXPECTED_SOURCE_ROWS = (0, 1, 2, 3, 4, 9, 16, 17, 19, 20, 21, 28, 35)
EXPECTED_B11 = (21, "tr2", "A^4")
EXPECTED_B13 = (35, "tr6", "1")


def theorem_record() -> dict:
    rec = json.loads(THEOREM_PATH.read_text())
    required = {
        "status": "characteristic_zero_degree12_local_orbit_theorem_certificate",
        "qq_theorem": True,
        "obstruction_nullity": 68,
        "degree12_leading_obstruction_rank": 68,
        "reduced_ambient_dimension": 78,
        "local_orbit_dimension": 10,
    }
    for key, expected in required.items():
        if rec.get(key) != expected:
            raise ValueError(
                f"frozen degree-12 theorem gate failed: {key}={rec.get(key)!r}, "
                f"expected {expected!r}"
            )
    eq = rec.get("equation_rank", {})
    if eq.get("lower_bound") != 13 or eq.get("upper_bound") != 13:
        raise ValueError("frozen degree-12 equation rank is not exactly 13")
    primes = tuple(int(x) for x in rec.get("primes", ()))
    if len(primes) != 15 or len(set(primes)) != 15:
        raise ValueError("frozen degree-12 theorem does not carry 15 distinct primes")
    return rec


def step1c_record() -> dict:
    rec = json.loads(STEP1C_PATH.read_text())
    if rec.get("status") != "step1c_source_identity_targets_corrected":
        raise ValueError("corrected step1c source-row metadata is missing")
    basis = rec.get("equation", {}).get("basis", ())
    rows = tuple(int(x["source_row"]) for x in basis)
    if rows != EXPECTED_SOURCE_ROWS:
        raise ValueError(f"unexpected independent source rows: {rows}")
    by_id = {x["basis_id"]: x for x in basis}
    for bid, expected in (("B11", EXPECTED_B11), ("B13", EXPECTED_B13)):
        row = by_id.get(bid)
        if row is None:
            raise ValueError(f"step1c metadata lacks {bid}")
        got = (
            int(row["source_row"]),
            str(row["generator"]),
            str(row["output_monomial_text"]),
        )
        if got != expected:
            raise ValueError(f"{bid} metadata changed: {got!r} != {expected!r}")
    return rec


def _prime_from_name(path: Path) -> int | None:
    match = re.fullmatch(r"fields_prime(\d+)\.json", path.name)
    return int(match.group(1)) if match else None


def discover_model_directory() -> tuple[Path, tuple[int, ...]]:
    """Find the one run directory containing every frozen theorem prime."""
    required = tuple(int(x) for x in theorem_record()["primes"])
    need = set(required)
    groups: dict[Path, set[int]] = {}
    for path in Path("runs").rglob("fields_prime*.json"):
        p = _prime_from_name(path)
        if p is not None:
            groups.setdefault(path.parent, set()).add(p)
    exact = [parent for parent, primes in groups.items() if need <= primes]
    if not exact:
        partial = sorted(
            (
                (len(need & primes), str(parent), sorted(need - primes))
                for parent, primes in groups.items()
            ),
            reverse=True,
        )[:10]
        raise FileNotFoundError(
            "no run directory contains all frozen degree-12 theorem primes; "
            f"best candidates={partial}"
        )
    # Deterministic preference: shortest path, then lexical order.
    parent = sorted(exact, key=lambda p: (len(str(p)), str(p)))[0]
    return parent, required


def load_theorem_equations() -> tuple[Path, tuple[int, ...], dict[int, list[list[int]]], list[str]]:
    model_dir, primes = discover_model_directory()
    records = load_records([model_dir / f"fields_prime{p}.json" for p in primes])
    equations: dict[int, list[list[int]]] = {}
    labels_ref = None
    for p in primes:
        ids, _, labels, rows, _ = obstruction_equations(records[p], p)
        if len(ids) != 78 or tuple(ids[:6]) != ("A", "B", "C", "D", "U", "V"):
            raise ValueError(f"theorem model at prime {p} has wrong reduced coordinates")
        if len(rows) != 36:
            raise ValueError(f"theorem model at prime {p} has {len(rows)} rows, expected 36")
        if labels_ref is None:
            labels_ref = list(labels)
        elif list(labels) != labels_ref:
            raise ValueError(f"obstruction ansatz labels changed at theorem prime {p}")
        equations[p] = rows
    step1c_record()  # structural row-label gate
    return model_dir, primes, equations, list(labels_ref or ())
