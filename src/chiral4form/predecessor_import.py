"""Import validated predecessor degree-12 modular coefficient certificates.

These certificates are an additional source of modular images.  Importing them
does not promote any characteristic-zero claim: the existing degree12-certify
pipeline remains the only gate for QQ candidate status.
"""
from __future__ import annotations

import json
import shutil
from pathlib import Path
from typing import Mapping

from .finite_field import solve_many
from .fitting import fields_from_json, fields_to_json
from .polynomial import Poly
from .registry import import_predecessor_registry
from .stress import coefficient_ids
from .weighted import stress_generators

GENERATOR_MAP = {
    "tr_tau": "tr1",
    "tr_tau2": "tr2",
    "tr_tau3": "tr3",
    "tr_tau4": "tr4",
    "tr_tau^2": "tr1*tr1",
    "tr_tau*tr_tau2": "tr1*tr2",
    "tr_tau2^2": "tr2*tr2",
    "tr_tau5": "tr5",
    "tr_tau*tr_tau3": "tr1*tr3",
    "tr_tau2*tr_tau3": "tr2*tr3",
    "tr_tau6": "tr6",
    "tr_tau*tr_tau4": "tr1*tr4",
    "tr_tau2*tr_tau4": "tr2*tr4",
    "tr_tau3^2": "tr3*tr3",
    "tr_tau^3": "tr1*tr1*tr1",
    "tr_tau^2*tr_tau2": "tr1*tr1*tr2",
    "tr_tau*tr_tau2^2": "tr1*tr2*tr2",
    "tr_tau2^3": "tr2*tr2*tr2",
}

EXPECTED_PREDECESSOR_PRIMES = (32693, 32713, 32717, 32719, 32749, 32771)


def _catalogue():
    return [
        {
            "id": g.id,
            "factors": list(g.factors),
            "leading_degree": g.leading_degree,
        }
        for g in stress_generators(12)
    ]


def full_model_from_certificate(certificate: Mapping, source_root: Path) -> dict:
    p = int(certificate["prime"])
    if not certificate.get("all_holdouts_passed"):
        raise ValueError(f"predecessor certificate {p} did not pass all holdouts")
    normalization = str(certificate.get("normalization", ""))
    if "tau=48*T" not in normalization.replace(" ", ""):
        raise ValueError(f"predecessor certificate {p} has unexpected normalization")

    registry = import_predecessor_registry(source_root, 12)
    ids = tuple(coefficient_ids(registry, 12))
    if len(ids) != 96:
        raise ValueError(f"expected 96 cumulative coordinates, got {len(ids)}")
    index = {name: i for i, name in enumerate(ids)}
    current_generators = {g.id for g in stress_generators(12)}
    if set(GENERATOR_MAP.values()) != current_generators:
        raise AssertionError("predecessor/current generator map is incomplete")

    raw = {g: [dict() for _ in ids] for g in sorted(current_generators)}
    seen_generators = set()

    for target in certificate["targets"]:
        if not target.get("holdout_passed", False):
            raise ValueError(f"predecessor target failed holdout: {target.get('id')}")
        old_generator = target["generator"]
        if old_generator not in GENERATOR_MAP:
            raise ValueError(f"unknown predecessor generator: {old_generator}")
        generator = GENERATOR_MAP[old_generator]
        seen_generators.add(generator)

        degree = int(target["field_degree"])
        expected_basis = list(registry.degree_bases[degree])
        basis = list(target["basis"])
        if basis != expected_basis:
            raise ValueError(
                f"basis ordering differs for degree {degree} in target {target['id']}"
            )
        coordinates = list(target["coordinates"])
        if len(coordinates) != len(basis):
            raise ValueError(f"coordinate length mismatch in target {target['id']}")

        powers = [0] * len(ids)
        for name in target["coefficient_monomial"]:
            if name not in index:
                raise ValueError(f"unknown coefficient variable {name}")
            powers[index[name]] += 1
        monomial = tuple(powers)

        for name, value in zip(basis, coordinates):
            j = index[name]
            raw[generator][j][monomial] = (
                raw[generator][j].get(monomial, 0) + int(value)
            ) % p

    if seen_generators != current_generators:
        missing = sorted(current_generators - seen_generators)
        raise ValueError(f"predecessor certificate lacks generators: {missing}")

    fields = {
        g: tuple(Poly(len(ids), terms, p) for terms in coordinate_terms)
        for g, coordinate_terms in raw.items()
    }
    model = fields_to_json(
        ids,
        fields,
        "validated_predecessor_degree12_modular_coefficient_certificate",
    )
    model["generator_catalogue"] = _catalogue()
    model["source_prime"] = p
    model["source_kind"] = "pinned_predecessor_interacting_degree12_certificate"
    model["source_certificate_engine_sha256"] = certificate.get("engine_sha256")
    model["source_degree12_sha256"] = certificate.get("degree12_sha256")
    model["source_all_holdouts_passed"] = True
    return model


def restrict_to_certified_degree10_orbit(model: Mapping, obstruction: Mapping) -> dict:
    ids96, fields96 = fields_from_json(model)
    if len(ids96) != 96:
        raise ValueError("expected a 96-coordinate full degree-12 model")

    expected_prefix = tuple(obstruction["coordinate_ids"])
    if tuple(ids96[:24]) != expected_prefix:
        raise ValueError("degree-10 coordinate prefix differs from frozen certificate")

    sample_poly = next(iter(next(iter(fields96.values()))))
    p = int(sample_poly.p)
    n = 78
    A, B, C, D, U, V = [Poly.variable(n, i, p) for i in range(6)]
    Z12 = [Poly.variable(n, 6 + i, p) for i in range(72)]
    zero = Poly(n, p=p)

    rows = [[int(x) % p for x in row] for row in obstruction["obstructions_integer_coordinates"]]
    if len(rows) != 12 or any(len(row) != 19 for row in rows):
        raise ValueError("malformed degree-10 obstruction certificate")
    leading = [row[:14] for row in rows]
    tails = [row[14:] for row in rows]

    lower = [A**4, A**2 * B, A * C, A * D, B**2]
    free_positions = {0: U, 12: V}
    dependent_positions = [j for j in range(14) if j not in free_positions]

    Adep = [[row[j] for j in dependent_positions] for row in leading]
    identity = [[int(i == j) for j in range(12)] for i in range(12)]
    inverse = solve_many(Adep, identity, p)

    rhs = []
    for i in range(12):
        value = zero - leading[i][0] * U - leading[i][12] * V
        for coeff, monomial in zip(tails[i], lower):
            value = value - coeff * monomial
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

    substitution = {
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
    for name, value in zip(ids96[10:24], z10):
        substitution[name] = value
    for name, value in zip(ids96[24:], Z12):
        substitution[name] = value
    values = [substitution[name] for name in ids96]

    reduced_fields = {}
    for gname, field in fields96.items():
        full = [component.substitute(values, n_out=n) for component in field]
        reduced_fields[gname] = tuple(
            [full[0], full[1], full[3], full[9], full[10], full[22]]
            + full[24:96]
        )

    reduced_ids = ("A", "B", "C", "D", "U", "V", *ids96[24:])
    out = fields_to_json(
        reduced_ids,
        reduced_fields,
        "predecessor_certificate_restricted_to_certified_degree10_orbit",
    )
    out["generator_catalogue"] = model.get("generator_catalogue", _catalogue())
    out["source_prime"] = p
    out["source_kind"] = model.get("source_kind")
    out["restriction"] = "frozen_global_six_parameter_degree10_orbit"
    return out


def import_predecessor_certificates(
    source_root: Path,
    existing_models: Path,
    output_dir: Path,
    obstruction_path: Path,
) -> dict:
    source_root = Path(source_root)
    existing_models = Path(existing_models)
    output_dir = Path(output_dir)
    obstruction_path = Path(obstruction_path)

    cert_dir = source_root / "results" / "stress_flow" / "certificates"
    if not cert_dir.is_dir():
        raise ValueError(f"certificate directory missing: {cert_dir}")
    obstruction = json.loads(obstruction_path.read_text())

    output_dir.mkdir(parents=True, exist_ok=True)
    existing_primes = []
    for path in sorted(existing_models.glob("fields_prime*.json")):
        target = output_dir / path.name
        shutil.copy2(path, target)
        existing_primes.append(int(path.stem.replace("fields_prime", "")))

    imported = []
    for p in EXPECTED_PREDECESSOR_PRIMES:
        cert_path = cert_dir / f"interacting_degree12_{p}.json"
        if not cert_path.exists():
            raise ValueError(f"missing predecessor certificate: {cert_path}")
        certificate = json.loads(cert_path.read_text())
        if int(certificate["prime"]) != p:
            raise ValueError(f"prime mismatch in {cert_path}")
        full = full_model_from_certificate(certificate, source_root)
        reduced = restrict_to_certified_degree10_orbit(full, obstruction)
        target = output_dir / f"fields_prime{p}.json"
        if target.exists():
            raise ValueError(f"prime collision at {p}: {target}")
        target.write_text(json.dumps(reduced, indent=2) + "\n")
        imported.append(p)

    return {
        "schema": 1,
        "status": "predecessor_degree12_certificates_imported",
        "existing_primes": sorted(existing_primes),
        "imported_primes": imported,
        "total_primes": len(existing_primes) + len(imported),
        "output_dir": str(output_dir),
    }
