"""Import validated predecessor degree-12 modular coefficient certificates.

The adapter is deliberately strict:
- predecessor certificates must have passed their own holdouts;
- their lower-degree vector fields must reduce exactly to the frozen degree-10
  characteristic-zero model modulo the certificate prime;
- reduced-model structural metadata is copied from a native reduced model
  rather than invented by the importer.

Importing these modular images never promotes a characteristic-zero claim.
"""
from __future__ import annotations

import copy
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


def _support_signature(record: Mapping) -> dict:
    """Coefficient-free polynomial support for structural comparison."""
    out = {}
    for generator, components in record["fields"].items():
        out[generator] = [
            tuple(sorted(tuple(term["powers"]) for term in component.get("terms", ())))
            for component in components
        ]
    return out


def _align_reduced_metadata(reduced: dict, native_reference: Mapping) -> dict:
    """Mirror optional theorem-neutral metadata from a native reduced model."""
    out = dict(reduced)

    # generator_catalogue is optional in the native reduced schema.  If native
    # models omit it, imported models must omit it as well.  If they carry it,
    # use the exact native metadata rather than a separately generated copy.
    if "generator_catalogue" in native_reference:
        out["generator_catalogue"] = copy.deepcopy(
            native_reference["generator_catalogue"]
        )
    else:
        out.pop("generator_catalogue", None)

    if "restriction" in native_reference:
        out["restriction"] = native_reference["restriction"]
    else:
        out.pop("restriction", None)

    out["metadata_alignment"] = {
        "status": "aligned_to_native_reduced_model",
        "reference_prime": int(native_reference["source_prime"]),
        "generator_catalogue_present": "generator_catalogue" in native_reference,
        "restriction": native_reference.get("restriction"),
    }
    return out


def _assert_same_reduced_structure(imported: Mapping, native_reference: Mapping) -> None:
    if imported["coordinate_ids"] != native_reference["coordinate_ids"]:
        raise ValueError("imported/native reduced coordinate ordering differs")
    if set(imported["fields"]) != set(native_reference["fields"]):
        raise ValueError("imported/native reduced generator sets differ")
    # Support equality is a strong check here because the nine native models
    # already exhibited the same support and the predecessor certificates are
    # being imported as additional images of that same map.
    if _support_signature(imported) != _support_signature(native_reference):
        raise ValueError("imported/native reduced polynomial support differs")


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


def verify_lower_sector_against_degree10(
    full_model: Mapping,
    frozen_degree10_model: Mapping,
) -> dict:
    """Prove imported d<=10 maps equal the frozen QQ model modulo this prime."""
    ids96, fields96 = fields_from_json(full_model)
    ids10, fields10 = fields_from_json(frozen_degree10_model)
    if tuple(ids96[:24]) != tuple(ids10):
        raise ValueError("imported full model degree-10 coordinate prefix differs")

    p = next(iter(next(iter(fields96.values())))).p
    if p is None:
        raise ValueError("imported model is not over a finite field")
    lower_generators = set(fields10)
    if not lower_generators <= set(fields96):
        raise ValueError("imported model is missing a frozen degree-10 generator")

    for generator, expected_field in fields10.items():
        actual_field = fields96[generator]
        for i in range(24):
            terms = {}
            for powers, coefficient in actual_field[i].terms.items():
                if any(powers[24:]):
                    raise ValueError(
                        f"{generator}/{ids10[i]} depends on degree-12 couplings"
                    )
                short = tuple(powers[:24])
                terms[short] = (terms.get(short, 0) + int(coefficient)) % p
            actual = Poly(24, terms, p)
            expected = Poly(24, dict(expected_field[i].terms), p)
            if actual != expected:
                raise ValueError(
                    f"predecessor degree<=10 normalization mismatch at prime {p}: "
                    f"{generator}/{ids10[i]}"
                )

    degree12_only = set(fields96) - lower_generators
    for generator in degree12_only:
        for i in range(24):
            if fields96[generator][i]:
                raise ValueError(
                    f"degree-12-only generator {generator} changes lower "
                    f"coordinate {ids10[i]} at prime {p}"
                )

    return {
        "status": "lower_sector_matches_frozen_degree10_exactly_mod_p",
        "prime": int(p),
        "lower_generator_count": len(lower_generators),
        "degree12_only_generator_count": len(degree12_only),
    }


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

    rows = [
        [int(x) % p for x in row]
        for row in obstruction["obstructions_integer_coordinates"]
    ]
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
    out["source_prime"] = p
    out["source_kind"] = model.get("source_kind")
    out["source_certificate_engine_sha256"] = model.get(
        "source_certificate_engine_sha256"
    )
    out["source_degree12_sha256"] = model.get("source_degree12_sha256")
    out["source_all_holdouts_passed"] = bool(
        model.get("source_all_holdouts_passed")
    )
    return out


def import_predecessor_certificates(
    source_root: Path,
    existing_models: Path,
    output_dir: Path,
    obstruction_path: Path,
    degree10_model_path: Path,
) -> dict:
    source_root = Path(source_root)
    existing_models = Path(existing_models)
    output_dir = Path(output_dir)
    obstruction_path = Path(obstruction_path)
    degree10_model_path = Path(degree10_model_path)

    cert_dir = source_root / "results" / "stress_flow" / "certificates"
    if not cert_dir.is_dir():
        raise ValueError(f"certificate directory missing: {cert_dir}")
    obstruction = json.loads(obstruction_path.read_text())
    frozen_degree10 = json.loads(degree10_model_path.read_text())

    native_paths = sorted(existing_models.glob("fields_prime*.json"))
    if not native_paths:
        raise ValueError("at least one native reduced model is required")
    native_reference = json.loads(native_paths[0].read_text())
    if "source_prime" not in native_reference:
        raise ValueError("native reduced reference lacks source_prime")

    # Native inputs must themselves share theorem-neutral metadata and structure.
    native_support = _support_signature(native_reference)
    for path in native_paths[1:]:
        record = json.loads(path.read_text())
        if record["coordinate_ids"] != native_reference["coordinate_ids"]:
            raise ValueError(f"native coordinate order differs: {path}")
        if set(record["fields"]) != set(native_reference["fields"]):
            raise ValueError(f"native generator set differs: {path}")
        if _support_signature(record) != native_support:
            raise ValueError(f"native polynomial support differs: {path}")
        if record.get("generator_catalogue") != native_reference.get(
            "generator_catalogue"
        ):
            raise ValueError(f"native generator catalogue metadata differs: {path}")
        if record.get("restriction") != native_reference.get("restriction"):
            raise ValueError(f"native restriction metadata differs: {path}")

    output_dir.mkdir(parents=True, exist_ok=True)
    existing_primes = []
    for path in native_paths:
        target = output_dir / path.name
        shutil.copy2(path, target)
        existing_primes.append(int(path.stem.replace("fields_prime", "")))

    imported = []
    compatibility = []
    for p in EXPECTED_PREDECESSOR_PRIMES:
        cert_path = cert_dir / f"interacting_degree12_{p}.json"
        if not cert_path.exists():
            raise ValueError(f"missing predecessor certificate: {cert_path}")
        certificate = json.loads(cert_path.read_text())
        if int(certificate["prime"]) != p:
            raise ValueError(f"prime mismatch in {cert_path}")

        full = full_model_from_certificate(certificate, source_root)
        lower_check = verify_lower_sector_against_degree10(full, frozen_degree10)
        reduced = restrict_to_certified_degree10_orbit(full, obstruction)
        reduced = _align_reduced_metadata(reduced, native_reference)
        _assert_same_reduced_structure(reduced, native_reference)

        reduced["import_verification"] = {
            "predecessor_holdouts_passed": True,
            "lower_sector_check": lower_check,
            "native_structure_reference_prime": int(
                native_reference["source_prime"]
            ),
        }

        target = output_dir / f"fields_prime{p}.json"
        if target.exists():
            raise ValueError(f"prime collision at {p}: {target}")
        target.write_text(json.dumps(reduced, indent=2) + "\n")
        imported.append(p)
        compatibility.append(lower_check)

    return {
        "schema": 2,
        "status": "predecessor_degree12_certificates_imported_and_normalized",
        "native_reference_prime": int(native_reference["source_prime"]),
        "existing_primes": sorted(existing_primes),
        "imported_primes": imported,
        "total_primes": len(existing_primes) + len(imported),
        "lower_sector_compatibility": compatibility,
        "output_dir": str(output_dir),
    }
