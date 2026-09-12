"""Cross-prime structural/normalization audits for modular vector-field models."""
from __future__ import annotations

import hashlib
import json
from pathlib import Path
from typing import Mapping, Sequence


def _stable_json(value) -> str:
    return json.dumps(value, sort_keys=True, separators=(",", ":"))


def semantic_hash(value) -> str:
    return hashlib.sha256(_stable_json(value).encode()).hexdigest()


def _model_signature(record: Mapping) -> dict:
    ids = tuple(record.get("coordinate_ids", ()))
    fields = record.get("fields", {})
    if not ids or not isinstance(fields, dict) or not fields:
        raise ValueError("invalid vector-field model")
    generators = tuple(sorted(fields))
    shapes = {}
    prime = None
    for name in generators:
        comps = fields[name]
        if len(comps) != len(ids):
            raise ValueError(f"generator {name} has wrong component count")
        nvars = {int(x.get("nvars", -1)) for x in comps}
        primes = {x.get("prime") for x in comps}
        if nvars != {len(ids)} or len(primes) != 1:
            raise ValueError(f"generator {name} has inconsistent polynomial metadata")
        p = next(iter(primes))
        prime = p if prime is None else prime
        if p != prime:
            raise ValueError("mixed prime fields inside one model")
        # Term exponents, not coefficients, are structural and should align.
        shapes[name] = [tuple(tuple(row["powers"]) for row in comp.get("terms", ())) for comp in comps]
    catalog = record.get("generator_catalogue")
    if catalog is not None:
        catalog = sorted(
            (
                row.get("id"),
                tuple(row.get("factors", ())),
                row.get("leading_degree"),
            )
            for row in catalog
        )
    return {
        "coordinate_ids": ids,
        "generators": generators,
        "generator_catalogue": catalog,
        "restriction": record.get("restriction"),
        "field_prime": prime,
        "term_support_hash": semantic_hash(shapes),
    }


def audit_model_records(records_by_prime: Mapping[int, Mapping]) -> dict:
    """Require identical cross-prime structural normalization metadata."""
    if len(records_by_prime) < 2:
        raise ValueError("at least two prime models required")
    signatures = {int(p): _model_signature(record) for p, record in records_by_prime.items()}
    reference_prime = sorted(signatures)[0]
    reference = signatures[reference_prime]
    comparisons = {}
    for p, sig in sorted(signatures.items()):
        reference_restriction = reference["restriction"]
        current_restriction = sig["restriction"]
        restriction_consistent = (
            reference_restriction is None
            or current_restriction is None
            or current_restriction == reference_restriction
        )
        checks = {
            "coordinate_ids": sig["coordinate_ids"] == reference["coordinate_ids"],
            "generators": sig["generators"] == reference["generators"],
            "generator_catalogue": sig["generator_catalogue"] == reference["generator_catalogue"],
            "restriction_metadata_consistent": restriction_consistent,
            # Supports can legitimately differ if coefficients vanish at exceptional primes,
            # so record but do not gate theorem status on support hashes.
            "term_support_same": sig["term_support_hash"] == reference["term_support_hash"],
            "declared_prime_matches": sig["field_prime"] == int(p),
        }
        comparisons[p] = checks
        if not all(checks[k] for k in ("coordinate_ids", "generators", "generator_catalogue", "restriction_metadata_consistent", "declared_prime_matches")):
            raise ValueError(f"normalization mismatch at prime {p}: {checks}")
    return {
        "status": "normalization_audit_passed",
        "reference_prime": reference_prime,
        "primes": sorted(signatures),
        "coordinate_count": len(reference["coordinate_ids"]),
        "generator_count": len(reference["generators"]),
        "comparisons": comparisons,
    }


def audit_registry_records(registries: Sequence[Mapping]) -> dict:
    """Require registry definitions (including graphs/products/order) to match exactly."""
    if not registries:
        return {"status": "not_supplied"}
    fingerprints = [semantic_hash(r) for r in registries]
    if len(set(fingerprints)) != 1:
        raise ValueError("registry definitions differ across inputs")
    return {"status": "registry_definitions_identical", "fingerprint": fingerprints[0], "count": len(registries)}


def load_records(paths: Sequence[Path]) -> dict[int, dict]:
    out = {}
    for path in paths:
        record = json.loads(Path(path).read_text())
        # Prefer explicit source_prime, otherwise infer from polynomial metadata.
        p = record.get("source_prime")
        if p is None:
            fields = record.get("fields", {})
            if not fields:
                raise ValueError(f"cannot infer prime from {path}")
            first = next(iter(fields.values()))[0]
            p = first.get("prime")
        p = int(p)
        if p in out:
            raise ValueError(f"duplicate prime model {p}")
        out[p] = record
    return out
