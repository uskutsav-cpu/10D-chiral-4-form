"""Resolve a complete invariant registry for higher-cutoff theorem code.

The bundled data/fixtures/low_degree_registry.json is intentionally only a
curated degree<=8 excerpt. Higher-degree proofs must use a registry containing
the actual degree-10/12 bases, either from a completed local run artifact or
from the pinned predecessor source snapshot.
"""
from __future__ import annotations

from pathlib import Path

from .provenance import atomic_json
from .registry import Registry, import_predecessor_registry

EXPECTED_DIMS = {4: 1, 6: 2, 8: 7, 10: 14, 12: 72}

DEFAULT_SOURCE_ROOT = Path(".cache/upstream-degree12-fresh-20260911-223942")

CANDIDATE_REGISTRIES = (
    Path("runs/degree12-generalized-selected-verify/registry.json"),
    Path("runs/degree12-generalized-universal/registry.json"),
    Path("runs/degree12-fresh/registry.json"),
)


def registry_has_cutoff(reg: Registry, cutoff: int) -> bool:
    needed = {d: n for d, n in EXPECTED_DIMS.items() if d <= cutoff}
    return all(
        d in reg.degree_bases and len(reg.degree_bases[d]) == n
        for d, n in needed.items()
    )


def _try_load(path: Path, cutoff: int):
    if not path.exists():
        return None
    try:
        reg = Registry.load(path)
    except Exception:
        return None
    return reg if registry_has_cutoff(reg, cutoff) else None


def resolve_registry(
    cutoff: int,
    *,
    explicit_path: str | Path | None = None,
    source_root: str | Path = DEFAULT_SOURCE_ROOT,
) -> tuple[Registry, dict]:
    if cutoff not in EXPECTED_DIMS:
        raise ValueError("supported resolved cutoffs are 4,6,8,10,12")

    candidates = []
    if explicit_path is not None:
        candidates.append(Path(explicit_path))
    candidates.extend(CANDIDATE_REGISTRIES)

    for path in candidates:
        reg = _try_load(path, cutoff)
        if reg is not None:
            return reg, {
                "status": "resolved_from_existing_registry",
                "path": str(path),
                "cutoff": cutoff,
                "dimensions": {
                    str(d): len(reg.degree_bases[d])
                    for d in EXPECTED_DIMS
                    if d <= cutoff
                },
            }

    source_root = Path(source_root)
    if source_root.exists():
        reg = import_predecessor_registry(source_root, cutoff)
        if not registry_has_cutoff(reg, cutoff):
            raise ValueError(
                f"pinned predecessor import does not supply expected dimensions through {cutoff}"
            )
        return reg, {
            "status": "resolved_from_pinned_predecessor",
            "path": str(source_root),
            "cutoff": cutoff,
            "dimensions": {
                str(d): len(reg.degree_bases[d])
                for d in EXPECTED_DIMS
                if d <= cutoff
            },
        }

    attempted = [str(x) for x in candidates] + [str(source_root)]
    raise FileNotFoundError(
        "No complete registry found. Tried: " + ", ".join(attempted)
    )


def freeze_resolved_registry(
    cutoff: int,
    output: str | Path,
    *,
    explicit_path: str | Path | None = None,
    source_root: str | Path = DEFAULT_SOURCE_ROOT,
) -> Path:
    reg, provenance = resolve_registry(
        cutoff,
        explicit_path=explicit_path,
        source_root=source_root,
    )
    output = Path(output)
    output.parent.mkdir(parents=True, exist_ok=True)

    payload = reg.to_json()
    metadata = dict(payload.get("metadata", {}))
    metadata["registry_resolution"] = provenance
    payload["metadata"] = metadata

    atomic_json(output, payload)
    return output
