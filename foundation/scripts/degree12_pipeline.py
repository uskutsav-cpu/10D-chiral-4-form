#!/usr/bin/env python3
"""Exact, bounded, resumable degree-12 invariant discovery.

This pipeline deliberately makes no exhaustive order-12 graph-count claim.
It generates exact non-isomorphic nauty shards, verifies every stored
canonical ID, constructs a memory-gated deterministic schedule, and stops
only after the known degree-12 scalar-space upper bound is attained:

    10 lower products + 62 connected primitive polynomial directions = 72.

The same 62 primitives add 60 generic functional directions to the 21 lower
ones. The remaining two are retained as polynomially independent degree-12
invariants and explicitly certified as cumulative Jacobian dependencies.
"""

from __future__ import annotations

import argparse
from collections import Counter, defaultdict
import hashlib
import json
import os
from pathlib import Path
import platform
import resource
import shutil
import statistics
import sys
import time

import numpy as np

sys.path.insert(0, os.path.join(os.path.dirname(__file__), "..", "src"))

from sdinv.catalog import (
    atomic_write_json,
    canonical_graph_id,
    default_manifest_path,
    iter_graph_shard,
    load_manifest,
    write_graph_shard,
)
from sdinv.checkpoint import load_checkpoint, load_checkpoint_payload, write_checkpoint
from sdinv.contract import (
    build_compact_derivative_basis,
    contraction_plan_profile,
    greedy_contraction_plan_cost,
    planned_value,
    value_and_jacobian_row,
)
from sdinv.forms import (
    check_star_squared,
    metric_signs,
    random_form,
    selfdual_projector,
    to_dense,
)
from sdinv.graphs import (
    graph_from_label,
    graph_from_record,
    iter_graphs_nauty,
    validate_graph,
)
from sdinv.modp import ALT_P, P, RankSieve


D = 10
VALENCE = 5
ORDER = 12
MAX_MULTIPLICITY = 4
MIN_SUPPORT_DEGREE = 2
MAX_SUPPORT_DEGREE = 5
THIRD_P = 32693
DEGREE12_PRODUCTS = 10
DEGREE12_PRIMITIVES = 62
DEGREE12_DIMENSION = 72
LOWER_FUNCTIONAL_RANK = 21
CUMULATIVE_FUNCTIONAL_RANK = 81
DEFAULT_SHARDS = 4096
DEFAULT_RESIDUES = [63]
DEFAULT_CATALOG_DIR = Path("work/degree12-catalog")
DEFAULT_SCHEDULE = Path("work/degree12-schedule.json")
DEFAULT_CHECKPOINT = Path("work/degree12-discovery.checkpoint.json")
DEFAULT_LOWER_RESULT = Path("results/10d_order8.json")
DEFAULT_ORDER10_RESULT = Path("results/10d_order10.json")
DEFAULT_RESULT = Path("results/10d_order12.json")
DEFAULT_BENCHMARKS = Path("results/degree12_benchmarks.json")
PIPELINE_SCHEMA = 1


class RSSLimitError(MemoryError):
    """Raised when observed process RSS crosses the configured hard ceiling."""


def _file_sha256(path):
    digest = hashlib.sha256()
    with open(path, "rb") as stream:
        for block in iter(lambda: stream.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def _semantic_result_payload(value):
    """Remove run-time-only measurements from a result recursively."""
    volatile = {
        "seconds",
        "elapsed_seconds",
        "evaluation_seconds",
        "validation_seconds",
        "peak_rss_bytes",
        "validation_metrics",
    }
    if isinstance(value, dict):
        return {
            key: _semantic_result_payload(item)
            for key, item in value.items()
            if key not in volatile
        }
    if isinstance(value, list):
        return [_semantic_result_payload(item) for item in value]
    return value


def semantic_result_sha256(result):
    encoded = json.dumps(
        _semantic_result_payload(result),
        sort_keys=True,
        separators=(",", ":"),
        ensure_ascii=True,
    ).encode()
    return hashlib.sha256(encoded).hexdigest()


def _engine_sha256():
    root = Path(__file__).resolve().parents[1]
    digest = hashlib.sha256()
    for relative in (
        "scripts/degree12_pipeline.py",
        "src/sdinv/catalog.py",
        "src/sdinv/checkpoint.py",
        "src/sdinv/contract.py",
        "src/sdinv/forms.py",
        "src/sdinv/graphs.py",
        "src/sdinv/modp.py",
    ):
        path = root / relative
        digest.update(relative.encode())
        digest.update(path.read_bytes())
    return digest.hexdigest()


def _peak_rss_bytes():
    value = resource.getrusage(resource.RUSAGE_SELF).ru_maxrss
    return int(value if sys.platform == "darwin" else value * 1024)


def _guard_rss(limit, context):
    observed = _peak_rss_bytes()
    if observed > int(limit):
        raise RSSLimitError(
            f"RSS safety limit exceeded after {context}: observed "
            f"{observed} bytes > {int(limit)} bytes")
    return observed


def _distribution(values):
    ordered = sorted(float(value) for value in values)
    if not ordered:
        return {"count": 0}

    def percentile(fraction):
        return ordered[round(fraction * (len(ordered) - 1))]

    return {
        "count": len(ordered),
        "min": ordered[0],
        "p50": percentile(0.50),
        "p90": percentile(0.90),
        "p99": percentile(0.99),
        "max": ordered[-1],
        "mean": statistics.fmean(ordered),
    }


def _physical_memory_bytes():
    try:
        return int(os.sysconf("SC_PHYS_PAGES")) * int(
            os.sysconf("SC_PAGE_SIZE"))
    except (AttributeError, OSError, ValueError):
        return None


def _hardware_metadata():
    architecture = platform.machine() or "unknown"
    return {
        "architecture": architecture,
        "processor": platform.processor() or architecture,
        "logical_cpu_count": os.cpu_count(),
        "memory_bytes": _physical_memory_bytes(),
        "operating_system": platform.platform(),
        "python": sys.version.split()[0],
        "numpy": np.__version__,
    }


def _executable_identity(command):
    resolved = shutil.which(command)
    if resolved is None:
        raise FileNotFoundError(f"executable not found on PATH: {command}")
    path = Path(resolved).resolve()
    return {"path": str(path), "sha256": _file_sha256(path)}


def shard_name(residue, modulus):
    return f"order12-shard-{residue:04d}-of-{modulus:04d}.jsonl.gz"


def _feature_signature(M):
    upper = M[np.triu_indices(M.shape[0], 1)]
    multiplicities = Counter(int(x) for x in upper if x)
    support_degrees = Counter(
        int(np.count_nonzero(M[row])) for row in range(M.shape[0]))
    return json.dumps(
        {
            "edge_multiplicity_histogram": [
                multiplicities.get(k, 0) for k in range(1, 5)
            ],
            "support_degree_histogram": [
                support_degrees.get(k, 0) for k in range(1, 6)
            ],
        },
        sort_keys=True,
        separators=(",", ":"),
    )


def _generation_identity(args, residue):
    geng = _executable_identity(args.geng)
    multig = _executable_identity(args.multig)
    return {
        "schema": PIPELINE_SCHEMA,
        "software": "nauty gtools (content-addressed executables)",
        "pipeline": (
            "geng -cq -d2 -D5 12 RESIDUE/MODULUS | "
            "multig -q -T -m4 -r5"
        ),
        "geng": geng,
        "multig": multig,
        "order": ORDER,
        "valence": VALENCE,
        "max_multiplicity": MAX_MULTIPLICITY,
        "connected_only": True,
        "min_support_degree": MIN_SUPPORT_DEGREE,
        "max_support_degree": MAX_SUPPORT_DEGREE,
        "residue": int(residue),
        "modulus": int(args.shards),
        "claim": "exact non-isomorphic shard; no exhaustive global count claim",
    }


def _audit_shard(path):
    """Verify logical hash, every canonical ID, and within-shard uniqueness."""
    started = time.perf_counter()
    seen = set()
    for record, _ in iter_graph_shard(path, verify=True):
        if record["id"] in seen:
            raise ValueError(
                f"duplicate canonical graph ID in {path}: {record['id']}")
        seen.add(record["id"])
    manifest_path = default_manifest_path(path)
    manifest = load_manifest(manifest_path)
    if len(seen) != manifest["count"]:
        raise ValueError("shard unique-ID count does not match manifest")
    manifest["unique_canonical_ids"] = len(seen)
    manifest["verification_seconds"] = round(
        time.perf_counter() - started, 6)
    atomic_write_json(manifest_path, manifest)
    return manifest


def _enrich_candidate(M):
    peak, total = greedy_contraction_plan_cost(M, D, VALENCE)
    return {
        "greedy_peak_work": int(peak),
        "greedy_total_work": int(total),
        "feature_signature": _feature_signature(M),
    }


def command_generate(args):
    catalog_dir = Path(args.catalog_dir)
    catalog_dir.mkdir(parents=True, exist_ok=True)
    for residue in args.residues:
        if not 0 <= residue < args.shards:
            raise ValueError("residue must satisfy 0 <= residue < shards")
        shard = catalog_dir / shard_name(residue, args.shards)
        manifest_path = default_manifest_path(shard)
        identity = _generation_identity(args, residue)
        if args.resume and shard.exists() and manifest_path.exists():
            manifest = load_manifest(manifest_path)
            if manifest["generation"] != identity:
                raise ValueError(f"existing shard identity mismatch: {shard}")
            manifest = _audit_shard(shard)
            print(
                f"verified shard {residue}/{args.shards}: "
                f"{manifest['count']} candidates, all canonical IDs unique",
                flush=True,
            )
            continue
        graphs = iter_graphs_nauty(
            ORDER,
            VALENCE,
            MAX_MULTIPLICITY,
            geng=args.geng,
            multig=args.multig,
            connected_only=True,
            min_degree=MIN_SUPPORT_DEGREE,
            max_degree=MAX_SUPPORT_DEGREE,
            residue=residue,
            modulus=args.shards,
        )
        manifest = write_graph_shard(
            shard, graphs, identity, enrich=_enrich_candidate)
        manifest = _audit_shard(shard)
        print(
            f"wrote shard {residue}/{args.shards}: {manifest['count']} "
            f"unique candidates, sha256 {manifest['logical_sha256']}, "
            f"{manifest['seconds']:.3f}s",
            flush=True,
        )


def _load_source_shards(catalog_dir, shards, residues):
    sources = []
    for residue in residues:
        path = Path(catalog_dir) / shard_name(residue, shards)
        manifest_path = default_manifest_path(path)
        manifest = load_manifest(manifest_path)
        generation = manifest["generation"]
        if (
            generation["residue"] != residue
            or generation["modulus"] != shards
            or generation["order"] != ORDER
        ):
            raise ValueError(f"shard identity mismatch: {path}")
        if manifest.get("unique_canonical_ids") != manifest["count"]:
            raise ValueError(
                f"shard lacks a complete unique-ID audit: {path}; "
                "rerun generate --resume")
        sources.append({
            "residue": int(residue),
            "path": str(path),
            "count": int(manifest["count"]),
            "logical_sha256": manifest["logical_sha256"],
            "compressed_bytes": int(manifest["compressed_bytes"]),
            "generation_seconds": float(manifest["seconds"]),
        })
    return sources


def _diversity_round_robin(candidates):
    groups = defaultdict(list)
    for item in candidates:
        groups[item["feature_signature"]].append(item)
    for items in groups.values():
        items.sort(key=lambda item: (
            item["plan"]["largest_pair_work"],
            item["plan"]["total_forward_reverse_work"],
            item["graph_id"],
        ))
    signatures = sorted(
        groups,
        key=lambda signature: (
            groups[signature][0]["plan"]["largest_pair_work"],
            groups[signature][0]["plan"]["total_forward_reverse_work"],
            signature,
        ),
    )
    ordered = []
    while signatures:
        remaining = []
        for signature in signatures:
            items = groups[signature]
            ordered.append(items.pop(0))
            if items:
                remaining.append(signature)
        signatures = remaining
    return ordered


def command_schedule(args):
    sources = _load_source_shards(
        args.catalog_dir, args.shards, args.residues)
    existing = None
    prefix = []
    prior_examined = 0
    prior_rejected = Counter()
    schedule_path = Path(args.schedule)
    if args.extend_existing and schedule_path.exists():
        with schedule_path.open() as stream:
            existing = json.load(stream)
        old_sources = existing["sources"]
        if [
            source["logical_sha256"] for source in sources[:len(old_sources)]
        ] != [
            source["logical_sha256"] for source in old_sources
        ]:
            raise ValueError(
                "extended schedule sources do not preserve the old prefix")
        prefix = list(existing["candidates"])
        prior_examined = int(existing["examined"])
        prior_rejected.update(existing["rejected"])
        if args.target_safe <= len(prefix):
            raise ValueError(
                "--target-safe must exceed the existing schedule length")

    started_all = time.perf_counter()
    eligible_new = []
    seen = set()
    examined = 0
    sources_to_scan = sources
    if existing is not None and args.extension_new_sources_only:
        if len(sources) <= len(existing["sources"]):
            raise ValueError(
                "--extension-new-sources-only requires at least one added "
                "source shard")
        sources_to_scan = sources[len(existing["sources"]):]
        examined = prior_examined
        seen.update(item["graph_id"] for item in prefix)
    planner_seconds = []
    rejected = Counter(prior_rejected)
    stop = False
    for source in sources_to_scan:
        for record, M in iter_graph_shard(source["path"], verify=True):
            if record["id"] in seen:
                raise ValueError(
                    f"duplicate canonical ID across shards: {record['id']}")
            seen.add(record["id"])
            examined += 1
            if examined <= prior_examined:
                continue
            if examined >= args.candidate_limit:
                stop = True
                break
            started = time.perf_counter()
            profile = contraction_plan_profile(M, D, VALENCE)
            seconds = time.perf_counter() - started
            planner_seconds.append(seconds)
            if profile["max_output_rank"] > args.max_output_rank:
                rejected["max_output_rank"] += 1
                continue
            if profile["max_pair_union_rank"] > args.max_pair_union_rank:
                rejected["max_pair_union_rank"] += 1
                continue
            if profile["estimated_peak_bytes"] > args.max_memory_bytes:
                rejected["memory_estimate"] += 1
                continue
            validate_graph(M, VALENCE, MAX_MULTIPLICITY, True)
            eligible_new.append({
                "graph_id": record["id"],
                "graph": record["label"],
                "upper_triangle": record["upper_triangle"],
                "residue": int(source["residue"]),
                "shard_index": int(record["shard_index"]),
                "greedy_plan_cost": [
                    int(record["greedy_peak_work"]),
                    int(record["greedy_total_work"]),
                ],
                "feature_signature": record["feature_signature"],
                "plan": profile,
                "planner_seconds": round(seconds, 9),
            })
            if len(prefix) + len(eligible_new) >= args.target_safe:
                stop = True
                break
            if examined % args.progress_every == 0:
                print(
                    f"planned {examined}, eligible "
                    f"{len(prefix) + len(eligible_new)}, "
                    f"peak RSS {_peak_rss_bytes() / 2**20:.1f} MiB",
                    flush=True,
                )
        if stop:
            break
    eligible = prefix + eligible_new
    if len(eligible) < args.minimum_safe:
        raise ValueError(
            f"only {len(eligible)} safe candidates found; need at least "
            f"{args.minimum_safe}; add shards or raise --candidate-limit")
    ordered = prefix + _diversity_round_robin(eligible_new)
    identity = {
        "schema": PIPELINE_SCHEMA,
        "order": ORDER,
        "shards": int(args.shards),
        "residues": [int(x) for x in args.residues],
        "source_logical_sha256": [
            source["logical_sha256"] for source in sources
        ],
        "source_counts": [source["count"] for source in sources],
        "engine_sha256": _engine_sha256(),
        "planner": "globally optimal subset dynamic program",
        "ordering": (
            "exact-cost-gated diversity round-robin by multiplicity/support "
            "histogram; exact cost and canonical ID within each group"
        ),
        "max_output_rank": int(args.max_output_rank),
        "max_pair_union_rank": int(args.max_pair_union_rank),
        "max_memory_bytes": int(args.max_memory_bytes),
        "candidate_limit": int(args.candidate_limit),
    }
    schedule = {
        "schema": PIPELINE_SCHEMA,
        "claim": (
            "Bounded exact schedule from verified non-isomorphic shards; "
            "not an exhaustive order-12 catalog."
        ),
        "identity": identity,
        "sources": sources,
        "examined": examined,
        "eligible": len(eligible),
        "rejected": dict(sorted(rejected.items())),
        "planner_seconds": _distribution(planner_seconds),
        "plan_width_histogram": dict(sorted(Counter(
            str(item["plan"]["max_output_rank"]) for item in eligible
        ).items())),
        "pair_union_rank_histogram": dict(sorted(Counter(
            str(item["plan"]["max_pair_union_rank"]) for item in eligible
        ).items())),
        "estimated_peak_bytes": _distribution(
            item["plan"]["estimated_peak_bytes"] for item in eligible),
        "seconds": round(time.perf_counter() - started_all, 6),
        "candidates": ordered,
    }
    if existing is not None:
        schedule["extended_from"] = {
            "schedule_sha256": _file_sha256(args.schedule),
            "engine_sha256": existing["identity"]["engine_sha256"],
            "candidate_prefix_count": len(prefix),
            "examined_prefix_count": prior_examined,
        }
    atomic_write_json(args.schedule, schedule)
    print(
        f"schedule wrote {len(ordered)} safe candidates after examining "
        f"{examined}: {args.schedule}",
        flush=True,
    )


def _preflight(mod):
    if check_star_squared(D, VALENCE, True, mod) != 1:
        raise ValueError("self-duality is not well posed")
    projector = selfdual_projector(D, VALENCE, True, mod)
    if not np.array_equal((projector @ projector) % mod, projector % mod):
        raise ValueError("self-dual projector is not idempotent")
    return projector


def _sample_dense(projector, mod, seed):
    compact = (
        projector @ random_form(
            D, VALENCE, np.random.default_rng(seed), mod)
    ) % mod
    return to_dense(compact, D, VALENCE, mod)


def _contexts(mod, seeds):
    projector = _preflight(mod)
    basis = build_compact_derivative_basis(
        D, VALENCE, projector, mod, independent=True)
    if basis.ncols != 126:
        raise ValueError(
            f"self-dual tangent dimension is {basis.ncols}, expected 126")
    return projector, basis, [
        _sample_dense(projector, mod, seed) for seed in seeds
    ]


def _load_generators(path):
    with open(path) as stream:
        return json.load(stream)["generators"]


def _lower_inventory(lower_result, order10_result):
    lower = _load_generators(lower_result)
    degree10 = _load_generators(order10_result)
    if [item["order"] for item in lower] != [4, 6, 6] + [8] * 6:
        raise ValueError("unexpected order-8 generator inventory")
    if len(degree10) != 12 or any(item["order"] != 10 for item in degree10):
        raise ValueError("unexpected order-10 generator inventory")
    return lower, degree10


def _product_rows(values, rows, mod):
    """The ten degree-12 lower-product values and exact gradients."""
    I4 = values["I4_1"]
    dI4 = rows["I4_1"]
    I61, I62 = values["I6_1"], values["I6_2"]
    dI61, dI62 = rows["I6_1"], rows["I6_2"]
    products = [
        (
            "I4_1^3",
            I4 * I4 * I4 % mod,
            3 * I4 * I4 * dI4 % mod,
        ),
        (
            "I6_1^2",
            I61 * I61 % mod,
            2 * I61 * dI61 % mod,
        ),
        (
            "I6_1*I6_2",
            I61 * I62 % mod,
            (I61 * dI62 + I62 * dI61) % mod,
        ),
        (
            "I6_2^2",
            I62 * I62 % mod,
            2 * I62 * dI62 % mod,
        ),
    ]
    for k in range(1, 7):
        name = f"I8_{k}"
        value8, row8 = values[name], rows[name]
        products.append((
            f"I4_1*{name}",
            I4 * value8 % mod,
            (I4 * row8 + value8 * dI4) % mod,
        ))
    if len(products) != DEGREE12_PRODUCTS:
        raise AssertionError("wrong degree-12 product inventory")
    return products


def _evaluate_graphs(items, F_dense, basis, mod, max_memory_bytes):
    values = {}
    rows = {}
    seconds = []
    for item in items:
        M = graph_from_label(item["graph"])
        started = time.perf_counter()
        scalar, row = value_and_jacobian_row(
            M,
            F_dense,
            basis,
            D,
            VALENCE,
            True,
            mod,
            backend="optimized",
            max_memory_bytes=max_memory_bytes,
        )
        seconds.append(time.perf_counter() - started)
        values[item["id"]] = int(scalar)
        rows[item["id"]] = row
    return values, rows, seconds


def _seed_discovery(
    lower,
    degree10,
    samples,
    basis,
    mod,
    max_memory_bytes,
):
    sample_values = []
    sample_rows = []
    evaluation_seconds = []
    for F_dense in samples:
        values, rows, seconds = _evaluate_graphs(
            lower, F_dense, basis, mod, max_memory_bytes)
        sample_values.append(values)
        sample_rows.append(rows)
        evaluation_seconds.extend(seconds)

    cumulative = RankSieve(basis.ncols, mod)
    rank_seconds = []
    for item in lower:
        started = time.perf_counter()
        if not cumulative.add(sample_rows[0][item["id"]]):
            raise ValueError(f"dependent lower generator: {item['id']}")
        rank_seconds.append(time.perf_counter() - started)
    values10, rows10, seconds10 = _evaluate_graphs(
        degree10, samples[0], basis, mod, max_memory_bytes)
    evaluation_seconds.extend(seconds10)
    for item in degree10:
        started = time.perf_counter()
        if not cumulative.add(rows10[item["id"]]):
            raise ValueError(f"dependent degree-10 generator: {item['id']}")
        rank_seconds.append(time.perf_counter() - started)
    if cumulative.rank != LOWER_FUNCTIONAL_RANK:
        raise ValueError("lower cumulative rank is not 21")

    product_samples = [
        _product_rows(values, rows, mod)
        for values, rows in zip(sample_values, sample_rows)
    ]
    product_names = [item[0] for item in product_samples[0]]
    if any(
        [item[0] for item in products] != product_names
        for products in product_samples[1:]
    ):
        raise ValueError("product order differs between samples")
    polynomial = RankSieve(len(samples) * basis.ncols, mod)
    product_seconds = []
    for index, name in enumerate(product_names):
        started = time.perf_counter()
        stacked = np.concatenate([
            products[index][2] for products in product_samples
        ])
        if not polynomial.add(stacked):
            raise ValueError(
                f"degree-12 product is dependent across samples: {name}")
        product_seconds.append(time.perf_counter() - started)
    if polynomial.rank != DEGREE12_PRODUCTS:
        raise ValueError("degree-12 product rank is not 10")
    return {
        "polynomial": polynomial,
        "cumulative": cumulative,
        "product_names": product_names,
        "evaluation_seconds": evaluation_seconds,
        "rank_seconds": rank_seconds,
        "product_seconds": product_seconds,
    }


def _load_schedule(path, require_current_engine=True):
    with open(path) as stream:
        schedule = json.load(stream)
    if schedule.get("schema") != PIPELINE_SCHEMA:
        raise ValueError("unsupported degree-12 schedule schema")
    if (
        require_current_engine
        and schedule["identity"]["engine_sha256"] != _engine_sha256()
    ):
        raise ValueError(
            "schedule engine hash differs; rebuild schedule after code changes")
    return schedule


def _discovery_identity(args, schedule, basis):
    return {
        "schema": PIPELINE_SCHEMA,
        "order": ORDER,
        "prime": int(args.prime),
        "seeds": [int(seed) for seed in args.seeds],
        "schedule_sha256": _file_sha256(args.schedule),
        "schedule_identity": schedule["identity"],
        "lower_result_sha256": _file_sha256(args.lower_result),
        "order10_result_sha256": _file_sha256(args.order10_result),
        "engine_sha256": _engine_sha256(),
        "max_memory_bytes": int(args.max_memory_bytes),
        "rss_limit_bytes": int(args.rss_limit_bytes),
        "coordinate_indices": basis.coordinate_indices,
        "target_products": DEGREE12_PRODUCTS,
        "target_primitives": DEGREE12_PRIMITIVES,
        "target_degree12_dimension": DEGREE12_DIMENSION,
        "target_cumulative_rank": CUMULATIVE_FUNCTIONAL_RANK,
    }


def _load_or_migrate_discovery_checkpoint(path, identity, schedule):
    try:
        return load_checkpoint(path, identity)
    except ValueError as original_error:
        payload = load_checkpoint_payload(path)
        old_identity = payload["identity"]
        state = payload["state"]
        lineage = schedule.get("extended_from")
        if lineage is None:
            raise original_error
        if (
            lineage["schedule_sha256"] != old_identity["schedule_sha256"]
            or lineage["engine_sha256"] != old_identity["engine_sha256"]
            or state["cursor"] > lineage["candidate_prefix_count"]
        ):
            raise original_error
        ignored = {"schedule_sha256", "schedule_identity", "engine_sha256"}
        old_stable = {
            key: value for key, value in old_identity.items()
            if key not in ignored
        }
        new_stable = {
            key: value for key, value in identity.items()
            if key not in ignored
        }
        if old_stable != new_stable:
            raise original_error
        for item in state["generators"]:
            index = item["schedule_index"]
            if (
                index >= lineage["candidate_prefix_count"]
                or schedule["candidates"][index]["graph_id"]
                != item["graph_id"]
            ):
                raise original_error
        write_checkpoint(path, identity, state)
        print(
            "migrated checkpoint onto a verified schedule extension with "
            f"{lineage['candidate_prefix_count']} identical prefix entries",
            flush=True,
        )
        return state


def _checkpoint(path, identity, state, polynomial, cumulative, started):
    durable = dict(state)
    durable["polynomial_sieve"] = polynomial.to_state()
    durable["cumulative_sieve"] = cumulative.to_state()
    durable["elapsed_seconds"] = (
        float(state["elapsed_seconds"]) + time.perf_counter() - started
    )
    checkpoint_started = time.perf_counter()
    write_checkpoint(path, identity, durable)
    return time.perf_counter() - checkpoint_started


def command_discover(args):
    if len(args.seeds) < 2:
        raise ValueError("discovery needs at least two generic samples")
    schedule = _load_schedule(args.schedule)
    lower, degree10 = _lower_inventory(
        args.lower_result, args.order10_result)
    _, basis, samples = _contexts(args.prime, args.seeds)
    identity = _discovery_identity(args, schedule, basis)
    checkpoint_path = Path(args.checkpoint)
    if checkpoint_path.exists():
        state = _load_or_migrate_discovery_checkpoint(
            checkpoint_path, identity, schedule)
        polynomial = RankSieve.from_state(state["polynomial_sieve"])
        cumulative = RankSieve.from_state(state["cumulative_sieve"])
        print(
            f"resumed at schedule index {state['cursor']}, "
            f"{len(state['generators'])}/62 primitives, "
            f"ranks {polynomial.rank}/72 and {cumulative.rank}/81",
            flush=True,
        )
    else:
        seeded = _seed_discovery(
            lower,
            degree10,
            samples,
            basis,
            args.prime,
            args.max_memory_bytes,
        )
        polynomial = seeded["polynomial"]
        cumulative = seeded["cumulative"]
        state = {
            "cursor": 0,
            "evaluated": 0,
            "skipped_polynomial_dependent": 0,
            "skipped_early_functional_dependent": 0,
            "memory_rejections": 0,
            "elapsed_seconds": 0.0,
            "generators": [],
            "functional_dependencies": [],
            "product_names": seeded["product_names"],
            "lower_evaluation_seconds": seeded["evaluation_seconds"],
            "evaluation_seconds": [],
            "rank_seconds": seeded["rank_seconds"],
            "product_seconds": seeded["product_seconds"],
            "checkpoint_seconds": [],
            "rss_bytes": [_peak_rss_bytes()],
            "complete": False,
            "polynomial_sieve": polynomial.to_state(),
            "cumulative_sieve": cumulative.to_state(),
        }
        write_checkpoint(checkpoint_path, identity, state)
        print(
            "seeded exact ranks: 10 degree-12 products across samples; "
            "21 cumulative lower primitives",
            flush=True,
        )

    started = time.perf_counter()
    prior_elapsed = float(state["elapsed_seconds"])
    last_checkpoint = time.monotonic()
    candidates = schedule["candidates"]
    try:
        while len(state["generators"]) < DEGREE12_PRIMITIVES:
            if state["cursor"] >= len(candidates):
                break
            if (
                args.max_evaluated is not None
                and state["evaluated"] >= args.max_evaluated
            ):
                break
            candidate = candidates[state["cursor"]]
            M = graph_from_record({
                "order": ORDER,
                "upper_triangle": candidate["upper_triangle"],
            })
            if canonical_graph_id(M) != candidate["graph_id"]:
                raise ValueError("scheduled canonical graph ID mismatch")
            rows = []
            values = []
            per_sample_seconds = []
            try:
                for F_dense in samples:
                    eval_started = time.perf_counter()
                    scalar, row = value_and_jacobian_row(
                        M,
                        F_dense,
                        basis,
                        D,
                        VALENCE,
                        True,
                        args.prime,
                        backend="optimized",
                        max_memory_bytes=args.max_memory_bytes,
                    )
                    per_sample_seconds.append(
                        time.perf_counter() - eval_started)
                    values.append(int(scalar))
                    rows.append(row)
                    _guard_rss(
                        args.rss_limit_bytes,
                        f"discovery graph {candidate['graph_id']} sample",
                    )
            except RSSLimitError:
                raise
            except MemoryError:
                next_state = dict(state)
                next_state["memory_rejections"] = (
                    int(state["memory_rejections"]) + 1)
                next_state["cursor"] = int(state["cursor"]) + 1
                next_state["evaluated"] = int(state["evaluated"]) + 1
                state = next_state
                continue

            staged_polynomial = RankSieve.from_state(polynomial.to_state())
            rank_started = time.perf_counter()
            polynomial_increment = staged_polynomial.add(
                np.concatenate(rows))
            rank_seconds = time.perf_counter() - rank_started
            accepted = False
            functional_increment = False
            staged_cumulative = cumulative
            next_state = dict(state)
            for key in (
                "generators",
                "functional_dependencies",
                "evaluation_seconds",
                "rank_seconds",
                "checkpoint_seconds",
                "rss_bytes",
            ):
                next_state[key] = list(state[key])
            if polynomial_increment:
                staged_cumulative = RankSieve.from_state(
                    cumulative.to_state())
                rank_started = time.perf_counter()
                functional_increment = staged_cumulative.add(rows[0])
                rank_seconds += time.perf_counter() - rank_started
                # Build the 60-dimensional functional extension first.
                if (
                    cumulative.rank < CUMULATIVE_FUNCTIONAL_RANK
                    and not functional_increment
                ):
                    next_state["skipped_early_functional_dependent"] = (
                        int(state["skipped_early_functional_dependent"]) + 1)
                elif (
                    cumulative.rank == CUMULATIVE_FUNCTIONAL_RANK
                    and functional_increment
                ):
                    raise ValueError(
                        "cumulative rank exceeded known upper bound 81")
                else:
                    accepted = True
            else:
                next_state["skipped_polynomial_dependent"] = (
                    int(state["skipped_polynomial_dependent"]) + 1)

            next_state["cursor"] = int(state["cursor"]) + 1
            next_state["evaluated"] = int(state["evaluated"]) + 1
            next_state["evaluation_seconds"].extend(per_sample_seconds)
            next_state["rank_seconds"].append(rank_seconds)
            next_state["rss_bytes"].append(_peak_rss_bytes())
            if accepted:
                item = {
                    "id": f"I12_{len(state['generators']) + 1}",
                    "order": ORDER,
                    "graph_id": candidate["graph_id"],
                    "graph": candidate["graph"],
                    "residue": candidate["residue"],
                    "shard_index": candidate["shard_index"],
                    "schedule_index": next_state["cursor"] - 1,
                    "feature_signature": candidate["feature_signature"],
                    "greedy_plan_cost": candidate["greedy_plan_cost"],
                    "optimal_plan": candidate["plan"],
                    "discovery_values": values,
                    "polynomial_rank": staged_polynomial.rank,
                    "functional_increment": bool(functional_increment),
                    "cumulative_rank": staged_cumulative.rank,
                    "evaluation_seconds": round(
                        sum(per_sample_seconds), 6),
                }
                next_state["generators"].append(item)
                if not functional_increment:
                    next_state["functional_dependencies"].append(item["id"])
                # One assignment commits the complete candidate transaction.
                state, polynomial, cumulative = (
                    next_state,
                    staged_polynomial,
                    staged_cumulative,
                )
                checkpoint_seconds = _checkpoint(
                    checkpoint_path,
                    identity,
                    state,
                    polynomial,
                    cumulative,
                    started,
                )
                state["checkpoint_seconds"].append(checkpoint_seconds)
                last_checkpoint = time.monotonic()
                print(
                    f"NEW {item['id']}: polynomial rank {polynomial.rank}/72, "
                    f"cumulative rank {cumulative.rank}/81, "
                    f"functional {'yes' if functional_increment else 'NO'}, "
                    f"{sum(per_sample_seconds):.3f}s",
                    flush=True,
                )
            else:
                state = next_state
            if not accepted and (
                state["evaluated"] % args.checkpoint_every == 0
                or time.monotonic() - last_checkpoint >= 60
            ):
                checkpoint_seconds = _checkpoint(
                    checkpoint_path,
                    identity,
                    state,
                    polynomial,
                    cumulative,
                    started,
                )
                state["checkpoint_seconds"].append(checkpoint_seconds)
                last_checkpoint = time.monotonic()
            if state["evaluated"] % args.progress_every == 0:
                print(
                    f"evaluated {state['evaluated']}, selected "
                    f"{len(state['generators'])}/62, ranks "
                    f"{polynomial.rank}/72 and {cumulative.rank}/81, "
                    f"peak RSS {_peak_rss_bytes() / 2**20:.1f} MiB",
                    flush=True,
                )
    finally:
        state["complete"] = (
            len(state["generators"]) == DEGREE12_PRIMITIVES
            and polynomial.rank == DEGREE12_DIMENSION
            and cumulative.rank == CUMULATIVE_FUNCTIONAL_RANK
            and len(state["functional_dependencies"]) == 2
        )
        _checkpoint(
            checkpoint_path,
            identity,
            state,
            polynomial,
            cumulative,
            started,
        )

    elapsed = prior_elapsed + time.perf_counter() - started
    print(
        f"discovery {'complete' if state['complete'] else 'partial'}: "
        f"{len(state['generators'])}/62, ranks "
        f"{polynomial.rank}/72 and {cumulative.rank}/81, "
        f"{state['evaluated']} evaluated, {elapsed:.1f}s",
        flush=True,
    )
    if not state["complete"] and state["cursor"] >= len(candidates):
        raise RuntimeError(
            "safe schedule exhausted before completion; generate or schedule "
            "more exact candidates")


def _boost_matrix(mod, t=7):
    inverse_t = pow(t, mod - 2, mod)
    half = pow(2, mod - 2, mod)
    c = (t + inverse_t) * half % mod
    s = (inverse_t - t) * half % mod
    L = np.eye(D, dtype=np.int64)
    L[0, 0] = c
    L[0, 1] = s
    L[1, 0] = s
    L[1, 1] = c
    eta = np.diag(metric_signs(D, True)).astype(np.int64) % mod
    if not np.array_equal(L.T @ eta @ L % mod, eta):
        raise ValueError("finite-field boost does not preserve the metric")
    return L


def _transform_form(F_dense, L, mod):
    transformed = F_dense
    for axis in range(VALENCE):
        transformed = np.tensordot(
            L, transformed, axes=(1, axis)) % mod
        transformed = np.moveaxis(transformed, 0, axis)
    return transformed


def _validate_selected_ids(generators):
    seen = set()
    for item in generators:
        M = graph_from_label(item["graph"])
        validate_graph(M, VALENCE, MAX_MULTIPLICITY, True)
        actual = canonical_graph_id(M)
        if actual != item["graph_id"]:
            raise ValueError(f"canonical ID mismatch for {item['id']}")
        if actual in seen:
            raise ValueError(f"duplicate selected graph: {item['id']}")
        seen.add(actual)


def command_validate(args):
    prior_result = None
    if args.selection_result:
        with open(args.selection_result) as stream:
            prior_result = json.load(stream)
        generators = prior_result["generators"]
        prior_discovery = prior_result["discovery"]
        discovery_identity = prior_discovery["identity"]
        discovery = {
            "complete": True,
            "generators": generators,
            "functional_dependencies": prior_discovery[
                "functional_dependencies"],
            "product_names": [
                item["id"] for item in prior_result["products"]
            ],
            "evaluated": prior_discovery["evaluated"],
            "skipped_polynomial_dependent": prior_discovery[
                "skipped_polynomial_dependent"],
            "skipped_early_functional_dependent": prior_discovery[
                "skipped_early_functional_dependent"],
            "memory_rejections": prior_discovery["memory_rejections"],
            "elapsed_seconds": prior_discovery["seconds"],
        }
    elif args.checkpoint:
        payload = load_checkpoint_payload(args.checkpoint)
        discovery_identity = payload["identity"]
        discovery = payload["state"]
        generators = discovery["generators"]
    else:
        raise ValueError("provide --checkpoint or --selection-result")
    if not discovery["complete"] or len(generators) != DEGREE12_PRIMITIVES:
        raise ValueError("validation requires a complete discovery checkpoint")
    if len(discovery["functional_dependencies"]) != 2:
        raise ValueError("discovery did not identify exactly two dependencies")
    _validate_selected_ids(generators)
    lower, degree10 = _lower_inventory(
        args.lower_result, args.order10_result)
    selected_graphs = [
        graph_from_label(item["graph"]) for item in generators
    ]
    validation = {}
    all_started = time.perf_counter()
    validation_eval_seconds = []
    validation_rank_seconds = []
    validation_product_seconds = []
    lorentz_seconds = []

    for prime in args.primes:
        prime_started = time.perf_counter()
        projector, basis, samples = _contexts(prime, args.seeds)
        per_sample = []
        selected_rows_by_sample = []
        product_rows_by_sample = []
        selected_values_by_sample = []
        for seed, F_dense in zip(args.seeds, samples):
            sample_started = time.perf_counter()
            lower_values, lower_rows, seconds = _evaluate_graphs(
                lower, F_dense, basis, prime, args.max_memory_bytes)
            validation_eval_seconds.extend(seconds)
            values10, rows10, seconds = _evaluate_graphs(
                degree10, F_dense, basis, prime, args.max_memory_bytes)
            validation_eval_seconds.extend(seconds)
            cumulative = RankSieve(basis.ncols, prime)
            for item in lower:
                rank_started = time.perf_counter()
                if not cumulative.add(lower_rows[item["id"]]):
                    raise ValueError(
                        f"lower dependency at prime {prime}, seed {seed}: "
                        f"{item['id']}")
                validation_rank_seconds.append(
                    time.perf_counter() - rank_started)
            for item in degree10:
                rank_started = time.perf_counter()
                if not cumulative.add(rows10[item["id"]]):
                    raise ValueError(
                        f"degree-10 dependency at prime {prime}, seed {seed}: "
                        f"{item['id']}")
                validation_rank_seconds.append(
                    time.perf_counter() - rank_started)
            if cumulative.rank != LOWER_FUNCTIONAL_RANK:
                raise ValueError("lower cumulative rank is not 21")

            product_started = time.perf_counter()
            products = _product_rows(lower_values, lower_rows, prime)
            validation_product_seconds.append(
                time.perf_counter() - product_started)
            if any(value == 0 for _, value, _ in products):
                raise ValueError(
                    f"zero product value at prime {prime}, seed {seed}")
            product_rows_by_sample.append(
                [row for _, _, row in products])

            selected_rows = []
            selected_values = []
            dependencies = []
            first_sixty_all_increment = True
            for index, (item, M) in enumerate(zip(generators, selected_graphs)):
                eval_started = time.perf_counter()
                scalar, row = value_and_jacobian_row(
                    M,
                    F_dense,
                    basis,
                    D,
                    VALENCE,
                    True,
                    prime,
                    backend="optimized",
                    max_memory_bytes=args.max_memory_bytes,
                )
                validation_eval_seconds.append(
                    time.perf_counter() - eval_started)
                _guard_rss(
                    args.rss_limit_bytes,
                    f"validation {item['id']} prime {prime} seed {seed}",
                )
                if scalar == 0:
                    raise ValueError(
                        f"{item['id']} vanished at prime {prime}, seed {seed}")
                selected_values.append(int(scalar))
                selected_rows.append(row)
                rank_started = time.perf_counter()
                increment = cumulative.add(row)
                validation_rank_seconds.append(
                    time.perf_counter() - rank_started)
                if not increment:
                    dependencies.append(item["id"])
                    if index < 60:
                        first_sixty_all_increment = False
            if cumulative.rank != CUMULATIVE_FUNCTIONAL_RANK:
                raise ValueError(
                    f"cumulative rank {cumulative.rank}, expected 81 at "
                    f"prime {prime}, seed {seed}")
            if dependencies != discovery["functional_dependencies"]:
                raise ValueError(
                    f"dependency IDs changed at prime {prime}, seed {seed}: "
                    f"{dependencies}")
            if not first_sixty_all_increment:
                raise ValueError("one of the first 60 primitives is dependent")
            selected_rows_by_sample.append(selected_rows)
            selected_values_by_sample.append(selected_values)
            per_sample.append({
                "seed": int(seed),
                "lower_rank": LOWER_FUNCTIONAL_RANK,
                "degree12_functional_increments": 60,
                "dependency_ids": dependencies,
                "cumulative_rank": cumulative.rank,
                "all_product_values_nonzero": True,
                "all_primitive_values_nonzero": True,
                "seconds": round(time.perf_counter() - sample_started, 6),
            })
            print(
                f"prime {prime}, seed {seed}: cumulative rank 81 "
                f"(+60; dependencies {dependencies})",
                flush=True,
            )

        polynomial = RankSieve(len(samples) * basis.ncols, prime)
        product_pivots = []
        primitive_pivots = []
        product_names = [item[0] for item in products]
        for index, name in enumerate(product_names):
            rank_started = time.perf_counter()
            if not polynomial.add(np.concatenate([
                rows[index] for rows in product_rows_by_sample
            ])):
                raise ValueError(
                    f"stacked product dependency at prime {prime}: {name}")
            validation_rank_seconds.append(
                time.perf_counter() - rank_started)
            product_pivots.append(name)
        for index, item in enumerate(generators):
            rank_started = time.perf_counter()
            if not polynomial.add(np.concatenate([
                rows[index] for rows in selected_rows_by_sample
            ])):
                raise ValueError(
                    f"stacked primitive dependency at prime {prime}: "
                    f"{item['id']}")
            validation_rank_seconds.append(
                time.perf_counter() - rank_started)
            primitive_pivots.append(item["id"])
        if polynomial.rank != DEGREE12_DIMENSION:
            raise ValueError(
                f"degree-12 polynomial rank {polynomial.rank}, expected 72")
        print(
            f"prime {prime}: stacked four-sample degree-12 rank 72",
            flush=True,
        )

        L = _boost_matrix(prime)
        transformed = _transform_form(samples[0], L, prime)
        lorentz_matches = []
        for item, M, original in zip(
            generators, selected_graphs, selected_values_by_sample[0]
        ):
            lorentz_started = time.perf_counter()
            transformed_value = planned_value(
                M,
                transformed,
                D,
                VALENCE,
                True,
                prime,
                max_memory_bytes=args.max_memory_bytes,
            )
            lorentz_seconds.append(time.perf_counter() - lorentz_started)
            _guard_rss(
                args.rss_limit_bytes,
                f"Lorentz check {item['id']} prime {prime}",
            )
            if transformed_value != original:
                raise ValueError(
                    f"Lorentz check failed for {item['id']} at prime {prime}")
            lorentz_matches.append(item["id"])
        print(
            f"prime {prime}: all 62 exact finite-field boost checks pass",
            flush=True,
        )
        validation[str(prime)] = {
            "samples": per_sample,
            "degree12_polynomial_space": {
                "method": (
                    "exact rank of gradients stacked across four independent "
                    "generic samples"
                ),
                "columns": len(samples) * basis.ncols,
                "product_pivots": product_pivots,
                "primitive_pivots": primitive_pivots,
                "rank": polynomial.rank,
                "upper_bound": DEGREE12_DIMENSION,
            },
            "lorentz_boost": {
                "seed": int(args.seeds[0]),
                "matched": lorentz_matches,
            },
            "seconds": round(time.perf_counter() - prime_started, 6),
        }

    if args.calculation_note:
        pdf = Path(args.calculation_note)
        note = {
            "path_description": "independent supplied calculation note",
            "sha256": _file_sha256(pdf),
            "pages": 8,
            "use": (
                "target dimensions and independent cross-check only; no code "
                "or candidate list copied"
            ),
        }
    elif prior_result is not None:
        note = prior_result["literature"]["independent_calculation_note"]
    else:
        raise ValueError(
            "first-time validation requires --calculation-note; "
            "selection-result revalidation reuses its committed fingerprint")
    functional_dependencies = discovery["functional_dependencies"]
    degree12_basis = [
        {
            "id": name,
            "kind": "lower_product",
            "expression": name,
        }
        for name in discovery["product_names"]
    ] + [
        {
            "id": item["id"],
            "kind": "connected_primitive",
            "graph_id": item["graph_id"],
            "graph": item["graph"],
            "functional_increment": item["functional_increment"],
        }
        for item in generators
    ]
    primitive_candidates = [
        {
            "id": item["id"],
            "order": item["order"],
            "functionally_independent_in_selection": True,
        }
        for item in lower + degree10
    ] + [
        {
            "id": item["id"],
            "order": ORDER,
            "functionally_independent_in_selection": (
                item["id"] not in functional_dependencies
            ),
            "graph_id": item["graph_id"],
        }
        for item in generators
    ]
    result = {
        "schema": PIPELINE_SCHEMA,
        "claim": (
            "Complete degree-12 scalar basis relative to the known "
            "Hilbert-series upper bound: 10 lower products and 62 connected "
            "primitive polynomial directions give exact rank 72. Through "
            "degree 12, 83 primitive candidates have cumulative generic "
            "Jacobian rank 81; the two explicitly listed degree-12 "
            "candidates add no new functional direction."
        ),
        "proof_scope": {
            "arithmetic": (
                "exact finite-field arithmetic; no floating-point rank or "
                "tolerance"
            ),
            "degree12_completeness": (
                "72 independent degree-12 polynomials attain the supplied "
                "Hilbert-series upper bound 72"
            ),
            "functional_independence": (
                "rank 81 attained at four deterministic generic samples over "
                "each of three primes"
            ),
            "graph_generation": (
                "direct exact non-isomorphic nauty shards with canonical-ID "
                "verification; no exhaustive order-12 catalog claim"
            ),
            "probabilistic_clause": (
                "generic finite-field ranks are lower bounds in characteristic "
                "zero; agreement over three primes and four samples guards "
                "against exceptional reductions"
            ),
        },
        "literature": {
            "source": "https://arxiv.org/abs/2509.14350v2",
            "independent_calculation_note": note,
            "targets": {
                "degree12_scalar_dimension": DEGREE12_DIMENSION,
                "degree12_lower_products": DEGREE12_PRODUCTS,
                "degree12_connected_primitives": DEGREE12_PRIMITIVES,
                "primitive_candidates_through_degree12": 83,
                "cumulative_functional_rank": CUMULATIVE_FUNCTIONAL_RANK,
            },
        },
        "discovery": {
            "identity": discovery_identity,
            "evaluated": discovery["evaluated"],
            "selected": len(generators),
            "skipped_polynomial_dependent": discovery[
                "skipped_polynomial_dependent"],
            "skipped_early_functional_dependent": discovery[
                "skipped_early_functional_dependent"],
            "memory_rejections": discovery["memory_rejections"],
            "degree12_polynomial_rank": DEGREE12_DIMENSION,
            "cumulative_functional_rank": CUMULATIVE_FUNCTIONAL_RANK,
            "functional_dependencies": functional_dependencies,
            "seconds": discovery["elapsed_seconds"],
        },
        "catalog": {
            "claim": (
                "verified exact non-isomorphic search shards; no exhaustive "
                "order-12 catalog count"
            ),
            "shards": [
                {
                    "residue": residue,
                    "modulus": discovery_identity[
                        "schedule_identity"]["shards"],
                    "count": count,
                    "logical_sha256": digest,
                }
                for residue, count, digest in zip(
                    discovery_identity["schedule_identity"]["residues"],
                    discovery_identity["schedule_identity"]["source_counts"],
                    discovery_identity["schedule_identity"][
                        "source_logical_sha256"],
                )
            ],
            "schedule_sha256": discovery_identity["schedule_sha256"],
        },
        "products": [
            {"id": name, "expression": name}
            for name in discovery["product_names"]
        ],
        "generators": generators,
        "degree12_basis": degree12_basis,
        "primitive_candidates_through_degree12": primitive_candidates,
        "validation": validation,
        "validation_metrics": {
            "value_jacobian_seconds": _distribution(
                validation_eval_seconds),
            "rank_seconds": _distribution(validation_rank_seconds),
            "product_formation_seconds": _distribution(
                validation_product_seconds),
            "lorentz_scalar_seconds": _distribution(lorentz_seconds),
        },
        "validation_seconds": round(time.perf_counter() - all_started, 6),
        "validation_engine_sha256": _engine_sha256(),
        "peak_rss_bytes": _peak_rss_bytes(),
    }
    if len(primitive_candidates) != 83:
        raise AssertionError("wrong primitive candidate inventory")
    if sum(
        item["functionally_independent_in_selection"]
        for item in primitive_candidates
    ) != CUMULATIVE_FUNCTIONAL_RANK:
        raise AssertionError("wrong functional selection count")
    atomic_write_json(args.out, result)
    print(
        f"validated order-12 result written to {args.out}; peak RSS "
        f"{_peak_rss_bytes() / 2**20:.1f} MiB",
        flush=True,
    )


def command_benchmark(args):
    schedule = _load_schedule(args.schedule, require_current_engine=False)
    payload = load_checkpoint_payload(args.checkpoint)
    discovery = payload["state"]
    if (
        _file_sha256(args.schedule)
        != payload["identity"]["schedule_sha256"]
        or schedule["identity"]["engine_sha256"]
        != payload["identity"]["engine_sha256"]
    ):
        raise ValueError(
            "historical schedule does not match the discovery checkpoint")
    with open(args.result) as stream:
        result = json.load(stream)
    sample = schedule["candidates"][:args.metadata_sample]
    canonical_seconds = []
    serialization_seconds = []
    for item in sample:
        M = graph_from_record({
            "order": ORDER,
            "upper_triangle": item["upper_triangle"],
        })
        started = time.perf_counter()
        if canonical_graph_id(M) != item["graph_id"]:
            raise ValueError("benchmark canonical ID mismatch")
        canonical_seconds.append(time.perf_counter() - started)
        started = time.perf_counter()
        json.dumps(item, sort_keys=True, separators=(",", ":"))
        serialization_seconds.append(time.perf_counter() - started)

    width_histogram = Counter(
        str(item["plan"]["max_output_rank"])
        for item in schedule["candidates"]
    )
    pair_histogram = Counter(
        str(item["plan"]["max_pair_union_rank"])
        for item in schedule["candidates"]
    )
    benchmark = {
        "schema": PIPELINE_SCHEMA,
        "hardware": _hardware_metadata(),
        "schedule_sha256": _file_sha256(args.schedule),
        "result_sha256": _file_sha256(args.result),
        "result_semantic_sha256": semantic_result_sha256(result),
        "catalog_generation": {
            "sources": schedule["sources"],
            "seconds": _distribution(
                source["generation_seconds"]
                for source in schedule["sources"]
            ),
            "candidate_counts": [
                source["count"] for source in schedule["sources"]
            ],
        },
        "canonical_seconds": _distribution(canonical_seconds),
        "serialization_seconds": _distribution(serialization_seconds),
        "planner_seconds": schedule["planner_seconds"],
        "value_jacobian_seconds": _distribution(
            discovery["evaluation_seconds"]),
        "lower_value_jacobian_seconds": _distribution(
            discovery["lower_evaluation_seconds"]),
        "product_formation_and_seed_seconds": _distribution(
            discovery["product_seconds"]),
        "rank_update_seconds": _distribution(discovery["rank_seconds"]),
        "checkpoint_seconds": _distribution(
            discovery["checkpoint_seconds"]),
        "validation_metrics": result["validation_metrics"],
        "rss_bytes": {
            "samples": _distribution(discovery["rss_bytes"]),
            "discovery_peak": max(discovery["rss_bytes"]),
            "validation_peak": result["peak_rss_bytes"],
        },
        "plan_width_histogram": dict(sorted(width_histogram.items())),
        "pair_union_rank_histogram": dict(sorted(pair_histogram.items())),
        "estimated_peak_bytes": _distribution(
            item["plan"]["estimated_peak_bytes"]
            for item in schedule["candidates"]
        ),
        "selection": {
            "examined_by_scheduler": schedule["examined"],
            "safe_scheduled": schedule["eligible"],
            "evaluated": discovery["evaluated"],
            "selected": len(discovery["generators"]),
            "polynomial_dependencies_skipped": discovery[
                "skipped_polynomial_dependent"],
            "early_functional_dependencies_skipped": discovery[
                "skipped_early_functional_dependent"],
        },
    }
    atomic_write_json(args.out, benchmark)
    print(f"benchmark written to {args.out}", flush=True)


def build_parser():
    parser = argparse.ArgumentParser()
    subparsers = parser.add_subparsers(dest="command", required=True)

    generate = subparsers.add_parser(
        "generate", help="generate verified exact non-isomorphic shards")
    generate.add_argument("--catalog-dir", default=str(DEFAULT_CATALOG_DIR))
    generate.add_argument("--geng", default="geng")
    generate.add_argument("--multig", default="multig")
    generate.add_argument("--shards", type=int, default=DEFAULT_SHARDS)
    generate.add_argument(
        "--residues", type=int, nargs="+", default=DEFAULT_RESIDUES)
    generate.add_argument("--resume", action="store_true")
    generate.set_defaults(function=command_generate)

    schedule = subparsers.add_parser(
        "schedule", help="build exact cost/memory/diversity schedule")
    schedule.add_argument("--catalog-dir", default=str(DEFAULT_CATALOG_DIR))
    schedule.add_argument("--schedule", default=str(DEFAULT_SCHEDULE))
    schedule.add_argument("--shards", type=int, default=DEFAULT_SHARDS)
    schedule.add_argument(
        "--residues", type=int, nargs="+", default=DEFAULT_RESIDUES)
    schedule.add_argument("--candidate-limit", type=int, default=5000)
    schedule.add_argument("--target-safe", type=int, default=160)
    schedule.add_argument("--minimum-safe", type=int, default=80)
    schedule.add_argument("--max-output-rank", type=int, default=7)
    schedule.add_argument("--max-pair-union-rank", type=int, default=9)
    schedule.add_argument(
        "--max-memory-bytes", type=int, default=2 * 1024**3)
    schedule.add_argument("--progress-every", type=int, default=250)
    schedule.add_argument(
        "--extend-existing", action="store_true",
        help="preserve the existing candidate order as a verified prefix")
    schedule.add_argument(
        "--extension-new-sources-only", action="store_true",
        help="when extending, plan only newly appended source shards")
    schedule.set_defaults(function=command_schedule)

    discover = subparsers.add_parser(
        "discover", help="run or resume exact two-sample discovery")
    discover.add_argument("--schedule", default=str(DEFAULT_SCHEDULE))
    discover.add_argument("--checkpoint", default=str(DEFAULT_CHECKPOINT))
    discover.add_argument("--lower-result", default=str(DEFAULT_LOWER_RESULT))
    discover.add_argument(
        "--order10-result", default=str(DEFAULT_ORDER10_RESULT))
    discover.add_argument("--prime", type=int, default=P)
    discover.add_argument(
        "--seeds", type=int, nargs="+", default=[20260729, 20260730])
    discover.add_argument(
        "--max-memory-bytes", type=int, default=2 * 1024**3)
    discover.add_argument(
        "--rss-limit-bytes", type=int, default=3 * 1024**3)
    discover.add_argument("--max-evaluated", type=int)
    discover.add_argument("--checkpoint-every", type=int, default=5)
    discover.add_argument("--progress-every", type=int, default=10)
    discover.set_defaults(function=command_discover)

    validate = subparsers.add_parser(
        "validate", help="validate ranks over three primes and four samples")
    validate_source = validate.add_mutually_exclusive_group()
    validate_source.add_argument(
        "--checkpoint", default=str(DEFAULT_CHECKPOINT))
    validate_source.add_argument(
        "--selection-result",
        help="revalidate explicit formulas from a prior result without work/",
    )
    validate.add_argument("--lower-result", default=str(DEFAULT_LOWER_RESULT))
    validate.add_argument(
        "--order10-result", default=str(DEFAULT_ORDER10_RESULT))
    validate.add_argument(
        "--primes", type=int, nargs="+", default=[P, ALT_P, THIRD_P])
    validate.add_argument(
        "--seeds", type=int, nargs="+",
        default=[20260729, 20260730, 20260731, 20260732])
    validate.add_argument(
        "--max-memory-bytes", type=int, default=2 * 1024**3)
    validate.add_argument(
        "--rss-limit-bytes", type=int, default=3 * 1024**3)
    validate.add_argument(
        "--calculation-note",
        help="independent calculation-note PDF (required with --checkpoint)",
    )
    validate.add_argument("--out", default=str(DEFAULT_RESULT))
    validate.set_defaults(function=command_validate)

    benchmark = subparsers.add_parser(
        "benchmark", help="assemble recorded and measured stage benchmarks")
    benchmark.add_argument("--schedule", default=str(DEFAULT_SCHEDULE))
    benchmark.add_argument("--checkpoint", default=str(DEFAULT_CHECKPOINT))
    benchmark.add_argument("--result", default=str(DEFAULT_RESULT))
    benchmark.add_argument("--metadata-sample", type=int, default=100)
    benchmark.add_argument("--out", default=str(DEFAULT_BENCHMARKS))
    benchmark.set_defaults(function=command_benchmark)
    return parser


if __name__ == "__main__":
    arguments = build_parser().parse_args()
    arguments.function(arguments)
