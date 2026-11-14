#!/usr/bin/env python3
"""Checkpointed exact degree-10 trace-invariant pipeline.

The generation stage is exhaustive and sharded by nauty's canonical
support-graph partition. Discovery and validation are separate commands so a
long contraction run never endangers the catalog and can resume exactly.
"""

from __future__ import annotations

import argparse
from collections import Counter
import hashlib
import json
import os
from pathlib import Path
import platform
import resource
import shutil
import sqlite3
import statistics
import sys
import tempfile
import time
import numpy as np

sys.path.insert(0, os.path.join(os.path.dirname(__file__), "..", "src"))

from sdinv.catalog import (atomic_write_json, iter_graph_shard, load_manifest,
                           verify_graph_shard, write_graph_shard)
from sdinv.checkpoint import (load_checkpoint, load_checkpoint_payload,
                              write_checkpoint)
from sdinv.contract import (build_compact_derivative_basis,
                            contraction_plan_cost,
                            greedy_contraction_plan_cost, jacobian_row,
                            value, value_and_jacobian_row)
from sdinv.forms import (check_star_squared, random_form, selfdual_projector,
                         to_dense)
from sdinv.catalog import canonical_graph_id
from sdinv.graphs import (graph_from_label, graph_from_record,
                          iter_graphs_nauty)
from sdinv.modp import ALT_P, P, RankSieve


D = 10
VALENCE = 5
ORDER = 10
MAX_MULTIPLICITY = 4
MIN_SUPPORT_DEGREE = 2
MAX_SUPPORT_DEGREE = 5
CATALOG_INDEX_SCHEMA = 1
DEFAULT_SHARDS = 64
DEFAULT_CATALOG_DIR = Path("work/degree10-catalog")
DEFAULT_SCHEDULE = Path("work/degree10-schedule.sqlite3")
DEFAULT_CHECKPOINT = Path("work/degree10-discovery.checkpoint.json")
DEFAULT_LOWER_RESULT = Path("results/10d_order8.json")
DEFAULT_RESULT = Path("results/10d_order10.json")
DEFAULT_BENCHMARKS = Path("results/degree10_benchmarks.json")
SCHEDULE_SCHEMA = 1
DISCOVERY_SCHEMA = 1


def shard_name(residue, modulus):
    return f"order10-shard-{residue:04d}-of-{modulus:04d}.jsonl.gz"


def generation_identity(args, residue):
    geng = _executable_identity(args.geng)
    multig = _executable_identity(args.multig)
    return {
        "software": "nauty gtools (content-addressed executables)",
        "pipeline": (
            "geng -cq -d2 -D5 10 RESIDUE/MODULUS | "
            "multig -q -T -m4 -r5"
        ),
        "geng": geng["path"],
        "geng_sha256": geng["sha256"],
        "multig": multig["path"],
        "multig_sha256": multig["sha256"],
        "order": ORDER,
        "valence": VALENCE,
        "max_multiplicity": MAX_MULTIPLICITY,
        "connected_only": True,
        "min_support_degree": MIN_SUPPORT_DEGREE,
        "max_support_degree": MAX_SUPPORT_DEGREE,
        "residue": residue,
        "modulus": args.shards,
        "schedule_cost": "deterministic greedy pair contraction v1",
    }


def _executable_identity(command):
    """Resolve the executable exactly as subprocess will and content-address it."""
    resolved = shutil.which(command)
    if resolved is None:
        raise FileNotFoundError(f"executable not found on PATH: {command}")
    path = Path(resolved).resolve()
    return {"path": str(path), "sha256": _file_sha256(path)}


def enrich_candidate(M):
    peak, total = greedy_contraction_plan_cost(M, D, VALENCE)
    return {
        "greedy_peak_work": int(peak),
        "greedy_total_work": int(total),
    }


def _manifest_matches(manifest, generation):
    return manifest.get("generation") == generation


def build_catalog_index(catalog_dir, shards):
    """Verify every shard and exact ID, then write a global catalog index."""
    catalog_dir = Path(catalog_dir)
    seen = set()
    entries = []
    peak_histogram = Counter()
    total = 0
    common_generation = None
    started = time.perf_counter()
    logical_catalog_digest = hashlib.sha256()
    for residue in range(shards):
        shard = catalog_dir / shard_name(residue, shards)
        manifest_path = shard.with_suffix(shard.suffix + ".manifest.json")
        manifest = load_manifest(manifest_path)
        generation = manifest["generation"]
        if generation.get("residue") != residue:
            raise ValueError(f"wrong residue identity for {shard}")
        if generation.get("modulus") != shards:
            raise ValueError(f"wrong shard modulus identity for {shard}")
        generation_common = {
            key: value for key, value in generation.items()
            if key != "residue"
        }
        if common_generation is None:
            common_generation = generation_common
        elif generation_common != common_generation:
            raise ValueError(f"inconsistent generation identity for {shard}")
        count = 0
        for record, _ in iter_graph_shard(shard, verify=True):
            graph_id = record["id"]
            if graph_id in seen:
                raise ValueError(
                    f"duplicate graph ID across shards: {graph_id}")
            seen.add(graph_id)
            peak_histogram[str(record["greedy_peak_work"])] += 1
            count += 1
        if count != manifest["count"]:
            raise ValueError(f"manifest count mismatch for {shard}")
        entry = {
            "residue": residue,
            "path": shard.name,
            "count": count,
            "logical_sha256": manifest["logical_sha256"],
            "compressed_bytes": manifest["compressed_bytes"],
        }
        entries.append(entry)
        logical_entry = {
            "residue": residue,
            "count": count,
            "logical_sha256": manifest["logical_sha256"],
        }
        logical_catalog_digest.update(
            json.dumps(
                logical_entry, sort_keys=True, separators=(",", ":")).encode())
        total += count

    index = {
        "schema": CATALOG_INDEX_SCHEMA,
        "claim": (
            "Exhaustive exact catalog of connected, loop-free, valence-5 "
            "order-10 multigraphs with edge multiplicity at most 4."
        ),
        "order": ORDER,
        "valence": VALENCE,
        "max_multiplicity": MAX_MULTIPLICITY,
        "shards": shards,
        "candidate_count": total,
        "unique_canonical_ids": len(seen),
        "catalog_sha256": logical_catalog_digest.hexdigest(),
        "catalog_sha256_definition": (
            "SHA-256 over ordered canonical JSON objects containing only "
            "residue, logical record count, and logical shard SHA-256"
        ),
        "greedy_peak_work_histogram": dict(
            sorted(peak_histogram.items(), key=lambda item: int(item[0]))),
        "entries": entries,
        "generation": common_generation,
        "verification_seconds": round(time.perf_counter() - started, 6),
    }
    atomic_write_json(catalog_dir / "catalog-index.json", index)
    return index


def command_generate(args):
    catalog_dir = Path(args.catalog_dir)
    catalog_dir.mkdir(parents=True, exist_ok=True)
    residues = (
        list(range(args.shards))
        if args.residue is None
        else [args.residue]
    )
    if any(not 0 <= residue < args.shards for residue in residues):
        raise ValueError("residue must satisfy 0 <= residue < shards")

    for residue in residues:
        shard = catalog_dir / shard_name(residue, args.shards)
        manifest_path = shard.with_suffix(shard.suffix + ".manifest.json")
        generation = generation_identity(args, residue)
        if args.resume and shard.exists() and manifest_path.exists():
            manifest = load_manifest(manifest_path)
            if not _manifest_matches(manifest, generation):
                raise ValueError(
                    f"existing shard identity mismatch: {shard}")
            try:
                verify_graph_shard(shard)
            except (EOFError, OSError, ValueError) as error:
                print(
                    f"shard {residue}/{args.shards}: existing shard failed "
                    f"verification ({error}); rebuilding",
                    flush=True,
                )
            else:
                print(
                    f"shard {residue}/{args.shards}: verified existing "
                    f"{manifest['count']} candidates",
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
            shard,
            graphs,
            generation,
            enrich=enrich_candidate,
        )
        print(
            f"shard {residue}/{args.shards}: wrote {manifest['count']} "
            f"candidates, {manifest['compressed_bytes']} bytes, "
            f"{manifest['seconds']:.3f}s",
            flush=True,
        )

    if args.residue is None:
        index = build_catalog_index(catalog_dir, args.shards)
        print(
            f"catalog: {index['candidate_count']} candidates, "
            f"{index['unique_canonical_ids']} unique IDs, "
            f"sha256 {index['catalog_sha256']}",
            flush=True,
        )


def load_catalog_index(catalog_dir):
    path = Path(catalog_dir) / "catalog-index.json"
    with path.open() as stream:
        index = json.load(stream)
    if index.get("schema") != CATALOG_INDEX_SCHEMA:
        raise ValueError(f"unsupported catalog index schema in {path}")
    if index["candidate_count"] != index["unique_canonical_ids"]:
        raise ValueError("catalog index does not certify unique graph IDs")
    return index


def _durable_replace(temporary, destination):
    temporary = Path(temporary)
    destination = Path(destination)
    with temporary.open("rb") as stream:
        os.fsync(stream.fileno())
    os.replace(temporary, destination)
    directory_fd = os.open(destination.parent, os.O_RDONLY)
    try:
        os.fsync(directory_fd)
    finally:
        os.close(directory_fd)


def command_schedule(args):
    """Build a compact on-disk global cost order from verified graph shards."""
    catalog_dir = Path(args.catalog_dir)
    index = load_catalog_index(catalog_dir)
    schedule = Path(args.schedule)
    schedule.parent.mkdir(parents=True, exist_ok=True)
    temporary = schedule.with_name(f".{schedule.name}.{os.getpid()}.tmp")
    if temporary.exists():
        temporary.unlink()
    connection = sqlite3.connect(temporary)
    try:
        connection.execute("PRAGMA journal_mode=DELETE")
        connection.execute("PRAGMA synchronous=FULL")
        connection.execute(
            "CREATE TABLE metadata (key TEXT PRIMARY KEY, value TEXT NOT NULL)")
        connection.execute(
            """
            CREATE TABLE candidates (
                id TEXT PRIMARY KEY,
                greedy_peak_work INTEGER NOT NULL,
                greedy_total_work INTEGER NOT NULL,
                label TEXT NOT NULL,
                upper_triangle BLOB NOT NULL,
                residue INTEGER NOT NULL,
                shard_index INTEGER NOT NULL
            )
            """
        )
        identity = {
            "schema": SCHEDULE_SCHEMA,
            "catalog_sha256": index["catalog_sha256"],
            "candidate_count": index["candidate_count"],
            "ordering": (
                "greedy_peak_work,greedy_total_work,canonical_graph_id"
            ),
        }
        connection.execute(
            "INSERT INTO metadata VALUES (?, ?)",
            ("identity", json.dumps(identity, sort_keys=True)),
        )
        inserted = 0
        for entry in index["entries"]:
            shard = catalog_dir / entry["path"]
            manifest = load_manifest(
                shard.with_suffix(shard.suffix + ".manifest.json"))
            if (
                manifest["count"] != entry["count"]
                or manifest["logical_sha256"] != entry["logical_sha256"]
            ):
                raise ValueError(
                    f"shard no longer matches catalog index: {shard}; "
                    "rerun generate --resume before schedule")
            rows = []
            for record, _ in iter_graph_shard(shard, verify=True):
                rows.append((
                    record["id"],
                    int(record["greedy_peak_work"]),
                    int(record["greedy_total_work"]),
                    record["label"],
                    bytes(record["upper_triangle"]),
                    int(entry["residue"]),
                    int(record["shard_index"]),
                ))
            with connection:
                connection.executemany(
                    "INSERT INTO candidates VALUES (?, ?, ?, ?, ?, ?, ?)",
                    rows,
                )
            inserted += len(rows)
            print(
                f"schedule: imported shard {entry['residue']}, "
                f"{inserted}/{index['candidate_count']}",
                flush=True,
            )
        if inserted != index["candidate_count"]:
            raise ValueError("schedule candidate count mismatch")
        connection.execute(
            """
            CREATE UNIQUE INDEX schedule_order
            ON candidates(greedy_peak_work, greedy_total_work, id)
            """
        )
        connection.commit()
        integrity = connection.execute("PRAGMA integrity_check").fetchone()[0]
        if integrity != "ok":
            raise ValueError(f"SQLite integrity check failed: {integrity}")
    except BaseException:
        connection.close()
        try:
            temporary.unlink()
        except FileNotFoundError:
            pass
        raise
    connection.close()
    _durable_replace(temporary, schedule)
    print(
        f"schedule: wrote {inserted} candidates to {schedule} "
        f"({schedule.stat().st_size} bytes)",
        flush=True,
    )


def _file_sha256(path):
    digest = hashlib.sha256()
    with open(path, "rb") as stream:
        for block in iter(lambda: stream.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def _engine_sha256():
    root = Path(__file__).resolve().parents[1]
    digest = hashlib.sha256()
    for relative in (
        "scripts/degree10_pipeline.py",
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


def _preflight(mod):
    if check_star_squared(D, VALENCE, True, mod) != 1:
        raise ValueError("self-duality is not well posed")
    projector = selfdual_projector(D, VALENCE, True, mod)
    if not np.array_equal((projector @ projector) % mod, projector % mod):
        raise ValueError("self-dual projector is not idempotent")
    return projector


def _sample_and_basis(mod, seed):
    projector = _preflight(mod)
    F_dense = _sample_dense(projector, mod, seed)
    basis = build_compact_derivative_basis(
        D, VALENCE, projector, mod, independent=True)
    if basis.ncols != 126:
        raise ValueError(
            f"self-dual tangent basis has {basis.ncols} columns, expected 126")
    return F_dense, basis


def _sample_dense(projector, mod, seed):
    rng = np.random.default_rng(seed)
    compact = (
        projector @ random_form(D, VALENCE, rng, mod)
    ) % mod
    return to_dense(compact, D, VALENCE, mod)


def _seed_lower_degree_sieve(sieve, F_dense, basis, lower_result, mod):
    with open(lower_result) as stream:
        result = json.load(stream)
    generators = result["generators"]
    if [item["order"] for item in generators] != [4, 6, 6] + [8] * 6:
        raise ValueError("unexpected lower-degree generator inventory")
    for item in generators:
        M = graph_from_label(item["graph"])
        row = jacobian_row(
            M, F_dense, basis, D, VALENCE, True, mod, backend="optimized")
        if not sieve.add(row):
            raise ValueError(
                f"lower-degree generator is dependent: {item['id']}")
    if sieve.rank != 9:
        raise ValueError(f"lower-degree rank is {sieve.rank}, expected 9")
    return generators


def _schedule_identity(connection):
    row = connection.execute(
        "SELECT value FROM metadata WHERE key='identity'").fetchone()
    if row is None:
        raise ValueError("schedule has no identity")
    identity = json.loads(row[0])
    if identity.get("schema") != SCHEDULE_SCHEMA:
        raise ValueError("unsupported schedule schema")
    return identity


def _peak_rss_bytes():
    value = resource.getrusage(resource.RUSAGE_SELF).ru_maxrss
    # macOS reports bytes; Linux and most BSDs report KiB.
    if sys.platform == "darwin":
        return int(value)
    return int(value) * 1024


def _discovery_identity(args, schedule_identity, basis):
    return {
        "schema": DISCOVERY_SCHEMA,
        "dimension": D,
        "form_degree": VALENCE,
        "order": ORDER,
        "prime": int(args.prime),
        "seed": int(args.seed),
        "backend": "optimized",
        "catalog_sha256": schedule_identity["catalog_sha256"],
        "candidate_count": schedule_identity["candidate_count"],
        "schedule_ordering": schedule_identity["ordering"],
        "discovery_ordering": args.schedule_mode,
        "max_greedy_peak": (
            int(args.max_greedy_peak)
            if args.schedule_mode == "hash-bounded"
            else None
        ),
        "lower_result_sha256": _file_sha256(args.lower_result),
        "engine_sha256": _engine_sha256(),
        "coordinate_indices": basis.coordinate_indices,
    }


def _initial_discovery_state(sieve):
    return {
        "cursor": {
            "greedy_peak_work": -1,
            "greedy_total_work": -1,
            "id": "",
        },
        "evaluated": 0,
        "elapsed_seconds": 0.0,
        "sieve": sieve.to_state(),
        "generators": [],
        "complete": False,
    }


def _checkpoint_state(path, identity, state, sieve, started):
    durable_state = dict(state)
    durable_state["sieve"] = sieve.to_state()
    durable_state["elapsed_seconds"] = (
        float(state["elapsed_seconds"]) + time.perf_counter() - started
    )
    write_checkpoint(path, identity, durable_state)


def _stage_discovery_candidate(state, sieve, row, cursor, make_generator):
    """Prepare one complete discovery transaction without mutating live state.

    If planning or metadata construction is interrupted, callers retain the
    previous consistent ``(state, sieve)`` pair. The returned pair can then be
    committed with one reference assignment before it is checkpointed.
    """
    staged_sieve = RankSieve.from_state(sieve.to_state())
    independent = staged_sieve.add(row)
    generator = (
        make_generator(staged_sieve.rank) if independent else None
    )
    staged_state = dict(state)
    staged_state["generators"] = list(state["generators"])
    staged_state["evaluated"] = int(state["evaluated"]) + 1
    staged_state["cursor"] = dict(cursor)
    if generator is not None:
        staged_state["generators"].append(generator)
    return staged_state, staged_sieve, generator


def _scheduled_rows(connection, cursor, batch_size, schedule_mode,
                    max_greedy_peak):
    if schedule_mode == "hash-bounded":
        return connection.execute(
            """
            SELECT greedy_peak_work, greedy_total_work, id, label,
                   upper_triangle, residue, shard_index
            FROM candidates
            WHERE id > ? AND greedy_peak_work <= ?
            ORDER BY id
            LIMIT ?
            """,
            (cursor["id"], max_greedy_peak, batch_size),
        ).fetchall()
    if schedule_mode != "cost":
        raise ValueError(f"unknown discovery schedule mode: {schedule_mode}")
    return connection.execute(
        """
        SELECT greedy_peak_work, greedy_total_work, id, label,
               upper_triangle, residue, shard_index
        FROM candidates
        WHERE greedy_peak_work > ?
           OR (greedy_peak_work = ? AND greedy_total_work > ?)
           OR (greedy_peak_work = ? AND greedy_total_work = ? AND id > ?)
        ORDER BY greedy_peak_work, greedy_total_work, id
        LIMIT ?
        """,
        (
            cursor["greedy_peak_work"],
            cursor["greedy_peak_work"],
            cursor["greedy_total_work"],
            cursor["greedy_peak_work"],
            cursor["greedy_total_work"],
            cursor["id"],
            batch_size,
        ),
    ).fetchall()


def command_discover(args):
    schedule = sqlite3.connect(f"file:{Path(args.schedule).resolve()}?mode=ro",
                               uri=True)
    schedule_identity = _schedule_identity(schedule)
    catalog_index = load_catalog_index(args.catalog_dir)
    if schedule_identity["catalog_sha256"] != catalog_index["catalog_sha256"]:
        raise ValueError("schedule/catalog identity mismatch")

    F_dense, basis = _sample_and_basis(args.prime, args.seed)
    identity = _discovery_identity(args, schedule_identity, basis)
    checkpoint_path = Path(args.checkpoint)
    if checkpoint_path.exists():
        state = load_checkpoint(checkpoint_path, identity)
        sieve = RankSieve.from_state(state["sieve"])
        if sieve.rank != 9 + len(state["generators"]):
            raise ValueError("checkpoint generator/rank mismatch")
        print(
            f"resumed after {state['evaluated']} candidates at rank "
            f"{sieve.rank}, with {len(state['generators'])} degree-10 "
            "generators",
            flush=True,
        )
    else:
        sieve = RankSieve(basis.ncols, args.prime)
        _seed_lower_degree_sieve(
            sieve, F_dense, basis, args.lower_result, args.prime)
        state = _initial_discovery_state(sieve)
        write_checkpoint(checkpoint_path, identity, state)
        print("seeded exact lower-degree rank 9", flush=True)

    started = time.perf_counter()
    prior_elapsed = float(state["elapsed_seconds"])
    since_checkpoint = 0
    committed = (state, sieve, None)
    try:
        while len(state["generators"]) < args.target_new:
            if (
                args.max_candidates is not None
                and state["evaluated"] >= args.max_candidates
            ):
                break
            rows = _scheduled_rows(
                schedule,
                state["cursor"],
                min(args.batch_size, (
                    args.max_candidates - state["evaluated"]
                    if args.max_candidates is not None
                    else args.batch_size
                )),
                args.schedule_mode,
                args.max_greedy_peak,
            )
            if not rows:
                break
            for peak, total, graph_id, label, upper, residue, shard_index in rows:
                M = graph_from_record({
                    "order": ORDER,
                    "upper_triangle": list(upper),
                })
                candidate_started = time.perf_counter()
                scalar, row = value_and_jacobian_row(
                    M,
                    F_dense,
                    basis,
                    D,
                    VALENCE,
                    True,
                    args.prime,
                    backend="optimized",
                )
                seconds = time.perf_counter() - candidate_started
                cursor = {
                    "greedy_peak_work": int(peak),
                    "greedy_total_work": int(total),
                    "id": graph_id,
                }

                def make_generator(rank):
                    exact_peak, exact_total = contraction_plan_cost(
                        M, D, VALENCE)
                    return {
                        "id": f"I10_{len(state['generators']) + 1}",
                        "order": ORDER,
                        "graph_id": graph_id,
                        "graph": label,
                        "residue": int(residue),
                        "shard_index": int(shard_index),
                        "greedy_plan_cost": [int(peak), int(total)],
                        "optimal_plan_cost": [
                            int(exact_peak), int(exact_total)],
                        "discovery_value": int(scalar),
                        "discovery_rank": rank,
                        "evaluation_seconds": round(seconds, 6),
                    }

                committed = _stage_discovery_candidate(
                    state, sieve, row, cursor, make_generator)
                state, sieve, generator = committed
                since_checkpoint += 1
                if generator is not None:
                    print(
                        f"NEW {generator['id']}: rank {sieve.rank}, "
                        f"candidate {state['evaluated']}, {seconds:.3f}s, "
                        f"{label}",
                        flush=True,
                    )
                    _checkpoint_state(
                        checkpoint_path, identity,
                        committed[0], committed[1], started)
                    state["elapsed_seconds"] = prior_elapsed
                    since_checkpoint = 0
                elif state["evaluated"] % args.progress_every == 0:
                    elapsed = prior_elapsed + time.perf_counter() - started
                    print(
                        f"evaluated {state['evaluated']}, rank {sieve.rank}, "
                        f"last {seconds:.3f}s, elapsed {elapsed:.1f}s, "
                        f"peak RSS {_peak_rss_bytes() / 2**20:.1f} MiB",
                        flush=True,
                    )
                if since_checkpoint >= args.checkpoint_every:
                    _checkpoint_state(
                        checkpoint_path, identity,
                        committed[0], committed[1], started)
                    state["elapsed_seconds"] = prior_elapsed
                    since_checkpoint = 0
                if len(state["generators"]) >= args.target_new:
                    break
    finally:
        state, sieve, _ = committed
        final_state = dict(state)
        final_state["complete"] = (
            len(final_state["generators"]) == args.target_new)
        committed = (final_state, sieve, None)
        state = final_state
        _checkpoint_state(
            checkpoint_path, identity, committed[0], committed[1], started)
        schedule.close()

    elapsed = prior_elapsed + time.perf_counter() - started
    print(
        f"discovery {'complete' if state['complete'] else 'partial'}: "
        f"{len(state['generators'])}/{args.target_new} new, "
        f"rank {sieve.rank}, {state['evaluated']} evaluated, "
        f"{elapsed:.1f}s, peak RSS {_peak_rss_bytes() / 2**20:.1f} MiB",
        flush=True,
    )


def _load_generators(path):
    with open(path) as stream:
        return json.load(stream)["generators"]


def _degree10_value_row(F_dense, lower_graphs, degree10_graphs, mod):
    lower_values = [
        value(M, F_dense, D, VALENCE, True, mod)
        for M in lower_graphs[:3]
    ]
    degree10_values = [
        value(M, F_dense, D, VALENCE, True, mod)
        for M in degree10_graphs
    ]
    products = [
        (lower_values[0] * lower_values[1]) % mod,
        (lower_values[0] * lower_values[2]) % mod,
    ]
    return degree10_values + products


def _validate_catalog_selection(generators, catalog_dir, schedule_path):
    index = load_catalog_index(catalog_dir)
    connection = sqlite3.connect(
        f"file:{Path(schedule_path).resolve()}?mode=ro", uri=True)
    try:
        schedule_identity = _schedule_identity(connection)
        if schedule_identity["catalog_sha256"] != index["catalog_sha256"]:
            raise ValueError("validation schedule/catalog identity mismatch")
        for item in generators:
            row = connection.execute(
                "SELECT label, residue, shard_index FROM candidates WHERE id=?",
                (item["graph_id"],),
            ).fetchone()
            if row is None or row != (
                item["graph"], item["residue"], item["shard_index"]):
                raise ValueError(
                    f"saved catalog location mismatch for {item['id']}")
    finally:
        connection.close()
    return index


def _validate_saved_graph_ids(generators):
    """Check explicit graph/ID pairs without requiring the full catalog."""
    for item in generators:
        M = graph_from_label(item["graph"])
        if canonical_graph_id(M) != item["graph_id"]:
            raise ValueError(f"saved canonical ID mismatch for {item['id']}")


def command_validate(args):
    prior_result = None
    if args.selection_result:
        with open(args.selection_result) as stream:
            prior_result = json.load(stream)
        generators = prior_result["generators"]
        discovery_identity = {
            "prime": prior_result["discovery"]["prime"],
            "seed": prior_result["discovery"]["seed"],
            "discovery_ordering": prior_result["discovery"]["ordering"],
            "max_greedy_peak": prior_result["discovery"]["max_greedy_peak"],
        }
        discovery_state = {
            "generators": generators,
            "evaluated": prior_result["discovery"]["evaluated"],
            "elapsed_seconds": prior_result["discovery"]["seconds"],
            "sieve": {"rows": [None] * 21},
        }
    elif args.discovery_checkpoint:
        payload = load_checkpoint_payload(args.discovery_checkpoint)
        discovery_identity = payload["identity"]
        discovery_state = payload["state"]
        generators = discovery_state["generators"]
    else:
        raise ValueError(
            "provide --discovery-checkpoint or --selection-result")
    if len(generators) != 12:
        raise ValueError(
            f"discovery has {len(generators)} generators, expected 12")
    if len(discovery_state["sieve"]["rows"]) != 21:
        raise ValueError("discovery checkpoint does not have rank 21")

    _validate_saved_graph_ids(generators)
    if args.skip_catalog_check:
        if prior_result is None:
            raise ValueError(
                "--skip-catalog-check requires --selection-result")
        index = prior_result["catalog"]
    else:
        index = _validate_catalog_selection(
            generators, args.catalog_dir, args.schedule)
    lower_generators = _load_generators(args.lower_result)
    lower_graphs = [
        graph_from_label(item["graph"]) for item in lower_generators]
    degree10_graphs = [
        graph_from_label(item["graph"]) for item in generators]
    jacobian_seeds = [int(seed) for seed in args.jacobian_seeds]
    value_seeds = [
        int(args.value_seed_start + offset)
        for offset in range(args.value_samples)
    ]
    validation = {}
    started_all = time.perf_counter()

    for prime in args.primes:
        prime_started = time.perf_counter()
        projector = _preflight(prime)
        basis = build_compact_derivative_basis(
            D, VALENCE, projector, prime, independent=True)
        if basis.ncols != 126:
            raise ValueError("self-dual tangent dimension is not 126")

        jacobian_runs = []
        cached_value_rows = {}
        for seed in jacobian_seeds:
            sample_started = time.perf_counter()
            F_dense = _sample_dense(projector, prime, seed)
            sieve = RankSieve(basis.ncols, prime)
            order_ranks = {}
            scalar_values = {}
            for item, M in zip(lower_generators, lower_graphs):
                scalar, row = value_and_jacobian_row(
                    M, F_dense, basis, D, VALENCE, True, prime,
                    backend="optimized")
                if not sieve.add(row):
                    raise ValueError(
                        f"{item['id']} dependent at prime {prime}, seed {seed}")
                scalar_values[item["id"]] = int(scalar)
                order_ranks[str(item["order"])] = sieve.rank
            if order_ranks != {"4": 1, "6": 3, "8": 9}:
                raise ValueError(
                    f"wrong lower rank pattern at prime {prime}, seed {seed}")

            selected_values = []
            for item, M in zip(generators, degree10_graphs):
                scalar, row = value_and_jacobian_row(
                    M, F_dense, basis, D, VALENCE, True, prime,
                    backend="optimized")
                if scalar == 0:
                    raise ValueError(
                        f"{item['id']} vanished at prime {prime}, seed {seed}")
                if not sieve.add(row):
                    raise ValueError(
                        f"{item['id']} dependent at prime {prime}, seed {seed}")
                scalar_values[item["id"]] = int(scalar)
                selected_values.append(int(scalar))
            if sieve.rank != 21:
                raise ValueError(
                    f"cumulative rank is {sieve.rank}, expected 21")

            products = [
                scalar_values["I4_1"] * scalar_values["I6_1"] % prime,
                scalar_values["I4_1"] * scalar_values["I6_2"] % prime,
            ]
            cached_value_rows[seed] = selected_values + products
            jacobian_runs.append({
                "seed": seed,
                "lower_rank_pattern": {"4": 1, "6": 3, "8": 9},
                "new_degree10": 12,
                "cumulative_rank": 21,
                "all_degree10_values_nonzero": True,
                "seconds": round(time.perf_counter() - sample_started, 6),
            })
            print(
                f"prime {prime}, seed {seed}: Jacobian rank 21 "
                f"(+12 at degree 10)",
                flush=True,
            )

        reference_checks = []
        reference_seed = jacobian_seeds[0]
        F_reference = _sample_dense(projector, prime, reference_seed)
        if prime == args.primes[0]:
            for item, M in zip(generators, degree10_graphs):
                optimized = value_and_jacobian_row(
                    M, F_reference, basis, D, VALENCE, True, prime,
                    backend="optimized")
                reference = value_and_jacobian_row(
                    M, F_reference, basis, D, VALENCE, True, prime,
                    backend="reference")
                matched = (
                    optimized[0] == reference[0]
                    and np.array_equal(optimized[1], reference[1])
                )
                if not matched:
                    raise ValueError(
                        f"backend mismatch for {item['id']}")
                reference_checks.append(item["id"])
            print(
                f"prime {prime}: all 12 optimized/reference backend checks pass",
                flush=True,
            )

        value_sieve = RankSieve(14, prime)
        pivot_seeds = []
        for seed in value_seeds:
            row = cached_value_rows.get(seed)
            if row is None:
                F_dense = _sample_dense(projector, prime, seed)
                row = _degree10_value_row(
                    F_dense, lower_graphs, degree10_graphs, prime)
            if value_sieve.add(row):
                pivot_seeds.append(seed)
        if value_sieve.rank != 14:
            raise ValueError(
                f"degree-10 value rank is {value_sieve.rank}, expected 14")
        print(
            f"prime {prime}: degree-10 value rank 14 across "
            f"{len(value_seeds)} samples",
            flush=True,
        )
        validation[str(prime)] = {
            "jacobian_samples": jacobian_runs,
            "optimized_reference_matches": reference_checks,
            "value_space": {
                "basis": [
                    item["id"] for item in generators
                ] + ["I4_1*I6_1", "I4_1*I6_2"],
                "sample_seeds": value_seeds,
                "pivot_seeds": pivot_seeds,
                "rank": value_sieve.rank,
                "dimension": 14,
            },
            "seconds": round(time.perf_counter() - prime_started, 6),
        }

    result = {
        "schema": 1,
        "claim": (
            "Complete degree-10 trace basis: 12 connected primitive "
            "directions plus the two lower products I4_1*I6_1 and "
            "I4_1*I6_2. Exact value-space rank 14 and cumulative generic "
            "Jacobian rank 21 through degree 10."
        ),
        "proof_scope": {
            "arithmetic": "exact finite fields; no floating-point rank test",
            "primitive_completeness": (
                "12 independent connected directions attain the published "
                "Hilbert-series upper bound of 12 new degree-10 primitives"
            ),
            "catalog_enumeration": (
                "exhaustive connected graph catalog generated and verified; "
                "Hilbert-guided discovery stopped after rank target"
            ),
            "not_claimed": (
                "This is not a claim that all 81 invariants form a complete "
                "polynomial generating set."
            ),
        },
        "literature_targets": {
            "degree10_scalar_space_dimension": 14,
            "lower_product_directions": 2,
            "new_primitive_directions": 12,
            "cumulative_jacobian_rank": 21,
        },
        "catalog": {
            "directory": str(args.catalog_dir),
            "candidate_count": index["candidate_count"],
            "unique_canonical_ids": index["unique_canonical_ids"],
            "shards": index["shards"],
            "catalog_sha256": index["catalog_sha256"],
        },
        "discovery": {
            "prime": discovery_identity["prime"],
            "seed": discovery_identity["seed"],
            "ordering": discovery_identity["discovery_ordering"],
            "max_greedy_peak": discovery_identity["max_greedy_peak"],
            "evaluated": discovery_state["evaluated"],
            "seconds": discovery_state["elapsed_seconds"],
            "lower_rank": 9,
            "new_degree10": 12,
            "cumulative_rank": 21,
        },
        "generators": generators,
        "degree10_basis": [
            {
                "id": item["id"],
                "kind": "connected_primitive",
                "graph_id": item["graph_id"],
                "graph": item["graph"],
            }
            for item in generators
        ] + [
            {
                "id": "I4_1*I6_1",
                "kind": "lower_product",
                "expression": "I4_1*I6_1",
            },
            {
                "id": "I4_1*I6_2",
                "kind": "lower_product",
                "expression": "I4_1*I6_2",
            },
        ],
        "validation": validation,
        "validation_seconds": round(time.perf_counter() - started_all, 6),
        "peak_rss_bytes": _peak_rss_bytes(),
    }
    atomic_write_json(args.out, result)
    print(
        f"validated result written to {args.out}; "
        f"peak RSS {_peak_rss_bytes() / 2**20:.1f} MiB",
        flush=True,
    )


def _distribution(values):
    values = sorted(float(value) for value in values)
    if not values:
        return {"count": 0}

    def percentile(fraction):
        return values[round(fraction * (len(values) - 1))]

    return {
        "count": len(values),
        "min": values[0],
        "p50": percentile(0.50),
        "p90": percentile(0.90),
        "p99": percentile(0.99),
        "max": values[-1],
        "mean": statistics.fmean(values),
    }


def _physical_memory_bytes():
    """Return installed physical memory when the platform exposes it."""
    try:
        return int(os.sysconf("SC_PHYS_PAGES")) * int(
            os.sysconf("SC_PAGE_SIZE"))
    except (AttributeError, OSError, ValueError):
        return None


def _hardware_metadata():
    """Describe the current benchmark host without embedding machine identity."""
    architecture = platform.machine() or "unknown"
    processor = platform.processor() or architecture
    return {
        "architecture": architecture,
        "processor": processor,
        "logical_cpu_count": os.cpu_count(),
        "memory_bytes": _physical_memory_bytes(),
        "operating_system": platform.platform(),
        "python": sys.version.split()[0],
        "numpy": np.__version__,
    }


def command_benchmark(args):
    payload = load_checkpoint_payload(args.discovery_checkpoint)
    checkpoint_identity = payload["identity"]
    generators = payload["state"]["generators"]
    if len(generators) != 12:
        raise ValueError("benchmark requires a complete discovery checkpoint")

    schedule = sqlite3.connect(
        f"file:{Path(args.schedule).resolve()}?mode=ro", uri=True)
    schedule_identity = _schedule_identity(schedule)
    catalog_index = load_catalog_index(args.catalog_dir)
    required_identity = {
        "prime": int(args.prime),
        "seed": int(args.seed),
        "catalog_sha256": catalog_index["catalog_sha256"],
        "candidate_count": catalog_index["candidate_count"],
        "schedule_ordering": schedule_identity["ordering"],
        "engine_sha256": _engine_sha256(),
        "discovery_ordering": "hash-bounded",
        "max_greedy_peak": int(args.max_greedy_peak),
    }
    for key, expected in required_identity.items():
        if checkpoint_identity.get(key) != expected:
            schedule.close()
            raise ValueError(
                f"benchmark checkpoint identity mismatch for {key}")
    if (
        schedule_identity["catalog_sha256"] != catalog_index["catalog_sha256"]
        or schedule_identity["candidate_count"]
        != catalog_index["candidate_count"]
    ):
        schedule.close()
        raise ValueError("benchmark schedule/catalog identity mismatch")

    graphs = [graph_from_label(item["graph"]) for item in generators]
    F_dense, basis = _sample_and_basis(args.prime, args.seed)
    sample_rows = schedule.execute(
        """
        SELECT id, upper_triangle
        FROM candidates
        WHERE greedy_peak_work <= ?
        ORDER BY id
        LIMIT ?
        """,
        (args.max_greedy_peak, args.catalog_sample),
    ).fetchall()
    schedule.close()
    sample_graphs = [
        graph_from_record({
            "order": ORDER,
            "upper_triangle": list(upper),
        })
        for _, upper in sample_rows
    ]

    canonical_seconds = []
    for (expected_id, _), M in zip(sample_rows, sample_graphs):
        started = time.perf_counter()
        actual_id = canonical_graph_id(M)
        canonical_seconds.append(time.perf_counter() - started)
        if actual_id != expected_id:
            raise ValueError("canonical benchmark found an ID mismatch")

    greedy_plan_seconds = []
    for M in sample_graphs:
        started = time.perf_counter()
        greedy_contraction_plan_cost(M, D, VALENCE)
        greedy_plan_seconds.append(time.perf_counter() - started)

    optimal_plan_seconds = []
    for M in sample_graphs[:args.optimal_plan_sample]:
        started = time.perf_counter()
        contraction_plan_cost(M, D, VALENCE)
        optimal_plan_seconds.append(time.perf_counter() - started)

    value_seconds = []
    scalar_values = []
    for M in graphs:
        started = time.perf_counter()
        scalar_values.append(
            value(M, F_dense, D, VALENCE, True, args.prime))
        value_seconds.append(time.perf_counter() - started)

    combined_seconds = []
    rows = []
    for expected_scalar, M in zip(scalar_values, graphs):
        started = time.perf_counter()
        scalar, row = value_and_jacobian_row(
            M, F_dense, basis, D, VALENCE, True, args.prime,
            backend="optimized")
        combined_seconds.append(time.perf_counter() - started)
        if scalar != expected_scalar:
            raise ValueError("value backends disagree during benchmark")
        rows.append(row)

    reference_seconds = []
    for expected_scalar, expected_row, M in zip(
            scalar_values, rows, graphs[:args.reference_sample]):
        started = time.perf_counter()
        scalar, row = value_and_jacobian_row(
            M, F_dense, basis, D, VALENCE, True, args.prime,
            backend="reference")
        reference_seconds.append(time.perf_counter() - started)
        if scalar != expected_scalar or not np.array_equal(row, expected_row):
            raise ValueError("reference backend mismatch during benchmark")

    rank_seconds = []
    sieve = RankSieve(basis.ncols, args.prime)
    for row in rows:
        started = time.perf_counter()
        sieve.add(row)
        rank_seconds.append(time.perf_counter() - started)

    checkpoint_seconds = []
    checkpoint_sizes = []
    checkpoint_identity = {"benchmark": 1, "prime": args.prime}
    with tempfile.TemporaryDirectory() as temporary:
        checkpoint_path = Path(temporary) / "checkpoint.json"
        checkpoint_state = {
            "cursor": args.catalog_sample,
            "sieve": sieve.to_state(),
            "generators": generators,
        }
        for _ in range(args.checkpoint_sample):
            started = time.perf_counter()
            write_checkpoint(
                checkpoint_path, checkpoint_identity, checkpoint_state)
            checkpoint_seconds.append(time.perf_counter() - started)
            checkpoint_sizes.append(checkpoint_path.stat().st_size)
        load_checkpoint(checkpoint_path, checkpoint_identity)

    index = load_catalog_index(args.catalog_dir)
    generation_seconds = 0.0
    compressed_bytes = 0
    for entry in index["entries"]:
        shard = Path(args.catalog_dir) / entry["path"]
        manifest = load_manifest(
            shard.with_suffix(shard.suffix + ".manifest.json"))
        generation_seconds += float(manifest["seconds"])
        compressed_bytes += int(manifest["compressed_bytes"])

    benchmark = {
        "schema": 1,
        "hardware": _hardware_metadata(),
        "prime": args.prime,
        "seed": args.seed,
        "catalog": {
            "candidates": index["candidate_count"],
            "shards": index["shards"],
            "logical_sha256": index["catalog_sha256"],
            "generation_enrichment_seconds_sum": generation_seconds,
            "global_verification_seconds": index["verification_seconds"],
            "compressed_bytes": compressed_bytes,
        },
        "stage_seconds": {
            "canonicalization": _distribution(canonical_seconds),
            "greedy_schedule_planning": _distribution(greedy_plan_seconds),
            "globally_optimal_planning": _distribution(
                optimal_plan_seconds),
            "invariant_value": _distribution(value_seconds),
            "optimized_value_and_reverse_jacobian": _distribution(
                combined_seconds),
            "reference_value_and_amputated_jacobian": _distribution(
                reference_seconds),
            "modular_rank_insertion": _distribution(rank_seconds),
            "atomic_checkpoint_write": _distribution(checkpoint_seconds),
        },
        "checkpoint_bytes": _distribution(checkpoint_sizes),
        "peak_rss_bytes": _peak_rss_bytes(),
        "discovery_comparison": {
            "cost_first_before_planner_optimization": {
                "evaluated": 1000,
                "new_degree10": 4,
                "seconds": 225.3,
            },
            "hash_stratified_before_planner_optimization": {
                "evaluated": 132,
                "new_degree10": 12,
                "seconds": 51.48565558299015,
            },
            "hash_stratified_final_engine": {
                "evaluated": payload["state"]["evaluated"],
                "new_degree10": len(generators),
                "seconds": payload["state"]["elapsed_seconds"],
            },
        },
        "optimization_comparison": {
            "globally_optimal_planning_p50_seconds": {
                "before_bitsets_and_cache": 0.09515591600211337,
                "after_bitsets_and_cache": _distribution(
                    optimal_plan_seconds)["p50"],
            },
            "optimized_value_and_jacobian_p50_seconds": {
                "before_bitsets_and_cache": 0.32989416700729635,
                "after_bitsets_and_cache": _distribution(
                    combined_seconds)["p50"],
            },
            "full_two_prime_validation_seconds": {
                "before_bitsets_and_cache": 116.306579,
                "after_bitsets_and_cache": 50.500609,
            },
        },
    }
    atomic_write_json(args.out, benchmark)
    print(f"benchmark written to {args.out}", flush=True)


def build_parser():
    parser = argparse.ArgumentParser()
    subparsers = parser.add_subparsers(dest="command", required=True)
    generate = subparsers.add_parser(
        "generate", help="generate exact checksummed graph shards")
    generate.add_argument("--catalog-dir", default=str(DEFAULT_CATALOG_DIR))
    generate.add_argument("--geng", default="geng")
    generate.add_argument("--multig", default="multig")
    generate.add_argument("--shards", type=int, default=DEFAULT_SHARDS)
    generate.add_argument("--residue", type=int)
    generate.add_argument("--resume", action="store_true")
    generate.set_defaults(function=command_generate)

    schedule = subparsers.add_parser(
        "schedule", help="build a disk-backed global cost schedule")
    schedule.add_argument("--catalog-dir", default=str(DEFAULT_CATALOG_DIR))
    schedule.add_argument("--schedule", default=str(DEFAULT_SCHEDULE))
    schedule.set_defaults(function=command_schedule)

    discover = subparsers.add_parser(
        "discover", help="run or resume exact degree-10 rank discovery")
    discover.add_argument("--catalog-dir", default=str(DEFAULT_CATALOG_DIR))
    discover.add_argument("--schedule", default=str(DEFAULT_SCHEDULE))
    discover.add_argument("--checkpoint", default=str(DEFAULT_CHECKPOINT))
    discover.add_argument("--lower-result", default=str(DEFAULT_LOWER_RESULT))
    discover.add_argument("--prime", type=int, default=P)
    discover.add_argument("--seed", type=int, default=20260729)
    discover.add_argument("--target-new", type=int, default=12)
    discover.add_argument("--max-candidates", type=int)
    discover.add_argument(
        "--schedule-mode", choices=("cost", "hash-bounded"),
        default="hash-bounded")
    discover.add_argument(
        "--max-greedy-peak", type=int, default=100_000_000,
        help="safety bound used by --schedule-mode=hash-bounded")
    discover.add_argument("--batch-size", type=int, default=32)
    discover.add_argument("--checkpoint-every", type=int, default=5)
    discover.add_argument("--progress-every", type=int, default=10)
    discover.set_defaults(function=command_discover)

    validate = subparsers.add_parser(
        "validate", help="validate the selected degree-10 basis exactly")
    validate.add_argument("--catalog-dir", default=str(DEFAULT_CATALOG_DIR))
    validate.add_argument("--schedule", default=str(DEFAULT_SCHEDULE))
    validate.add_argument(
        "--discovery-checkpoint")
    validate.add_argument(
        "--selection-result",
        help="revalidate generators already saved in a committed result")
    validate.add_argument(
        "--skip-catalog-check", action="store_true",
        help="revalidate explicit graphs without a regenerated local catalog")
    validate.add_argument("--lower-result", default=str(DEFAULT_LOWER_RESULT))
    validate.add_argument("--primes", type=int, nargs="+",
                          default=[P, ALT_P])
    validate.add_argument("--jacobian-seeds", type=int, nargs="+",
                          default=[20260729, 20260730, 20260731])
    validate.add_argument("--value-seed-start", type=int, default=20260801)
    validate.add_argument("--value-samples", type=int, default=16)
    validate.add_argument("--out", default=str(DEFAULT_RESULT))
    validate.set_defaults(function=command_validate)

    benchmark = subparsers.add_parser(
        "benchmark", help="measure degree-10 pipeline stages")
    benchmark.add_argument("--catalog-dir", default=str(DEFAULT_CATALOG_DIR))
    benchmark.add_argument("--schedule", default=str(DEFAULT_SCHEDULE))
    benchmark.add_argument(
        "--discovery-checkpoint",
        default="work/degree10-final.checkpoint.json")
    benchmark.add_argument("--prime", type=int, default=P)
    benchmark.add_argument("--seed", type=int, default=20260729)
    benchmark.add_argument("--max-greedy-peak", type=int, default=100_000_000)
    benchmark.add_argument("--catalog-sample", type=int, default=500)
    benchmark.add_argument("--optimal-plan-sample", type=int, default=100)
    benchmark.add_argument("--reference-sample", type=int, default=12)
    benchmark.add_argument("--checkpoint-sample", type=int, default=20)
    benchmark.add_argument("--out", default=str(DEFAULT_BENCHMARKS))
    benchmark.set_defaults(function=command_benchmark)
    return parser


if __name__ == "__main__":
    arguments = build_parser().parse_args()
    arguments.function(arguments)
