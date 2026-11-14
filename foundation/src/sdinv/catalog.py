"""Streaming, checksummed contraction-graph catalogs.

Large degree catalogs must never be accumulated in one Python list or one
giant JSON document. A catalog is a set of deterministic gzip JSONL shards.
Each record contains an exact canonical ID plus the unambiguous upper triangle
of its adjacency matrix. Each shard has an atomic manifest over the logical
uncompressed record stream.
"""

from __future__ import annotations

import gzip
import hashlib
import json
import os
from pathlib import Path
import tempfile
import time

from .graphs import canonical, graph_from_record, graph_label, graph_to_record


CATALOG_SCHEMA = 1


def atomic_write_json(path, data):
    """Durably replace path with JSON data."""
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    fd, temporary = tempfile.mkstemp(
        prefix=f".{path.name}.", suffix=".tmp", dir=path.parent)
    try:
        with os.fdopen(fd, "w") as stream:
            json.dump(data, stream, indent=2, sort_keys=True)
            stream.write("\n")
            stream.flush()
            os.fsync(stream.fileno())
        os.replace(temporary, path)
        directory_fd = os.open(path.parent, os.O_RDONLY)
        try:
            os.fsync(directory_fd)
        finally:
            os.close(directory_fd)
    except BaseException:
        try:
            os.unlink(temporary)
        except FileNotFoundError:
            pass
        raise


def canonical_graph_id(M):
    certificate = canonical(M)
    if not isinstance(certificate, bytes):
        certificate = (
            "exact-tuple-v1:" + ",".join(str(int(x)) for x in certificate)
        ).encode()
    return hashlib.sha256(certificate).hexdigest()


def candidate_record(M, **extra):
    record = graph_to_record(M)
    record.update({
        "id": canonical_graph_id(M),
        "label": graph_label(M),
    })
    record.update(extra)
    return record


def default_manifest_path(shard_path):
    shard_path = Path(shard_path)
    return shard_path.with_suffix(shard_path.suffix + ".manifest.json")


def write_graph_shard(shard_path, graphs, generation, manifest_path=None,
                      enrich=None):
    """Stream graphs to one atomic gzip JSONL shard and write its manifest."""
    shard_path = Path(shard_path)
    manifest_path = (
        Path(manifest_path) if manifest_path
        else default_manifest_path(shard_path))
    shard_path.parent.mkdir(parents=True, exist_ok=True)
    fd, temporary = tempfile.mkstemp(
        prefix=f".{shard_path.name}.", suffix=".tmp", dir=shard_path.parent)
    os.close(fd)

    digest = hashlib.sha256()
    count = 0
    first_id = last_id = None
    started = time.perf_counter()
    try:
        with open(temporary, "wb") as raw:
            with gzip.GzipFile(
                    filename="", mode="wb", fileobj=raw, mtime=0) as stream:
                for M in graphs:
                    extra = enrich(M) if enrich else {}
                    if set(extra) & {"id", "label", "order",
                                     "upper_triangle", "shard_index"}:
                        raise ValueError("catalog enrichment uses a reserved key")
                    record = candidate_record(
                        M, shard_index=count, **extra)
                    encoded = (
                        json.dumps(record, sort_keys=True, separators=(",", ":"))
                        + "\n"
                    ).encode()
                    stream.write(encoded)
                    digest.update(encoded)
                    count += 1
                    first_id = first_id or record["id"]
                    last_id = record["id"]
            raw.flush()
            os.fsync(raw.fileno())
        os.replace(temporary, shard_path)
    except BaseException:
        try:
            os.unlink(temporary)
        except FileNotFoundError:
            pass
        raise

    manifest = {
        "schema": CATALOG_SCHEMA,
        "format": "gzip-jsonl",
        "path": shard_path.name,
        "count": count,
        "logical_sha256": digest.hexdigest(),
        "compressed_bytes": shard_path.stat().st_size,
        "first_id": first_id,
        "last_id": last_id,
        "generation": generation,
        "seconds": round(time.perf_counter() - started, 6),
    }
    atomic_write_json(manifest_path, manifest)
    return manifest


def load_manifest(path):
    with open(path) as stream:
        manifest = json.load(stream)
    if manifest.get("schema") != CATALOG_SCHEMA:
        raise ValueError(f"unsupported catalog schema in {path}")
    return manifest


def iter_graph_shard(shard_path, manifest_path=None, verify=True):
    """Yield (record, matrix), optionally verifying count and logical hash."""
    shard_path = Path(shard_path)
    manifest_path = (
        Path(manifest_path) if manifest_path
        else default_manifest_path(shard_path))
    manifest = load_manifest(manifest_path) if verify else None
    digest = hashlib.sha256()
    count = 0
    with gzip.open(shard_path, "rb") as stream:
        for encoded in stream:
            if not encoded.strip():
                continue
            digest.update(encoded)
            record = json.loads(encoded)
            M = graph_from_record(record)
            if verify and canonical_graph_id(M) != record["id"]:
                raise ValueError(f"canonical ID mismatch in {shard_path}")
            count += 1
            yield record, M
    if verify:
        if count != manifest["count"]:
            raise ValueError(f"record count mismatch in {shard_path}")
        if digest.hexdigest() != manifest["logical_sha256"]:
            raise ValueError(f"logical SHA-256 mismatch in {shard_path}")


def verify_graph_shard(shard_path, manifest_path=None):
    count = sum(
        1 for _ in iter_graph_shard(
            shard_path, manifest_path=manifest_path, verify=True))
    return count
