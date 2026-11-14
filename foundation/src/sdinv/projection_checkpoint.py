"""Durable, resumable checkpointing for atlas-column projections.

Motivation, from a measured failure mode rather than a hypothetical: the
degree-12 projection evaluates 72 atlas columns x 80 samples per prime, and
its cost per column *grows* as the columns get structurally larger --

    cols  1-12:   76 s   ( 6.3 s/col)
    cols 13-24:  609 s   (50.8 s/col)
    cols 25-36: 1086 s   (90.5 s/col)

with resident memory climbing past the 1.5 GB working ceiling on an 8 GB
machine. A single OOM kill at column 60 therefore destroys roughly an hour of
work that nothing can recover, because the original script held everything in
memory and wrote once at the end.

Design
------
One immutable file per completed unit, plus a small manifest. Units are never
rewritten, so a crash can corrupt at most the unit being written, and that
unit simply fails its checksum on resume and is recomputed.

    checkpoints/p12/
        manifest.json
        prime_32749/sample_000/column_000.json
                               column_001.json

Every unit records the provenance needed to reject a stale resume: source
commit, atlas hash, basis hash, evaluator version. A resume against a
different atlas is refused rather than silently mixed.

Manifest writes are atomic: temp file in the same directory, flush, fsync,
os.replace, then fsync the parent directory. os.replace is atomic on POSIX,
so a reader never observes a half-written manifest.
"""

import hashlib
import json
import os
import resource
import sys
import time
from pathlib import Path

CHECKPOINT_SCHEMA = 1


def _checksum(payload):
    """Stable checksum over the unit's semantic content."""
    body = json.dumps(payload, sort_keys=True, separators=(",", ":"))
    return hashlib.sha256(body.encode()).hexdigest()


def atomic_write_json(path, payload):
    """Write JSON atomically: temp + flush + fsync + replace + dir fsync."""
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    tmp = path.with_suffix(path.suffix + f".tmp.{os.getpid()}")
    with tmp.open("w") as handle:
        json.dump(payload, handle, sort_keys=True, indent=1)
        handle.flush()
        os.fsync(handle.fileno())
    os.replace(tmp, path)
    try:
        fd = os.open(str(path.parent), os.O_RDONLY)
        try:
            os.fsync(fd)
        finally:
            os.close(fd)
    except OSError:
        pass          # directory fsync unsupported on some filesystems


def _source_of(obj):
    try:
        import inspect
        return inspect.getsource(obj)
    except (OSError, TypeError):
        return f"<unavailable:{getattr(obj, '__name__', obj)!r}>"


def block_fingerprint():
    """Fingerprint of the shared tensor machinery every evaluator depends on.

    An evaluator's own source is not sufficient: almost all of them are thin
    wrappers over `composite_n1050`, `_raise_axes`, `mod_einsum` and the
    bracket engine, so a change THERE changes every cached value while leaving
    each evaluator's own text untouched. Hashing the shared modules closes
    that hole. It over-invalidates on comment edits, which is the right way to
    be wrong.
    """
    from . import index_symmetry_ops, modp, stress
    body = "".join(_source_of(m) for m in (stress, modp, index_symmetry_ops))
    return hashlib.sha256(body.encode()).hexdigest()[:32]


def evaluator_fingerprint(evaluator, include_blocks=True):
    """Semantic key for one evaluator: its own source plus the shared blocks.

    This replaces a hand-maintained version string. A hand-maintained string
    only invalidates when someone remembers to bump it, and the case where it
    is forgotten -- a formula edited without a bump -- silently serves stale
    values that look entirely plausible. Deriving the key from source removes
    the remembering.
    """
    body = _source_of(evaluator)
    if include_blocks:
        body += block_fingerprint()
    return hashlib.sha256(body.encode()).hexdigest()[:32]


def peak_rss_mb():
    """Peak resident set size in MiB.

    The unit of ru_maxrss is platform-dependent: bytes on Darwin/BSD,
    kibibytes on Linux. Deciding by MAGNITUDE is wrong -- a genuine ~1 GB
    peak on macOS reads as 1.03e9 bytes, which sits just under a 1<<30
    threshold and is then divided as if it were KiB, reporting ~1e6 MB.
    That is exactly how this function first failed its own test. Decide by
    platform instead, where the answer is unambiguous.
    """
    usage = resource.getrusage(resource.RUSAGE_SELF).ru_maxrss
    if sys.platform == "darwin":
        return usage / (1024 * 1024)          # bytes -> MiB
    return usage / 1024                        # kibibytes -> MiB


def current_rss_mb():
    try:
        with open(f"/proc/{os.getpid()}/statm") as handle:
            return int(handle.read().split()[1]) * os.sysconf("SC_PAGE_SIZE") / 1048576
    except OSError:
        import subprocess
        out = subprocess.run(["ps", "-o", "rss=", "-p", str(os.getpid())],
                             capture_output=True, text=True).stdout.strip()
        return float(out) / 1024 if out else 0.0


class ProjectionCheckpoint:
    """Column-level checkpoint store for one projection run."""

    def __init__(self, root, identity, max_rss_mb=1500):
        self.root = Path(root)
        self.identity = dict(identity)
        self.max_rss_mb = max_rss_mb
        self.root.mkdir(parents=True, exist_ok=True)
        self.manifest_path = self.root / "manifest.json"
        self.commit_drift = None
        self.manifest = self._load_manifest()

    def _load_manifest(self):
        if not self.manifest_path.exists():
            return {"schema": CHECKPOINT_SCHEMA, "identity": self.identity,
                    "units": {}, "created": time.time()}
        with self.manifest_path.open() as handle:
            manifest = json.load(handle)
        stored = manifest.get("identity", {})
        # Only fields that actually INVALIDATE stored values may block a
        # resume. `source_commit` is provenance, not invalidation: committing a
        # documentation change, or adding a new evaluator alongside the
        # existing ones, does not alter any value already computed. Treating it
        # as blocking discarded a completed two-prime projection the first time
        # this ran, which is the opposite of what a checkpoint is for.
        #
        # `evaluator_version` is the field that exists to invalidate, and it
        # must be bumped whenever an existing evaluator's OUTPUT changes.
        # `block_fingerprint` is deliberately NOT here. It is enforced per-unit
        # instead: a change to the shared tensor machinery should invalidate
        # the units it actually affects, not detonate the whole store, and the
        # per-unit key already does that precisely.
        for key in ("atlas_sha256", "basis_sha256", "quotient_sha256",
                    "evaluator_version", "modular_backend",
                    "prime", "seed_base", "degree"):
            if key in self.identity and key in stored:
                if stored[key] != self.identity[key]:
                    raise ValueError(
                        f"checkpoint identity mismatch on {key!r}: stored "
                        f"{stored[key]!r} != current {self.identity[key]!r}. "
                        f"Refusing to resume against different inputs.")
        drift = (stored.get("source_commit"), self.identity.get("source_commit"))
        if drift[0] and drift[1] and drift[0] != drift[1]:
            self.commit_drift = drift          # recorded, not fatal
        return manifest

    def unit_path(self, prime, sample, column):
        return (self.root / f"prime_{prime}" / f"sample_{sample:03d}"
                / f"column_{column:03d}.json")

    def key(self, prime, sample, column):
        return f"{prime}/{sample:03d}/{column:03d}"

    def load_unit(self, prime, sample, column, fingerprint=None):
        """Return a completed unit's value, or None if it cannot be trusted.

        A unit is rejected -- and therefore recomputed -- when it is absent,
        unreadable, fails its checksum (the partially-written-at-crash case),
        or was produced by different code than the caller is running now.

        The `fingerprint` is the per-unit semantic key. Global identity alone
        is too coarse: it cannot distinguish "the atlas is unchanged but this
        one formula was reimplemented", which is the common case while a
        classification is being developed. Passing the evaluator's source
        fingerprint makes a formula change invalidate exactly that formula's
        cached values and nothing else.
        """
        path = self.unit_path(prime, sample, column)
        if not path.exists():
            return None
        try:
            with path.open() as handle:
                record = json.load(handle)
        except (json.JSONDecodeError, OSError):
            return None
        stored = record.pop("checksum", None)
        if stored is None or _checksum(record) != stored:
            return None
        if fingerprint is not None and record.get("fingerprint") != fingerprint:
            return None
        return record

    def save_unit(self, prime, sample, column, invariant_id, value, seconds,
                  fingerprint=None):
        record = {
            "schema": CHECKPOINT_SCHEMA,
            "prime": int(prime), "sample": int(sample),
            "column": int(column), "invariant_id": invariant_id,
            "value": int(value), "seconds": round(float(seconds), 3),
            "identity": self.identity,
            "fingerprint": fingerprint,
            "peak_rss_mb": round(peak_rss_mb(), 1),
        }
        record["checksum"] = _checksum(record)
        atomic_write_json(self.unit_path(prime, sample, column), record)
        self.manifest["units"][self.key(prime, sample, column)] = {
            "invariant_id": invariant_id, "seconds": record["seconds"]}

    def flush_manifest(self):
        self.manifest["updated"] = time.time()
        atomic_write_json(self.manifest_path, self.manifest)

    def completed(self):
        return set(self.manifest.get("units", {}))

    def over_memory_limit(self):
        return current_rss_mb() > self.max_rss_mb
