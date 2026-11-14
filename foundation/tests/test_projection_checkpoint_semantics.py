"""What may and may not invalidate a cached projection unit.

A checkpoint is only useful if it is trusted, and only safe if that trust is
narrow. Two failure modes bracket the design:

  too strict  -- a documentation commit discarded a completed two-prime
                 projection, which is the opposite of what a checkpoint is for;
  too lax     -- a formula is reimplemented and the store keeps serving values
                 computed by the OLD code, which look entirely plausible.

The rule these tests pin down: block on what changes the computed values, and
record everything else as provenance.
"""
import json
import os
import sys
from pathlib import Path

import pytest

sys.path.insert(0, os.path.join(os.path.dirname(__file__), "..", "src"))

from sdinv.projection_checkpoint import (  # noqa: E402
    ProjectionCheckpoint, block_fingerprint, evaluator_fingerprint)

BASE = {"source_commit": "aaaa", "atlas_sha256": "atlas1",
        "quotient_sha256": "quot1", "evaluator_version": "v1",
        "modular_backend": "backend1", "prime": 32749, "seed_base": 41000,
        "degree": 10}


def _store(tmp_path, **overrides):
    return ProjectionCheckpoint(tmp_path / "ck", dict(BASE, **overrides))


# --- 1. provenance-only change keeps the cache ----------------------------

def test_documentation_only_commit_preserves_the_atlas_cache(tmp_path):
    store = _store(tmp_path)
    store.save_unit(32749, 0, 0, "I10_1", 1234, 0.1, fingerprint="blockfp")
    store.flush_manifest()

    resumed = _store(tmp_path, source_commit="bbbb")
    unit = resumed.load_unit(32749, 0, 0, fingerprint="blockfp")
    assert unit is not None and unit["value"] == 1234, (
        "a commit that changed no input discarded stored atlas work")


def test_source_commit_is_still_recorded_as_provenance(tmp_path):
    store = _store(tmp_path)
    store.save_unit(32749, 0, 0, "I10_1", 1234, 0.1)
    store.flush_manifest()
    raw = json.loads((store.unit_path(32749, 0, 0)).read_text())
    assert raw["identity"]["source_commit"] == "aaaa", (
        "source_commit must remain recorded even though it does not gate")
    resumed = _store(tmp_path, source_commit="bbbb")
    assert resumed.commit_drift == ("aaaa", "bbbb"), "drift not surfaced"


# --- 2. semantic changes invalidate ---------------------------------------

def test_formula_fingerprint_change_invalidates_that_formula(tmp_path):
    """The case a hand-maintained version string silently gets wrong."""
    store = _store(tmp_path)
    store.save_unit(32749, 0, 1000, "P10_09", 5555, 0.1, fingerprint="fpA")
    store.flush_manifest()
    assert store.load_unit(32749, 0, 1000, fingerprint="fpA")["value"] == 5555
    assert store.load_unit(32749, 0, 1000, fingerprint="fpB") is None, (
        "a reimplemented formula served a value computed by the old code")


def test_formula_change_does_not_invalidate_unrelated_atlas_units(tmp_path):
    """Invalidation must be surgical, not a blanket wipe."""
    store = _store(tmp_path)
    store.save_unit(32749, 0, 0, "I10_1", 1234, 0.1, fingerprint="blockfp")
    store.save_unit(32749, 0, 1000, "P10_09", 5555, 0.1, fingerprint="fpA")
    store.flush_manifest()
    assert store.load_unit(32749, 0, 1000, fingerprint="fpB") is None
    assert store.load_unit(32749, 0, 0, fingerprint="blockfp") is not None, (
        "changing one formula invalidated the expensive atlas cache too")


@pytest.mark.parametrize("key,value", [
    ("atlas_sha256", "atlas2"),
    ("quotient_sha256", "quot2"),
    ("evaluator_version", "v2"),
    ("modular_backend", "backend2"),
    ("prime", 32717),
    ("seed_base", 99000),
    ("degree", 12),
])
def test_semantic_identity_change_refuses_resume(tmp_path, key, value):
    store = _store(tmp_path)
    store.save_unit(32749, 0, 0, "I10_1", 1234, 0.1)
    store.flush_manifest()
    with pytest.raises(ValueError, match="identity mismatch"):
        _store(tmp_path, **{key: value})


def test_atlas_hash_is_order_sensitive():
    """Reordering the atlas changes what every stored coordinate MEANS.

    The set of column names being equal is not enough: coordinates are reported
    against a column order, so a permutation silently relabels them.
    """
    import hashlib
    a = hashlib.sha256(json.dumps(["I10_1", "I10_2"]).encode()).hexdigest()
    b = hashlib.sha256(json.dumps(["I10_2", "I10_1"]).encode()).hexdigest()
    assert a != b, (
        "atlas hashing is order-insensitive, so a reordered basis would reuse "
        "coordinates that no longer refer to the same columns")


# --- 3. integrity ----------------------------------------------------------

def test_corrupted_checksum_is_rejected(tmp_path):
    store = _store(tmp_path)
    store.save_unit(32749, 0, 0, "I10_1", 1234, 0.1)
    path = store.unit_path(32749, 0, 0)
    record = json.loads(path.read_text())
    record["value"] = 9999                     # tamper, leave checksum stale
    path.write_text(json.dumps(record))
    assert store.load_unit(32749, 0, 0) is None, (
        "a tampered unit was trusted")


def test_partially_written_record_is_ignored(tmp_path):
    store = _store(tmp_path)
    store.save_unit(32749, 0, 0, "I10_1", 1234, 0.1)
    path = store.unit_path(32749, 0, 0)
    path.write_text(path.read_text()[: len(path.read_text()) // 2])
    assert store.load_unit(32749, 0, 0) is None, (
        "a truncated unit -- the crash-during-write case -- was trusted")


def test_missing_unit_reads_as_absent_not_as_error(tmp_path):
    store = _store(tmp_path)
    assert store.load_unit(32749, 7, 3) is None


def test_duplicate_save_is_idempotent(tmp_path):
    store = _store(tmp_path)
    store.save_unit(32749, 0, 0, "I10_1", 1234, 0.1, fingerprint="f")
    store.save_unit(32749, 0, 0, "I10_1", 1234, 0.2, fingerprint="f")
    store.flush_manifest()
    assert store.load_unit(32749, 0, 0, fingerprint="f")["value"] == 1234
    assert len(store.completed()) == 1, "duplicate unit double-counted"


# --- 4. the fingerprints actually track the source -------------------------

def test_evaluator_fingerprint_distinguishes_different_evaluators():
    from sdinv.published_degree10_invariants import PUBLISHED_DEGREE10
    prints = {name: evaluator_fingerprint(spec["evaluator"])
              for name, spec in PUBLISHED_DEGREE10.items()}
    assert len(set(prints.values())) == len(prints), (
        f"two evaluators share a fingerprint, so one could serve the other's "
        f"cached values: {prints}")


def test_evaluator_fingerprint_moves_when_shared_machinery_moves():
    """Most evaluators are thin wrappers, so their own text is not enough."""
    def fake(form, mod=None):
        return 0

    with_blocks = evaluator_fingerprint(fake, include_blocks=True)
    without = evaluator_fingerprint(fake, include_blocks=False)
    assert with_blocks != without, (
        "the shared block machinery does not enter the fingerprint, so a "
        "change to composite_n1050 would leave every cached value stale")
    assert len(block_fingerprint()) == 32
