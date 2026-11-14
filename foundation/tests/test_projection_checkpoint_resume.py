"""Checkpoint durability: interrupt, resume, and reject corruption.

Built in response to a measured failure mode, not a hypothetical one: the
degree-12 projection's per-column cost grows (6.3 -> 50.8 -> 90.5 s/col) while
RSS climbs past the 1.5 GB ceiling, so an OOM kill late in a prime destroys
about an hour of unrecoverable work.
"""

import json
import os
import sys
import tempfile

import pytest
from pathlib import Path

sys.path.insert(0, os.path.join(os.path.dirname(__file__), "..", "src"))

from sdinv.projection_checkpoint import (  # noqa: E402
    ProjectionCheckpoint, atomic_write_json, peak_rss_mb)

IDENTITY = {"source_commit": "abc123", "atlas_sha256": "deadbeef",
            "basis_sha256": "cafe", "evaluator_version": "1"}


def _store(root, identity=None):
    return ProjectionCheckpoint(root, identity or IDENTITY)


def test_unit_roundtrip_and_completion_tracking():
    with tempfile.TemporaryDirectory() as d:
        cp = _store(d)
        cp.save_unit(32749, 0, 5, "I12_5", 12345, 1.5)
        cp.flush_manifest()
        got = cp.load_unit(32749, 0, 5)
        assert got is not None and got["value"] == 12345
        assert "32749/000/005" in cp.completed()


def test_resume_skips_completed_and_recomputes_missing():
    with tempfile.TemporaryDirectory() as d:
        first = _store(d)
        for col in range(4):
            first.save_unit(32749, 0, col, f"I12_{col}", 100 + col, 0.1)
        first.flush_manifest()

        # simulate a crash: a fresh store re-reads what survived
        second = _store(d)
        assert len(second.completed()) == 4
        for col in range(4):
            assert second.load_unit(32749, 0, col)["value"] == 100 + col
        assert second.load_unit(32749, 0, 4) is None       # never written


def test_corrupt_unit_is_treated_as_absent_not_trusted():
    with tempfile.TemporaryDirectory() as d:
        cp = _store(d)
        cp.save_unit(32749, 0, 1, "I12_1", 777, 0.2)
        path = cp.unit_path(32749, 0, 1)
        record = json.loads(path.read_text())
        record["value"] = 999                    # tamper, leave checksum stale
        path.write_text(json.dumps(record))
        assert cp.load_unit(32749, 0, 1) is None, (
            "a tampered unit must be recomputed, never trusted")


def test_truncated_unit_is_treated_as_absent():
    """The partially-written-at-crash case."""
    with tempfile.TemporaryDirectory() as d:
        cp = _store(d)
        cp.save_unit(32749, 0, 2, "I12_2", 42, 0.1)
        path = cp.unit_path(32749, 0, 2)
        path.write_text(path.read_text()[: len(path.read_text()) // 2])
        assert cp.load_unit(32749, 0, 2) is None


def test_atlas_hash_mismatch_refuses_resume():
    with tempfile.TemporaryDirectory() as d:
        cp = _store(d)
        cp.save_unit(32749, 0, 0, "I12_0", 1, 0.1)
        cp.flush_manifest()
        try:
            _store(d, {**IDENTITY, "atlas_sha256": "different"})
        except ValueError as exc:
            assert "identity mismatch" in str(exc)
            return
        raise AssertionError("resume against a different atlas was allowed")


def test_merge_order_does_not_affect_result():
    with tempfile.TemporaryDirectory() as d1, tempfile.TemporaryDirectory() as d2:
        a = _store(d1)
        for col in (0, 1, 2):
            a.save_unit(32749, 0, col, f"I12_{col}", col * 7, 0.1)
        a.flush_manifest()
        b = _store(d2)
        for col in (2, 0, 1):
            b.save_unit(32749, 0, col, f"I12_{col}", col * 7, 0.1)
        b.flush_manifest()
        assert a.completed() == b.completed()
        for col in (0, 1, 2):
            assert a.load_unit(32749, 0, col)["value"] == \
                   b.load_unit(32749, 0, col)["value"]


def test_atomic_write_leaves_no_temp_files():
    with tempfile.TemporaryDirectory() as d:
        target = Path(d) / "manifest.json"
        atomic_write_json(target, {"a": 1})
        assert target.exists()
        assert not [p for p in Path(d).iterdir() if ".tmp." in p.name]


def test_peak_rss_is_reported_in_megabytes():
    value = peak_rss_mb()
    assert 1.0 < value < 100000.0, f"implausible peak RSS {value} MB"


def test_commit_drift_does_not_block_resume_but_atlas_change_does(tmp_path):
    """A new commit must not discard valid work; a new atlas must.

    The first version treated `source_commit` as invalidating, so committing a
    documentation change threw away a completed two-prime projection. Only
    fields that actually change the computed values may block a resume.
    """
    base = {"source_commit": "aaaa", "atlas_sha256": "atlas1",
            "evaluator_version": "v1", "prime": 32749, "degree": 10}
    store = ProjectionCheckpoint(tmp_path / "ck", base)
    store.save_unit(32749, 0, 0, "I10_1", 1234, 0.1)
    store.flush_manifest()

    moved = dict(base, source_commit="bbbb")
    resumed = ProjectionCheckpoint(tmp_path / "ck", moved)
    unit = resumed.load_unit(32749, 0, 0)
    assert unit is not None and unit["value"] == 1234, (
        "a commit that changes no input discarded stored units")
    assert resumed.commit_drift == ("aaaa", "bbbb"), "drift not recorded"

    for key, value in (("atlas_sha256", "atlas2"),
                       ("evaluator_version", "v2")):
        with pytest.raises(ValueError, match="identity mismatch"):
            ProjectionCheckpoint(tmp_path / "ck", dict(base, **{key: value}))
