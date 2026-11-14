"""Durable, identity-checked checkpoints for long exact computations."""

from __future__ import annotations

import hashlib
import json
from pathlib import Path

from .catalog import atomic_write_json


CHECKPOINT_SCHEMA = 1


def _canonical_json(data):
    return json.dumps(
        data, sort_keys=True, separators=(",", ":"), ensure_ascii=True
    ).encode()


def write_checkpoint(path, identity, state):
    """Atomically persist state together with a checksum and run identity."""
    payload = {
        "schema": CHECKPOINT_SCHEMA,
        "identity": identity,
        "state": state,
    }
    envelope = {
        "payload": payload,
        "payload_sha256": hashlib.sha256(_canonical_json(payload)).hexdigest(),
    }
    atomic_write_json(path, envelope)
    return envelope


def load_checkpoint_payload(path, expected_identity=None):
    """Load and verify the complete checkpoint payload."""
    path = Path(path)
    with path.open() as stream:
        envelope = json.load(stream)
    payload = envelope.get("payload")
    if not isinstance(payload, dict):
        raise ValueError(f"missing checkpoint payload in {path}")
    if payload.get("schema") != CHECKPOINT_SCHEMA:
        raise ValueError(f"unsupported checkpoint schema in {path}")
    actual = hashlib.sha256(_canonical_json(payload)).hexdigest()
    if actual != envelope.get("payload_sha256"):
        raise ValueError(f"checkpoint checksum mismatch in {path}")
    if (
        expected_identity is not None
        and payload.get("identity") != expected_identity
    ):
        raise ValueError(f"checkpoint run identity mismatch in {path}")
    return payload


def load_checkpoint(path, expected_identity=None):
    """Load checkpoint state, rejecting corruption or incompatible settings."""
    return load_checkpoint_payload(path, expected_identity)["state"]
