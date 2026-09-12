"""Typed access to the frozen Paper-2 baseline."""

from __future__ import annotations

import json
from dataclasses import dataclass
from pathlib import Path


@dataclass(frozen=True)
class DegreeRecord:
    degree: int
    full_dim: int
    static_stress_dim: int
    dynamic_reachable_dim: int
    quotient_dim: int
    new_forcing_dim: int | None = None

    def validate(self) -> None:
        if self.degree <= 0 or self.degree % 2:
            raise ValueError(f"invalid even field degree: {self.degree}")
        dims = [self.full_dim, self.static_stress_dim, self.dynamic_reachable_dim, self.quotient_dim]
        if any(x < 0 for x in dims):
            raise ValueError(f"negative dimension at degree {self.degree}")
        if self.static_stress_dim > self.dynamic_reachable_dim:
            raise ValueError("static stress span cannot exceed declared dynamic closure")
        if self.dynamic_reachable_dim > self.full_dim:
            raise ValueError("reachable space cannot exceed full invariant space")
        if self.full_dim - self.dynamic_reachable_dim != self.quotient_dim:
            raise ValueError("quotient dimension is inconsistent with full - reachable")
        if self.new_forcing_dim is not None and self.new_forcing_dim > self.full_dim:
            raise ValueError("forcing dimension cannot exceed full invariant space")


@dataclass(frozen=True)
class Baseline:
    status: str
    source_repository: str
    source_ref: str
    records: tuple[DegreeRecord, ...]

    def validate(self) -> None:
        if not self.source_repository:
            raise ValueError("missing source repository")
        if not self.source_ref:
            raise ValueError("missing source ref")
        degrees = [r.degree for r in self.records]
        if degrees != sorted(set(degrees)):
            raise ValueError("degree records must be unique and sorted")
        for record in self.records:
            record.validate()


def load_baseline(path: str | Path) -> Baseline:
    payload = json.loads(Path(path).read_text())
    records = tuple(
        DegreeRecord(
            degree=int(r["degree"]),
            full_dim=int(r["full_dim"]),
            static_stress_dim=int(r["static_stress_dim"]),
            dynamic_reachable_dim=int(r["dynamic_reachable_dim"]),
            quotient_dim=int(r["quotient_dim"]),
            new_forcing_dim=(None if r.get("new_forcing_dim") is None else int(r["new_forcing_dim"])),
        )
        for r in payload["degrees"]
    )
    baseline = Baseline(
        status=str(payload["status"]),
        source_repository=str(payload["source_repository"]),
        source_ref=str(payload["source_ref"]),
        records=records,
    )
    baseline.validate()
    return baseline
