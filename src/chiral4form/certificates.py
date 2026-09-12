"""Certificate helpers for quotient and obstruction claims."""

from __future__ import annotations

from dataclasses import dataclass
from typing import Sequence

from .finite_field import annihilator_basis, matvec, matrix_rank, normalize_projective


@dataclass(frozen=True)
class QuotientCertificate:
    prime: int
    ambient_dim: int
    reachable_rank: int
    quotient_dim: int
    annihilators: tuple[tuple[int, ...], ...]


def build_quotient_certificate(
    reachable_rows: Sequence[Sequence[int]], ambient_dim: int, prime: int
) -> QuotientCertificate:
    if ambient_dim <= 0:
        raise ValueError("ambient_dim must be positive")
    if any(len(row) != ambient_dim for row in reachable_rows):
        raise ValueError("reachable row has wrong ambient dimension")
    rank = matrix_rank(reachable_rows, prime) if reachable_rows else 0
    anns = annihilator_basis(reachable_rows, prime) if reachable_rows else [
        [1 if i == j else 0 for i in range(ambient_dim)] for j in range(ambient_dim)
    ]
    normalized = tuple(tuple(normalize_projective(v, prime)) for v in anns)
    cert = QuotientCertificate(
        prime=prime,
        ambient_dim=ambient_dim,
        reachable_rank=rank,
        quotient_dim=ambient_dim - rank,
        annihilators=normalized,
    )
    verify_quotient_certificate(reachable_rows, cert)
    return cert


def verify_quotient_certificate(
    reachable_rows: Sequence[Sequence[int]], cert: QuotientCertificate
) -> None:
    if cert.quotient_dim != cert.ambient_dim - cert.reachable_rank:
        raise ValueError("certificate quotient dimension mismatch")
    if len(cert.annihilators) != cert.quotient_dim:
        raise ValueError("wrong number of annihilators")
    for ell in cert.annihilators:
        if len(ell) != cert.ambient_dim:
            raise ValueError("annihilator dimension mismatch")
        if any(matvec(reachable_rows, ell, cert.prime)):
            raise ValueError("annihilator does not vanish on reachable rows")
