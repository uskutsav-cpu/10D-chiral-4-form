"""Certificates certify the supplied matrix, not its physical interpretation."""
from __future__ import annotations
from dataclasses import dataclass
from typing import Sequence
from .finite_field import (_check_prime_field_modulus, annihilator_basis, determinant,
                           independent_rows, matvec, matrix_rank, rref)

@dataclass(frozen=True)
class QuotientCertificate:
    prime: int
    ambient_dim: int
    reachable_rank: int  # legacy name; rank of the supplied row matrix only
    quotient_dim: int
    annihilators: tuple[tuple[int,...],...]
    pivot_rows: tuple[int,...] = ()
    pivot_columns: tuple[int,...] = ()
    minor_determinant: int = 1
    interpretation: str = "supplied_row_span_only"

    def verify(self, rows) -> bool:
        """Boolean convenience interface; the free verifier raises on failure."""
        try:
            verify_quotient_certificate(rows, self)
            return True
        except (ValueError, TypeError, IndexError):
            return False


def build_quotient_certificate(reachable_rows: Sequence[Sequence[int]], ambient_dim: int,
                               prime: int) -> QuotientCertificate:
    _check_prime_field_modulus(prime)
    if ambient_dim < 1 or any(len(row)!=ambient_dim for row in reachable_rows):
        raise ValueError("invalid ambient dimension or malformed row")
    rr,piv = rref(reachable_rows,prime,ncols=ambient_dim)
    pr = independent_rows(reachable_rows,prime)
    minor = [[reachable_rows[i][j] for j in piv] for i in pr]
    cert = QuotientCertificate(prime,ambient_dim,len(piv),ambient_dim-len(piv),
        tuple(map(tuple,annihilator_basis(reachable_rows,prime,ncols=ambient_dim))),
        tuple(pr),tuple(piv),determinant(minor,prime))
    verify_quotient_certificate(reachable_rows,cert)
    return cert


def verify_quotient_certificate(reachable_rows: Sequence[Sequence[int]], cert: QuotientCertificate) -> None:
    p,n = cert.prime,cert.ambient_dim
    _check_prime_field_modulus(p)
    if n < 1 or any(len(row)!=n for row in reachable_rows):
        raise ValueError("malformed row matrix")
    if not 0 <= cert.reachable_rank <= n:
        raise ValueError("invalid rank")
    if cert.quotient_dim != n-cert.reachable_rank:
        raise ValueError("incorrect codimension")
    if any(len(row)!=n for row in cert.annihilators):
        raise ValueError("malformed annihilator")
    if len(cert.annihilators)!=cert.quotient_dim or matrix_rank(cert.annihilators,p)!=cert.quotient_dim:
        raise ValueError("annihilators are not a basis")
    if matrix_rank(reachable_rows,p)!=cert.reachable_rank:
        raise ValueError("claimed rank does not match the supplied rows")
    for ell in cert.annihilators:
        if any(matvec(reachable_rows,ell,p)):
            raise ValueError("functional does not annihilate row space")
    if len(cert.pivot_rows)!=cert.reachable_rank or len(cert.pivot_columns)!=cert.reachable_rank:
        raise ValueError("missing rank witness")
    if any(i<0 or i>=len(reachable_rows) for i in cert.pivot_rows) or any(j<0 or j>=n for j in cert.pivot_columns):
        raise ValueError("invalid pivot index")
    minor = [[reachable_rows[i][j] for j in cert.pivot_columns] for i in cert.pivot_rows]
    det = determinant(minor,p)
    if det == 0 or det != cert.minor_determinant:
        raise ValueError("nonvanishing minor certificate failed")
