from fractions import Fraction

from chiral4form.modular_reconstruction import (
    profile_consensus_reconstruct,
    reduce_fraction,
    reconstruct_asymmetric,
)


def residues(q, primes):
    return {p: reduce_fraction(q, p) for p in primes}


def test_asymmetric_reconstruct_denominator_heavy():
    primes = [1000003, 1000033, 1000037, 1000039]
    q = Fraction(7, 100_000_019)
    r = residues(q, primes)
    out = reconstruct_asymmetric(
        r,
        primes[:3],
        numerator_bound=100,
        denominator_bound=200_000_000,
        verification_primes=primes[3:],
    )
    assert out.status == "verified"
    assert Fraction(out.value) == q


def test_profile_consensus_recovers_asymmetric_fraction():
    primes = [1000003, 1000033, 1000037, 1000039]
    q = Fraction(11, 999_983)
    r = residues(q, primes)
    out = profile_consensus_reconstruct(r, primes[:3], verification_primes=primes[3:])
    assert out["status"] == "verified_unique_consensus"
    assert Fraction(out["value"]) == q


def test_sage_style_height_bound():
    from chiral4form.modular_reconstruction import sage_style_echelon_height_bound
    rows = [[Fraction(1), Fraction(1, 3)], [Fraction(0), Fraction(5, 7)]]
    out = sage_style_echelon_height_bound(rows, source_height=10, modulus_product=10**9, ncols=2)
    assert out["proved"] is True
    assert out["common_denominator"] == 21
