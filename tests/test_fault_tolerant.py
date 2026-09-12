from fractions import Fraction

from chiral4form.fault_tolerant import gaussian_rational_candidate, leave_k_out_consensus
from chiral4form.modular_reconstruction import reduce_fraction


def test_gaussian_candidate_for_clean_images():
    primes = [1009, 1013, 1019, 1021]
    q = Fraction(17, 23)
    residues = {p: reduce_fraction(q, p) for p in primes}
    assert gaussian_rational_candidate(residues) == q


def test_leave_one_out_tolerates_one_corrupted_prime():
    primes = [1009, 1013, 1019, 1021, 1031]
    q = Fraction(29, 31)
    residues = {p: reduce_fraction(q, p) for p in primes}
    residues[1019] = (residues[1019] + 7) % 1019
    out = leave_k_out_consensus(residues, max_bad_primes=1)
    assert out["status"] == "verified_support_candidate"
    assert Fraction(out["value"]) == q
    assert 1019 in out["bad_primes"]
