"""Literal degree-10 candidates in arXiv:2509.14350v2, section 4.1.4.

Original TeX labels I101-I1012 and the colored PDF page 25 determine the
brackets unambiguously: candidate 4 symmetrizes four slots, candidate 9 three;
candidate 10 antisymmetrizes a trailing pair, and candidates 11/12 apply a
trailing-triple antisymmetrizer on each of their last three N1050 factors.
Every red operation follows the intrinsic five-index black antisymmetrizer.
All operations are normalized, and each contraction edge carries one metric.

PUBLISHED_DEGREE10 selects these literal source readings. Earlier incomplete
readings are preserved through LEGACY_PUBLISHED_DEGREE10 and
LEGACY_READING_VARIANTS for reproducing archived evidence. They are not alternate unresolved readings of the source.
The paper calls its twelve expressions candidates, without proving that this
particular list is a basis. This module does not assert a span/intersection rank.
"""

import numpy as np

from .index_symmetry_ops import BLACK, RED, BracketOp, BracketProgram
from .modp import P, inv, mod_einsum
from .stress import (
    _raise_axes,
    composite_n1050,
    composite_n4125,
    five_form_moment,
    matrix_trace_power,
)

SOURCE = {
    "paper": "Some remarks on invariants",
    "journal": "J. Phys. A 59 (2026) 065203",
    "doi": "10.1088/1751-8121/ae3bb8",
    "arxiv": "2509.14350v2",
    "equation": "(4.24)",
    "journal_pdf_page": 17,
    "arxiv_pdf_page": 25,
    "transcription_method": "original arXiv v2 TeX and independently rendered colored PDF page 25",
    "source_reading_version": "literal-red-brackets-2026-09-07",
}


def p10_01_trM5(five_form, mod=P, backend="optimized"):
    """I^(1)_10 = tr M^5."""
    _, mixed = five_form_moment(five_form, mod, backend)
    return int(matrix_trace_power(mixed, 5, mod) % mod)


def p10_02_mm_m_m_n4125(five_form, mod=P, backend="optimized"):
    """I^(2)_10 = (MM)_{mu1}^{nu1} M_{mu2}^{nu2} M_{mu3}^{nu3}
                  N^(4125)_{nu1 nu2 nu3}{}^{mu1 mu2 mu3}.

    composite_n4125 returns the all-lower 6-index tensor, so the last three
    indices are raised to match the source index placement.
    """
    _, mixed = five_form_moment(five_form, mod, backend)
    mm = (mixed @ mixed) % mod
    n = composite_n4125(five_form, mod, backend)
    n_raised = _raise_axes(n, (3, 4, 5), mod)
    # a=mu1 d=nu1, b=mu2 e=nu2, c=mu3 f=nu3 ; N indices (nu1,nu2,nu3,mu1,mu2,mu3)
    #
    # mod_einsum, NOT np.einsum. A bare int64 einsum over four operands with
    # entries up to p multiplies before summing: p^4 ~ 1.15e18 accumulated over
    # 10^6 terms wraps past 2^63 silently. That is the documented trap in this
    # repository's README and it was hit here first time; the homogeneity test
    # is what exposed it.
    out = mod_einsum("ad,be,cf,defabc->", [mm, mixed, mixed, n_raised], mod)
    return int(out % mod)


def p10_03_m3_n1050_n1050(five_form, mod=P, backend="optimized"):
    """I^(3)_10 = M_{mu1}^{nu1} M_{mu2}^{nu2} M_{mu3}^{nu3}
                  ( N_(1050)^{[mu1mu2mu3,a1a2]a3} N^(1050)_{[nu1nu2nu3,a1a2]a3} )

    All brackets in this candidate are BLACK, and `composite_n1050` already
    performs the five-index antisymmetrisation [abc,de], so no additional
    bracket program is required. The source is arXiv v2 section 4.1.4,
    candidate 3, original TeX label I103 and colored PDF page 25.
    """
    _, mixed = five_form_moment(five_form, mod, backend)
    n_low = composite_n1050(five_form, mod, backend)
    n_up = _raise_axes(n_low, (0, 1, 2, 3, 4, 5), mod)
    # a=mu1 d=nu1, b=mu2 e=nu2, c=mu3 f=nu3, xyz = a1 a2 a3
    return int(mod_einsum(
        "ad,be,cf,abcxyz,defxyz->",
        [mixed, mixed, mixed, n_up, n_low], mod) % mod)


def p10_06_mm_m_n1050_n4125(five_form, mod=P, backend="optimized"):
    """I^(6)_10 = (MM)_{nu1}^{mu1} M_{nu2}^{mu2}
                  N^(1050)_{[mu1mu2rho1,rho2rho3]rho4}
                  N_(4125)^{rho1rho2rho3,rho4nu1nu2}

    Black brackets only; the five-index antisymmetrisation is already carried
    by composite_n1050; see arXiv v2 original TeX label I106, PDF page 25.
    """
    _, mixed = five_form_moment(five_form, mod, backend)
    mm = (mixed @ mixed) % mod
    n1050 = composite_n1050(five_form, mod, backend)
    n4125 = composite_n4125(five_form, mod, backend)
    n4125_up = _raise_axes(n4125, (0, 1, 2, 3, 4, 5), mod)
    # (MM)_{nu1}^{mu1} -> [n1, m1];  N1050_{m1 m2 r1 r2 r3 r4}
    # N4125^{r1 r2 r3 r4 n1 n2}
    return int(mod_einsum(
        "am,bn,mnpqrs,pqrsab->",
        [mm, mixed, n1050, n4125_up], mod) % mod)


def p10_07_n1050_mm_n1050_n1050(five_form, mod=P, backend="optimized"):
    """I^(7)_10 = N^(1050)_{[rho1rho2rho3,rho4rho5]}{}^{mu} (MM)_{mu}{}^{kappa}
                  ( N_(1050)^{[rho1rho2rho3,}{}_{a1a2]a3}
                    N_(1050)^{[rho4rho5}{}_{kappa,a1a3]a2} )

    Black brackets only. Note the alpha ordering differs between the two inner
    N factors -- a1a2]a3 against a1a3]a2 -- which is carried by the einsum
    labels, not by any additional symmetrisation.

    INDEX PLACEMENT. Every contracted edge must have exactly ONE raised end;
    an edge joining two equally-placed slots contracts with delta instead of
    eta and is not Lorentz invariant. This function originally raised all six
    axes on BOTH inner N factors, so the three alpha edges (a1,a2,a3) joined
    two raised slots. The result was rotation invariant -- delta and eta agree
    on the spatial block -- but NOT boost invariant, and the degree-10
    projection correctly reported it as `not_in_atlas_span` on all six primes
    because a non-scalar cannot lie in the span of genuine invariants.

    The placement that fixes it turns on `mm`. `mixed` is M_{mu}{}^{nu}, so the
    matrix product `mixed @ mixed` contracts the up index of the first with the
    down index of the second, giving (MM)_{mu}{}^{kappa} -- slot 0 DOWN, slot 1
    UP, the opposite of what the original docstring asserted. Hence the outer N
    carries mu UP (axis 5 raised) and the second inner N carries kappa DOWN
    (axis 2 NOT raised).

    Which end of an alpha edge carries the metric is arbitrary: raising
    (0,1,2,3,4,5)/(0,1) instead of (0,1,2)/(0,1,3,4,5) gives bit-identical
    values, which is the check that this placement is right rather than merely
    boost invariant by accident. `test_p10_07_alpha_edge_placement_is_free`
    asserts it.
    """
    _, mixed = five_form_moment(five_form, mod, backend)
    mm = (mixed @ mixed) % mod
    n = composite_n1050(five_form, mod, backend)
    # N_{r1 r2 r3 r4 r5}{}^{m} (MM)_{m}{}^{k} N^{r1 r2 r3}{}_{a1 a2 a3}
    #                                          N^{r4 r5}{}_{k}{}^{a1 a3 a2}
    n_outer = _raise_axes(n, (5,), mod)
    n_a = _raise_axes(n, (0, 1, 2), mod)
    n_b = _raise_axes(n, (0, 1, 3, 4, 5), mod)
    return int(mod_einsum(
        "pqrstm,mk,pqrxyz,stkxzy->",
        [n_outer, mm, n_a, n_b], mod) % mod)


def _outer(left, right, mod):
    """Modular outer product; axes of `left` first, then axes of `right`."""
    return np.tensordot(left % mod, right % mod, axes=0) % mod


def p10_05_mm_m_antisym_n1050_n1050(five_form, mod=P, backend="optimized"):
    """I^(5)_10 = (MM)_{[nu1}{}^{[mu1} M_{nu2]}{}^{mu2]}
                  N^(1050)_{[rho1rho2rho3,rho4 mu1] mu2}
                  N_(1050)^{[rho1rho2rho3,rho4 nu1] nu2}

    Two BLACK antisymmetrisations, both on the M-product indices: one over the
    nu pair (nu1, nu2), one over the mu pair (mu1, mu2). Neither is supplied by
    `composite_n1050` -- that routine antisymmetrises only its own first five
    axes -- so both are applied explicitly through `BracketProgram`.

    Index placement is fixed by consistency, not by the extraction order of the
    stacked scripts (which is unreliable). `mm` and `mixed` are M_{a}{}^{b},
    slot 0 DOWN and slot 1 UP, so slot 1 contracts into the all-lower
    `composite_n1050` and slot 0 into the raised copy. Every contracted edge
    then has exactly one raised end, which the boost test confirms.
    """
    _, mixed = five_form_moment(five_form, mod, backend)
    mm = (mixed @ mixed) % mod
    n_low = composite_n1050(five_form, mod, backend)
    n_up = _raise_axes(n_low, (0, 1, 2, 3, 4, 5), mod)

    # slots: 0 = nu1, 1 = mu1, 2 = nu2, 3 = mu2
    block = _outer(mm, mixed, mod)
    program = BracketProgram(
        ops=[BracketOp("antisym", (0, 2), BLACK, True,
                       "eq (4.24) I^(5)_10 black bracket [nu1 ... nu2]"),
             BracketOp("antisym", (1, 3), BLACK, True,
                       "eq (4.24) I^(5)_10 black bracket [mu1 ... mu2]")],
        source="eq (4.24) I^(5)_10")
    block = program.apply(block, mod)

    # block[nu1, mu1, nu2, mu2] N_{r1r2r3r4 mu1 mu2} N^{r1r2r3r4 nu1 nu2}
    return int(mod_einsum(
        "nmvp,abcdmp,abcdnv->", [block, n_low, n_up], mod) % mod)


def p10_08_n1050_mm_antisym_n1050_n1050(five_form, mod=P, backend="optimized"):
    """I^(8)_10 = N^(1050)_{[rho1rho2rho3,rho4 nu] mu} M^{[nu}{}_{rho5}
                  M^{mu]}{}_{kappa}
                  ( N_(1050)^{[rho1rho2rho3,}{}_{a1a2]a3}
                    N_(1050)^{[rho4 rho5 kappa, a1a3]a2} )

    Structurally I^(7) with the single factor (MM)_{mu}{}^{kappa} split into
    two M factors carrying a BLACK antisymmetrisation over [nu ... mu], and the
    outer N's fifth index renamed rho5 -> nu accordingly. The inner pair is
    identical to I^(7), including the alpha reordering a1a2]a3 against a1a3]a2.

    The raise pattern is derived by the same edge rule that fixed I^(7): the
    outer N carries nu and mu UP (axes 4,5) because both meet the DOWN slot 0
    of an M_{a}{}^{b}; the second inner N carries rho4 UP but rho5 and kappa
    DOWN because those meet the UP slot 1 of an M.
    """
    _, mixed = five_form_moment(five_form, mod, backend)
    n = composite_n1050(five_form, mod, backend)
    outer = _raise_axes(n, (4, 5), mod)
    n_a = _raise_axes(n, (0, 1, 2), mod)
    n_b = _raise_axes(n, (0, 3, 4, 5), mod)

    # slots: 0 = nu, 1 = rho5, 2 = mu, 3 = kappa
    block = _outer(mixed, mixed, mod)
    program = BracketProgram(
        ops=[BracketOp("antisym", (0, 2), BLACK, True,
                       "eq (4.24) I^(8)_10 black bracket [nu ... mu]")],
        source="eq (4.24) I^(8)_10")
    block = program.apply(block, mod)

    return int(mod_einsum(
        "pqrsnm,ntmk,pqrxyz,stkxzy->",
        [outer, block, n_a, n_b], mod) % mod)


def p10_04_mm_m_red_n1050_n1050(five_form, mod=P, backend="optimized",
                                red_reading="all4"):
    """Literal candidate 4: normalized red symmetrization over mu,nu,rho,lambda.

The first-five black operations are present in composite_n1050. The red
operation acts on the four covariant free slots of their contracted pair.
The optional pairs reading is historical and is not the displayed source."""
    _, mixed = five_form_moment(five_form, mod, backend)
    mm = (mixed @ mixed) % mod
    n_low = composite_n1050(five_form, mod, backend)
    n_up = _raise_axes(n_low, (0, 1, 2, 3), mod)

    # T_{mu nu rho lambda}; alpha edges carry exactly one raised end
    tensor = mod_einsum("abcdmn,abcdrl->mnrl", [n_low, n_up], mod)

    if red_reading == "all4":
        program = BracketProgram(
            ops=[BracketOp("sym", (0, 1, 2, 3), RED, True,
                           "eq (4.24) I^(4)_10 red bracket (mu ... rho]lambda)")],
            source="eq (4.24) I^(4)_10 reading all4")
        tensor = program.apply(tensor, mod)
    elif red_reading == "pairs":
        # symmetrise the (mu nu) block against the (rho lambda) block
        swapped = np.transpose(tensor, (2, 3, 0, 1))
        tensor = ((tensor + swapped) % mod
                  * inv(2 % mod, mod)) % mod
    else:
        raise ValueError(f"unknown red_reading {red_reading!r}")

    mm_up = _raise_axes(mm, (0,), mod)
    m_up = _raise_axes(mixed, (0,), mod)
    return int(mod_einsum("mnrl,mn,rl->", [tensor, mm_up, m_up], mod) % mod)


def p10_04_pairs(five_form, mod=P, backend="optimized"):
    """Historical pair-exchange reading of candidate 4; not the literal source."""
    return p10_04_mm_m_red_n1050_n1050(five_form, mod, backend, "pairs")


def p10_09_red_n1050_m_n1050_n1050_n1050(five_form, mod=P, backend="optimized",
                                         red_reading="all3"):
    """Literal candidate 9: normalized red symmetrization over nu,rho,lambda.

The red parenthesis begins after the first five-index bracket closes, so
neither kappa nor mu participates. The all3 default is read directly from
original TeX and colored PDF. rholambda is a historical non-source reading."""
    _, mixed = five_form_moment(five_form, mod, backend)
    n = composite_n1050(five_form, mod, backend)
    n_a = n                                        # all lower
    n_b = _raise_axes(n, (0, 1, 2, 3), mod)
    n_c = _raise_axes(n, (5,), mod)
    n_d = _raise_axes(n, (0, 1, 2, 3, 4, 5), mod)
    m_up = _raise_axes(mixed, (0,), mod)           # M^{kappa mu}

    # slots of T: 0 = nu (down), 1 = mu (up), 2 = rho (down), 3 = lambda (down)
    tensor = mod_einsum("abcdkn,km,abcdrl->nmrl", [n_a, m_up, n_b], mod)
    if red_reading not in ("all3", "rholambda"):
        raise ValueError(f"unknown red_reading {red_reading!r}")
    # Preserve the historical two-slot alternative explicitly. Original TeX
    # resolves the displayed red group as the three slots (nu,rho,lambda).
    slots = (0, 2, 3) if red_reading == "all3" else (2, 3)
    program = BracketProgram(
        ops=[BracketOp("sym", slots, RED, True,
                       "eq (4.24) I^(9)_10 red bracket (nu ... rho]lambda)")],
        source=f"eq (4.24) I^(9)_10 reading {red_reading}")
    tensor = program.apply(tensor, mod)

    return int(mod_einsum(
        "nmrl,efghmn,efghrl->", [tensor, n_c, n_d], mod) % mod)


# The legacy helper is kept for archived readings. Source red operations use
# _source_red on covariant N1050 BEFORE the required axes are raised.
def _source_red(tensor, slots, mod, label):
    return BracketProgram(
        ops=[BracketOp("antisym", slots, RED, True, label)], source=label
    ).apply(tensor, mod)


def _nested_pair(tensor, mod, reading, label):
    """Historical none/pair program; not the current source triple program."""
    if reading == "outer":
        return tensor
    if reading != "nested":
        raise ValueError(f"unknown reading {reading!r}")
    program = BracketProgram(
        ops=[BracketOp("antisym", (4, 5), BLACK, True, label)], source=label)
    return program.apply(tensor, mod)


def p10_10_five_n1050(five_form, mod=P, backend="optimized", reading="source"):
    """Literal candidate 10: red Alt(mu1,mu2) on the first N1050 factor.

The source default applies the pair on all-covariant N1050 before raising.
outer and nested retain the archived none/pair programs, respectively."""
    n = composite_n1050(five_form, mod, backend)
    n1 = (_source_red(n, (4, 5), mod, "I^(10)_10 red [mu1 mu2]")
          if reading == "source" else
          _nested_pair(n, mod, reading, "legacy I^(10)_10 pair"))
    n2 = _raise_axes(n, (0, 1, 2), mod)
    n3 = _raise_axes(n, (0, 3, 4, 5), mod)
    n4 = _raise_axes(n, (4, 5), mod)
    n5 = _raise_axes(n, (0, 1, 2, 3, 4, 5), mod)
    return int(mod_einsum(
        "abcdmn,abcxyz,dpqxzy,efghpq,efghmn->",
        [n1, n2, n3, n4, n5], mod) % mod)


def p10_11_five_n1050(five_form, mod=P, backend="optimized", reading="source"):
    """Literal candidate 11: Q,Q,T,T,T with T=Alt_(last three slots) Q.

Red groups are (nu1,nu2,nu3), (lambda2,lambda3,nu1), and
(lambda1,lambda3,nu3) on factors C,D,E, respectively. All three normalized
red operations occur after intrinsic first-five black antisymmetrization.
outer/nested preserve archived no-red/pair-on-C-only readings."""
    n = composite_n1050(five_form, mod, backend)
    red = (_source_red(n, (3, 4, 5), mod, "I^(11)_10 red trailing triples")
           if reading == "source" else n)
    a = n
    b = _raise_axes(n, (3, 4, 5), mod)
    c = (_raise_axes(red, (0, 1, 2), mod) if reading == "source" else
         _nested_pair(_raise_axes(n, (0, 1, 2), mod), mod, reading,
                      "legacy I^(11)_10 pair on C only"))
    d = _raise_axes(red, (0, 1, 5), mod)
    e = _raise_axes(red, (0, 1, 2, 3, 4, 5), mod)
    return int(mod_einsum(
        "abcxyz,mnoxyz,mnopqr,abstup,cqtsur->",
        [a, b, c, d, e], mod) % mod)


def p10_12_five_n1050(five_form, mod=P, backend="optimized", reading="source"):
    """Literal candidate 12: Q,Q,T,T,T with T=Alt_(last three slots) Q.

Red groups are (nu1,nu2,nu3), (rho2,lambda2,lambda3), and
(nu3,lambda1,lambda3) on factors C,D,E, respectively. Apply them before
raising any of the participating slots. outer/nested are legacy readings."""
    n = composite_n1050(five_form, mod, backend)
    red = (_source_red(n, (3, 4, 5), mod, "I^(12)_10 red trailing triples")
           if reading == "source" else n)
    a = n
    b = _raise_axes(n, (3, 4, 5), mod)
    c = (_raise_axes(red, (0, 1, 2), mod) if reading == "source" else
         _nested_pair(_raise_axes(n, (0, 1, 2), mod), mod, reading,
                      "legacy I^(12)_10 pair on C only"))
    d = _raise_axes(red, (0, 1, 3), mod)
    e = _raise_axes(red, (0, 1, 2, 3, 4, 5), mod)
    return int(mod_einsum(
        "abcxyz,mnoxyz,mnopqr,pasbtu,cqtrsu->",
        [a, b, c, d, e], mod) % mod)


def _nested_variant(fn):
    def wrapped(five_form, mod=P, backend="optimized"):
        return fn(five_form, mod, backend, reading="nested")
    wrapped.__name__ = fn.__name__ + "_nested"
    wrapped.__doc__ = f"{fn.__name__} under the archived AMB-02 nested program."
    return wrapped


p10_10_nested = _nested_variant(p10_10_five_n1050)
p10_11_nested = _nested_variant(p10_11_five_n1050)
p10_12_nested = _nested_variant(p10_12_five_n1050)


PUBLISHED_DEGREE10 = {
    "P10_01": {"source_label": "I^(1)_10", "evaluator": p10_01_trM5,
                "implemented": True, "source_verified": True,
                "brackets": 'black only',
                "source_reading": "literal-red-brackets-2026-09-07"},
    "P10_02": {"source_label": "I^(2)_10", "evaluator": p10_02_mm_m_m_n4125,
                "implemented": True, "source_verified": True,
                "brackets": 'black only',
                "source_reading": "literal-red-brackets-2026-09-07"},
    "P10_03": {"source_label": "I^(3)_10", "evaluator": p10_03_m3_n1050_n1050,
                "implemented": True, "source_verified": True,
                "brackets": 'black only',
                "source_reading": "literal-red-brackets-2026-09-07"},
    "P10_04": {"source_label": "I^(4)_10", "evaluator": p10_04_mm_m_red_n1050_n1050,
                "implemented": True, "source_verified": True,
                "brackets": 'RED sym(mu,nu,rho,lambda)',
                "source_reading": "literal-red-brackets-2026-09-07"},
    "P10_05": {"source_label": "I^(5)_10", "evaluator": p10_05_mm_m_antisym_n1050_n1050,
                "implemented": True, "source_verified": True,
                "brackets": 'black only',
                "source_reading": "literal-red-brackets-2026-09-07"},
    "P10_06": {"source_label": "I^(6)_10", "evaluator": p10_06_mm_m_n1050_n4125,
                "implemented": True, "source_verified": True,
                "brackets": 'black only',
                "source_reading": "literal-red-brackets-2026-09-07"},
    "P10_07": {"source_label": "I^(7)_10", "evaluator": p10_07_n1050_mm_n1050_n1050,
                "implemented": True, "source_verified": True,
                "brackets": 'black only',
                "source_reading": "literal-red-brackets-2026-09-07"},
    "P10_08": {"source_label": "I^(8)_10", "evaluator": p10_08_n1050_mm_antisym_n1050_n1050,
                "implemented": True, "source_verified": True,
                "brackets": 'black only',
                "source_reading": "literal-red-brackets-2026-09-07"},
    "P10_09": {"source_label": "I^(9)_10", "evaluator": p10_09_red_n1050_m_n1050_n1050_n1050,
                "implemented": True, "source_verified": True,
                "brackets": 'RED sym(nu,rho,lambda)',
                "source_reading": "literal-red-brackets-2026-09-07"},
    "P10_10": {"source_label": "I^(10)_10", "evaluator": p10_10_five_n1050,
                "implemented": True, "source_verified": True,
                "brackets": 'RED antisym(mu1,mu2) on first N1050',
                "source_reading": "literal-red-brackets-2026-09-07"},
    "P10_11": {"source_label": "I^(11)_10", "evaluator": p10_11_five_n1050,
                "implemented": True, "source_verified": True,
                "brackets": 'RED antisym(last three slots) on N1050 factors C,D,E',
                "source_reading": "literal-red-brackets-2026-09-07"},
    "P10_12": {"source_label": "I^(12)_10", "evaluator": p10_12_five_n1050,
                "implemented": True, "source_verified": True,
                "brackets": 'RED antisym(last three slots) on N1050 factors C,D,E',
                "source_reading": "literal-red-brackets-2026-09-07"},
}

# Descriptive metadata retained independently of the source-reading selection.
PUBLISHED_DEGREE10['P10_01'].update({'formula': 'tr M^5', 'blocks': {'M': 5}})
PUBLISHED_DEGREE10['P10_02'].update({'formula': '(MM)M M N^(4125)', 'blocks': {'M': 4, 'N4125': 1}})
PUBLISHED_DEGREE10['P10_03'].update({'formula': 'M M M (N^(1050) N^(1050))', 'blocks': {'M': 3, 'N1050': 2}})
PUBLISHED_DEGREE10['P10_06'].update({'formula': '(MM) M N^(1050) N^(4125)', 'blocks': {'M': 3, 'N1050': 1, 'N4125': 1}})
PUBLISHED_DEGREE10['P10_04'].update({'formula': '(MM) M red( N^(1050) N^(1050) )', 'blocks': {'M': 3, 'N1050': 2}})
PUBLISHED_DEGREE10['P10_05'].update({'formula': '(MM) M [black nu, black mu] N^(1050) N^(1050)', 'blocks': {'M': 3, 'N1050': 2}})
PUBLISHED_DEGREE10['P10_08'].update({'formula': 'N^(1050) M M [black] (N^(1050) N^(1050))', 'blocks': {'M': 2, 'N1050': 3}})
PUBLISHED_DEGREE10['P10_09'].update({'formula': 'red( N^(1050) M N^(1050) ) N^(1050) N^(1050)', 'blocks': {'M': 1, 'N1050': 4}})
PUBLISHED_DEGREE10['P10_10'].update({'formula': 'five N^(1050), red trailing pair on first factor', 'blocks': {'N1050': 5}})
PUBLISHED_DEGREE10['P10_11'].update({'formula': 'five N^(1050), red trailing triples on factors C,D,E', 'blocks': {'N1050': 5}})
PUBLISHED_DEGREE10['P10_12'].update({'formula': 'five N^(1050), red trailing triples on factors C,D,E', 'blocks': {'N1050': 5}})
PUBLISHED_DEGREE10['P10_07'].update({'formula': 'N^(1050) (MM) (N^(1050) N^(1050))', 'blocks': {'M': 2, 'N1050': 3}})


NOT_IMPLEMENTED = {}

# Historical reading helpers are kept for reproducing archived evidence.
def p10_09_rholambda(five_form, mod=P, backend="optimized"):
    """Historical two-slot reading of candidate 9; not the literal source."""
    return p10_09_red_n1050_m_n1050_n1050_n1050(
        five_form, mod, backend, red_reading="rholambda")


LEGACY_READING_VARIANTS = {
    "P10_04": ("AMB-01", "pairs", p10_04_pairs),
    "P10_09": ("AMB-01", "rholambda", p10_09_rholambda),
    "P10_10": ("AMB-02", "nested", p10_10_nested),
    "P10_11": ("AMB-02", "nested", p10_11_nested),
    "P10_12": ("AMB-02", "nested", p10_12_nested),
}

# Original arXiv v2 TeX labels I101-I1012 and colored PDF page 25 identify
# the candidates that require red operations after the intrinsic black ones.
BRACKET_STAGES = {
    "P10_03": "black only", "P10_04": "RED present: (mu ... rho]lambda)",
    "P10_05": "black only", "P10_06": "black only",
    "P10_07": "black only", "P10_08": "black only",
    "P10_09": "RED present: (nu ... rho]lambda)",
    "P10_10": "RED present: trailing pair on first N1050",
    "P10_11": "RED present: trailing triples on C,D,E",
    "P10_12": "RED present: trailing triples on C,D,E",
}


def evaluate_implemented(five_form, mod=P, backend="optimized"):
    return {name: spec["evaluator"](five_form, mod, backend)
            for name, spec in PUBLISHED_DEGREE10.items()}


# Explicitly preserve the pre-audit primary readings for archived reproduction.
def _legacy_outer(fn):
    def wrapped(five_form, mod=P, backend="optimized"):
        return fn(five_form, mod, backend, reading="outer")
    wrapped.__name__ = fn.__name__ + "_legacy_outer"
    return wrapped


LEGACY_PUBLISHED_DEGREE10 = {name: dict(spec) for name, spec in PUBLISHED_DEGREE10.items()}
for _name, _fn in (("P10_10", p10_10_five_n1050),
                   ("P10_11", p10_11_five_n1050),
                   ("P10_12", p10_12_five_n1050)):
    LEGACY_PUBLISHED_DEGREE10[_name]["evaluator"] = _legacy_outer(_fn)
    LEGACY_PUBLISHED_DEGREE10[_name]["brackets"] = (
        "legacy outer: intrinsic first-five black only; no red trailing pair"
        if _name == "P10_10" else
        "legacy outer: intrinsic first-five black only; no red trailing triples on C,D,E")
for _spec in LEGACY_PUBLISHED_DEGREE10.values():
    _spec["source_verified"] = False
    _spec["source_reading"] = "legacy-pre-audit"

# Source colors/extent are resolved by the original TeX. Kept empty for callers
# that formerly enumerated unresolved source ambiguities.
AMBIGUITY_VARIANTS = {}
