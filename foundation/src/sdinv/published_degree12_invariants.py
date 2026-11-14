"""The three published degree-12 invariants of equation (4.25).

Source of record
----------------
Cederwall, Hutomo, Kuzenko, Lechner, Sorokin, "Some remarks on invariants",
J. Phys. A 59 (2026) 065203, doi:10.1088/1751-8121/ae3bb8.
Cross-checked against arXiv:2509.14350v2.

Equation (4.25), arXiv v2 PDF page 26, journal PDF page 18, transcribed
directly from the rendered source (not from memory):

    I^(1)_12 = tr M^6

    I^(2)_12 = (M^3)^{mu nu} N^(4125)_{mu a1 a2, nu b1 b2} M^{a1 b1} M^{a2 b2}
             = (M^3)^{mu nu} (N^(4125) M M)_{mu nu}

    I^(3)_12 = (N^(4125) M M)^{mu nu} (N^(4125) M M)_{mu nu}

The paper introduces them as invariants that "have appeared in certain models
of the non-linear self-dual 5-form theory".

Implementation notes
--------------------
These are built from the repository's existing, tested primitives rather than
re-derived: `n4125_mm` is exactly (N^(4125) M M)_{mu nu}, and `symmetric_inner`
performs the metric-correct contraction of two symmetric rank-2 tensors. No
convention is re-chosen here; if the repository convention differed from the
source it would be recorded, not silently adjusted.
"""

import numpy as np

from .modp import P, inv
from .stress import (
    _metric,
    five_form_moment,
    matrix_trace_power,
    n4125_mm,
    symmetric_inner,
)

SOURCE = {
    "paper": "Some remarks on invariants",
    "journal": "J. Phys. A 59 (2026) 065203",
    "doi": "10.1088/1751-8121/ae3bb8",
    "arxiv": "2509.14350v2",
    "equation": "(4.25)",
    "arxiv_pdf_page": 26,
    "journal_pdf_page": 18,
}


def _m_lower_and_mixed(five_form, mod, backend):
    return five_form_moment(five_form, mod, backend)


def p12_01_trM6(five_form, mod=P, backend="optimized"):
    """I^(1)_12 = tr M^6."""
    _, mixed = _m_lower_and_mixed(five_form, mod, backend)
    return int(matrix_trace_power(mixed, 6, mod) % mod)


def _m3_lower(five_form, mod, backend):
    """(M^3)_{mu nu} with both indices down."""
    lower, mixed = _m_lower_and_mixed(five_form, mod, backend)
    m3_mixed = (mixed @ mixed % mod) @ mixed % mod       # M^mu_nu cubed
    eta_lower = _metric(mod)[1]
    return (m3_mixed @ eta_lower) % mod


def p12_02_m3_n4125mm(five_form, mod=P, backend="optimized"):
    """I^(2)_12 = (M^3)^{mu nu} (N^(4125) M M)_{mu nu}."""
    m3 = _m3_lower(five_form, mod, backend)
    r = n4125_mm(five_form, mod, backend)
    return int(symmetric_inner(m3, r, mod) % mod)


def p12_03_n4125mm_squared(five_form, mod=P, backend="optimized"):
    """I^(3)_12 = (N^(4125) M M)^{mu nu} (N^(4125) M M)_{mu nu}."""
    r = n4125_mm(five_form, mod, backend)
    return int(symmetric_inner(r, r, mod) % mod)


PUBLISHED_DEGREE12 = {
    "P12_01": {
        "source_label": "I^(1)_12",
        "formula": "tr M^6",
        "evaluator": p12_01_trM6,
        "field_degree": 12,
        "blocks": {"M": 6},
    },
    "P12_02": {
        "source_label": "I^(2)_12",
        "formula": "(M^3)^{mu nu} (N^(4125) M M)_{mu nu}",
        "evaluator": p12_02_m3_n4125mm,
        "field_degree": 12,
        "blocks": {"M": 5, "N4125": 1},
    },
    "P12_03": {
        "source_label": "I^(3)_12",
        "formula": "(N^(4125) M M)^{mu nu} (N^(4125) M M)_{mu nu}",
        "evaluator": p12_03_n4125mm_squared,
        "field_degree": 12,
        "blocks": {"M": 4, "N4125": 2},
    },
}


def evaluate_all(five_form, mod=P, backend="optimized"):
    return {name: spec["evaluator"](five_form, mod, backend)
            for name, spec in PUBLISHED_DEGREE12.items()}
