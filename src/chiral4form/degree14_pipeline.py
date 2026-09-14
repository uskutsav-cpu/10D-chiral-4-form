"""Degree-14 falsification pipeline schema.

No degree-14 basis is invented. The pipeline becomes executable only after
an explicit degree-14 basis/map is installed into the registry.
"""
from __future__ import annotations

from .obstruction_ideal import degree14_probe_readiness


def degree14_project_plan(registry):
    readiness = degree14_probe_readiness(registry)
    return {
        "schema": 1,
        "status": readiness["status"],
        "readiness": readiness,
        "pipeline": [
            "install explicit degree-14 invariant basis with provenance",
            "compute/fit the pure-stress degree-14 coefficient map",
            "certify the characteristic-zero map (not modular fit only)",
            "derive degree-14 obstruction candidates",
            "prolong the lifted lower obstruction ideal to the degree-14 ring",
            "reduce every degree-14 obstruction candidate modulo that lower ideal",
            "record the quotient dimension/new generator count",
        ],
        "headline_question": (
            "Does degree 14 contain any genuinely new obstruction generator "
            "modulo prolongations of the lower all-orders ideal?"
        ),
        "success_signal": "new_generator_count = 0",
        "claim_boundary": (
            "A zero degree-14 quotient is evidence of stabilization, not an "
            "all-orders finite-generation proof."
        ),
    }
