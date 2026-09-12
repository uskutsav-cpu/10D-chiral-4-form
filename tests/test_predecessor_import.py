import copy

from chiral4form.predecessor_import import (
    GENERATOR_MAP,
    _align_reduced_metadata,
)

def _record(catalogue_marker=False):
    record = {
        "source_prime": 30203,
        "coordinate_ids": ["A"],
        "fields": {
            "tr1": [
                {
                    "nvars": 1,
                    "prime": 30203,
                    "terms": [],
                }
            ]
        },
        "restriction": "certified global six-parameter degree-10 pure-stress orbit",
    }
    if catalogue_marker:
        record["generator_catalogue"] = [
            {"id": "tr1", "factors": ["tr1"], "leading_degree": 4}
        ]
    return record

def test_predecessor_generator_map_complete_and_unique():
    assert len(GENERATOR_MAP) == 18
    assert len(set(GENERATOR_MAP.values())) == 18

def test_metadata_alignment_omits_catalogue_when_native_omits_it():
    native = _record(False)
    imported = _record(True)
    imported["source_prime"] = 32693
    aligned = _align_reduced_metadata(imported, native)
    assert "generator_catalogue" not in aligned
    assert aligned["restriction"] == native["restriction"]
    assert aligned["metadata_alignment"]["reference_prime"] == 30203

def test_metadata_alignment_copies_catalogue_when_native_has_it():
    native = _record(True)
    imported = _record(False)
    imported["source_prime"] = 32693
    aligned = _align_reduced_metadata(imported, native)
    assert aligned["generator_catalogue"] == native["generator_catalogue"]
    assert aligned["generator_catalogue"] is not native["generator_catalogue"]
    assert aligned["restriction"] == native["restriction"]
