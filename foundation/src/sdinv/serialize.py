"""Strict edge-list JSON interchange, retaining both historical formats.

Vertex order is preserved: isomorphism of odd-form contraction graphs may
change the overall sign, so a canonical topology ID is not a formula ID.
"""

import json
from numbers import Integral

import numpy as np

from .graphs import graph_from_label, graph_from_record


def _integer(value, name, minimum=0):
    if isinstance(value, bool) or not isinstance(value, Integral):
        raise ValueError(f"{name} must be an integer")
    if value < minimum or value > np.iinfo(np.int64).max:
        raise ValueError(f"{name} is out of range")
    return int(value)


def to_record(matrix):
    matrix = np.asarray(matrix)
    if matrix.ndim != 2 or matrix.shape[0] != matrix.shape[1] or not len(matrix):
        raise ValueError("adjacency matrix must be nonempty and square")
    if matrix.dtype.kind not in "iu" or np.any(matrix < 0):
        raise ValueError("edge multiplicities must be nonnegative integers")
    if not np.array_equal(matrix, matrix.T) or np.any(np.diag(matrix)):
        raise ValueError("adjacency matrix must be symmetric with zero diagonal")
    n = len(matrix)
    return {"n": n, "edges": [[i, j, _integer(matrix[i, j], "multiplicity", 1)]
            for i in range(n) for j in range(i + 1, n) if matrix[i, j]]}


def from_record(record):
    if not isinstance(record, dict):
        raise ValueError("graph record must be an object")
    if "order" in record and "upper_triangle" in record:
        n = _integer(record["order"], "order", 1)
        for value in record["upper_triangle"]:
            _integer(value, "multiplicity")
        return graph_from_record({"order": n, "upper_triangle": record["upper_triangle"]})
    n = _integer(record.get("n"), "n", 1)
    if not isinstance(record.get("edges"), list):
        raise ValueError("edges must be a list")
    matrix = np.zeros((n, n), dtype=np.int64)
    for edge in record["edges"]:
        if not isinstance(edge, (list, tuple)) or len(edge) != 3:
            raise ValueError("each edge must contain i, j, multiplicity")
        i, j = (_integer(v, "vertex") for v in edge[:2])
        mult = _integer(edge[2], "multiplicity", 1)
        if not 0 <= i < j < n or matrix[i, j]:
            raise ValueError("invalid or repeated edge")
        matrix[i, j] = matrix[j, i] = mult
    return matrix


def loads(data):
    if isinstance(data, str):
        data = data.strip()
        if data.startswith("n"):
            matrix = graph_from_label(data)
            to_record(matrix)
            return matrix
        data = json.loads(data)
    return from_record(data)


def dumps(matrix):
    return json.dumps(to_record(matrix), separators=(",", ":"))
