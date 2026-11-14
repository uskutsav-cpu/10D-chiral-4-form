"""Versioned registry for the verified low-degree five-form basis.

The registry consumes the committed graph representatives verbatim.  It does
not replace them with aesthetically simpler contractions unless an exact
change of basis is computed elsewhere.

Degree 12 is represented by an explicit extension interface.  The ten
product directions are already determined by the verified lower-degree
generators, while the 62 primitive slots intentionally have no formula until
the separate degree-12 result is imported.
"""

from dataclasses import dataclass
import json
from pathlib import Path

import numpy as np

from .contract import value
from .graphs import (
    graph_from_label,
    graph_from_record,
    graph_to_record,
    validate_graph,
)
from .interaction import invariant_value_and_derivative


@dataclass(frozen=True)
class InvariantItem:
    """One homogeneous basis element."""

    id: str
    degree: int
    kind: str
    graph: str | None = None
    graph_record: dict | None = None
    factors: tuple[str, ...] = ()
    source: str | None = None

    def __post_init__(self):
        if self.kind == "graph" and not self.graph:
            raise ValueError(f"graph item {self.id!r} has no graph formula")
        if self.kind == "graph_record":
            if not self.graph_record:
                raise ValueError(
                    f"graph-record item {self.id!r} has no graph record")
            matrix = graph_from_record(self.graph_record)
            if matrix.shape[0] != self.degree:
                raise ValueError(
                    f"{self.id}: graph order and degree do not agree")
            validate_graph(matrix, valence=5, max_mult=4)
        if self.kind == "product" and not self.factors:
            raise ValueError(f"product item {self.id!r} has no factors")
        if self.kind not in {
            "graph", "graph_record", "product", "placeholder"
        }:
            raise ValueError(f"unknown invariant kind: {self.kind!r}")


class InvariantRegistry:
    """Exact homogeneous bases and recursive value evaluation."""

    def __init__(self, degrees, metadata=None):
        self._degrees = {
            int(degree): tuple(items)
            for degree, items in degrees.items()
        }
        self.metadata = dict(metadata or {})
        all_items = [
            item for items in self._degrees.values() for item in items]
        self._by_id = {item.id: item for item in all_items}
        if len(self._by_id) != len(all_items):
            raise ValueError("invariant IDs must be globally unique")
        for degree, items in self._degrees.items():
            if any(item.degree != degree for item in items):
                raise ValueError("registry degree does not match item degree")
            for item in items:
                if item.kind == "product":
                    try:
                        factor_degree = sum(
                            self._by_id[factor].degree
                            for factor in item.factors)
                    except KeyError as exc:
                        raise ValueError(
                            f"unknown product factor: {exc.args[0]!r}") from exc
                    if factor_degree != degree:
                        raise ValueError(
                            f"{item.id}: product degree {factor_degree}, "
                            f"expected {degree}")

    @property
    def degrees(self):
        return tuple(sorted(self._degrees))

    def basis(self, degree):
        try:
            return self._degrees[int(degree)]
        except KeyError as exc:
            raise KeyError(f"degree {degree} is not registered") from exc

    def item(self, item_id):
        return self._by_id[item_id]

    def evaluate_item(self, item_id, five_form, mod, cache=None):
        cache = {} if cache is None else cache
        if item_id in cache:
            return cache[item_id]
        item = self.item(item_id)
        if item.kind == "placeholder":
            raise RuntimeError(
                f"{item.id} is a degree-{item.degree} import placeholder")
        if item.kind == "graph":
            result = value(
                graph_from_label(item.graph),
                five_form,
                10,
                5,
                True,
                mod,
            )
        elif item.kind == "graph_record":
            result = value(
                graph_from_record(item.graph_record),
                five_form,
                10,
                5,
                True,
                mod,
            )
        else:
            result = 1
            for factor in item.factors:
                result = (
                    result
                    * self.evaluate_item(factor, five_form, mod, cache)
                ) % mod
        cache[item_id] = int(result) % mod
        return cache[item_id]

    def evaluate_degree(self, degree, five_form, mod):
        cache = {}
        return [
            self.evaluate_item(item.id, five_form, mod, cache)
            for item in self.basis(degree)
        ]

    def evaluate_item_with_gradient(self, item_id, five_form, mod, cache=None):
        """Evaluate an invariant and its paper-normalized form derivative.

        The derivative is the covariant antisymmetric tensor
        ``dI/dLambda^{mu(5)}``.  Product items use the exact Leibniz rule.
        """
        cache = {} if cache is None else cache
        if item_id in cache:
            return cache[item_id]
        item = self.item(item_id)
        if item.kind == "placeholder":
            raise RuntimeError(
                f"{item.id} is a degree-{item.degree} import placeholder")
        if item.kind in {"graph", "graph_record"}:
            matrix = (
                graph_from_label(item.graph)
                if item.kind == "graph"
                else graph_from_record(item.graph_record)
            )
            scalar, gradient = invariant_value_and_derivative(
                matrix, five_form, 10, 5, True, mod)
        else:
            factors = [
                self.evaluate_item_with_gradient(
                    factor, five_form, mod, cache)
                for factor in item.factors
            ]
            scalar = 1
            for factor_value, _ in factors:
                scalar = scalar * factor_value % mod
            gradient = np.zeros((10,) * 5, dtype=np.int64)
            for index, (_, factor_gradient) in enumerate(factors):
                coefficient = 1
                for other, (factor_value, _) in enumerate(factors):
                    if other != index:
                        coefficient = coefficient * factor_value % mod
                gradient = (
                    gradient + coefficient * factor_gradient
                ) % mod
        cache[item_id] = (
            int(scalar) % mod,
            np.asarray(gradient, dtype=np.int64) % mod,
        )
        return cache[item_id]

    def with_degree12_primitives(self, primitives):
        """Return a copy with exactly 62 imported degree-12 primitives.

        Each primitive must be a concrete ``InvariantItem`` of degree 12.
        The record representation keeps the import independent of graph-label
        parsing conventions.
        """
        primitives = tuple(primitives)
        if len(primitives) != 62:
            raise ValueError("the degree-12 import requires 62 primitives")
        if any(item.degree != 12 or item.kind != "graph_record"
               for item in primitives):
            raise ValueError(
                "degree-12 primitives must use concrete degree-12 "
                "graph_record formulas")
        degrees = dict(self._degrees)
        degrees[12] = degree12_product_items() + primitives
        return InvariantRegistry(degrees, self.metadata)


def _load_json(path):
    with Path(path).open() as stream:
        return json.load(stream)


def load_verified_registry(repository_root):
    """Load the exact committed bases through degree 10."""
    root = Path(repository_root)
    lower_path = root / "results" / "10d_order8.json"
    degree10_path = root / "results" / "10d_order10.json"
    lower = _load_json(lower_path)
    upper = _load_json(degree10_path)

    lower_items = {
        item["id"]: InvariantItem(
            id=item["id"],
            degree=int(item["order"]),
            kind="graph",
            graph=item["graph"],
            source=str(lower_path.relative_to(root)),
        )
        for item in lower["generators"]
    }
    degree10_graphs = {
        item["id"]: InvariantItem(
            id=item["id"],
            degree=10,
            kind="graph",
            graph=item["graph"],
            source=str(degree10_path.relative_to(root)),
        )
        for item in upper["generators"]
    }
    degrees = {
        4: (lower_items["I4_1"],),
        6: (lower_items["I6_1"], lower_items["I6_2"]),
        8: tuple(lower_items[f"I8_{index}"] for index in range(1, 7))
        + (
            InvariantItem(
                id="I4_1^2",
                degree=8,
                kind="product",
                factors=("I4_1", "I4_1"),
                source=str(lower_path.relative_to(root)),
            ),
        ),
        10: tuple(
            degree10_graphs[f"I10_{index}"] for index in range(1, 13)
        ) + (
            InvariantItem(
                id="I4_1*I6_1",
                degree=10,
                kind="product",
                factors=("I4_1", "I6_1"),
                source=str(degree10_path.relative_to(root)),
            ),
            InvariantItem(
                id="I4_1*I6_2",
                degree=10,
                kind="product",
                factors=("I4_1", "I6_2"),
                source=str(degree10_path.relative_to(root)),
            ),
        ),
    }
    return InvariantRegistry(degrees, metadata={
        "schema": 1,
        "verified_through_degree": 10,
        "sources": [
            str(lower_path.relative_to(root)),
            str(degree10_path.relative_to(root)),
        ],
    })


def degree12_product_items():
    """The complete ten-dimensional product space at degree 12."""
    items = [
        InvariantItem(
            id="I4_1^3",
            degree=12,
            kind="product",
            factors=("I4_1", "I4_1", "I4_1"),
            source="degree12-interface",
        ),
        InvariantItem(
            id="I6_1^2",
            degree=12,
            kind="product",
            factors=("I6_1", "I6_1"),
            source="degree12-interface",
        ),
        InvariantItem(
            id="I6_1*I6_2",
            degree=12,
            kind="product",
            factors=("I6_1", "I6_2"),
            source="degree12-interface",
        ),
        InvariantItem(
            id="I6_2^2",
            degree=12,
            kind="product",
            factors=("I6_2", "I6_2"),
            source="degree12-interface",
        ),
    ]
    items.extend([
        InvariantItem(
            id=f"I4_1*I8_{index}",
            degree=12,
            kind="product",
            factors=("I4_1", f"I8_{index}"),
            source="degree12-interface",
        )
        for index in range(1, 7)
    ])
    return tuple(items)


def degree12_placeholder_items():
    """Return 62 labeled primitive slots for schema/adapter tests."""
    return tuple(
        InvariantItem(
            id=f"I12_primitive_{index:02d}",
            degree=12,
            kind="placeholder",
            source="awaiting-degree12-import",
        )
        for index in range(1, 63)
    )


def load_verified_registry_through_degree12(repository_root):
    """Load the committed 72-element degree-12 basis into the registry."""
    root = Path(repository_root)
    degree12_path = root / "results" / "10d_order12.json"
    raw = _load_json(degree12_path)
    generators = raw.get("generators", ())
    if len(generators) != 62:
        raise ValueError("degree-12 artifact must contain 62 generators")
    primitives = tuple(
        InvariantItem(
            id=item["id"],
            degree=12,
            kind="graph_record",
            graph_record=graph_to_record(graph_from_label(item["graph"])),
            source=str(degree12_path.relative_to(root)),
        )
        for item in generators
    )
    extended = load_verified_registry(root).with_degree12_primitives(
        primitives)
    artifact_basis = tuple(item["id"] for item in raw["degree12_basis"])
    loaded_basis = tuple(item.id for item in extended.basis(12))
    if loaded_basis != artifact_basis:
        raise ValueError(
            "registry degree-12 ordering differs from committed artifact")
    metadata = dict(extended.metadata)
    metadata.update({
        "verified_through_degree": 12,
        "degree12_source": str(degree12_path.relative_to(root)),
        "degree12_dimension": 72,
    })
    return InvariantRegistry(
        {degree: extended.basis(degree) for degree in extended.degrees},
        metadata,
    )
