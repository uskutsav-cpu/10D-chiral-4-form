"""Formal interacting stress-tensor traces through five-form degree 12.

The physical tensor is rescaled to ``tau = 48*T``.  This changes only the
normalization of arbitrary flow-coefficient functions and keeps all formal
coefficient maps integral:

``tau_mu^nu = M_mu^nu - 25 (D.D)_mu^nu
               + sum_d (d-2) V_d delta_mu^nu``,

where ``D=dV/dLambda^upper`` and ``V_d`` is homogeneous of field degree
``d``.  Sparse polynomials retain both the interaction-coefficient monomial
and the five-form degree, so truncation never confuses coupling order with
field order.
"""

from dataclasses import dataclass

import numpy as np

from .forms import metric_signs
from .graphs import graph_from_label, graph_from_record
from .interaction import invariant_value_derivative_hvp
from .modp import mod_einsum
from .stress import five_form_moment


MAX_FIELD_DEGREE = 12
INTERACTION_GRADIENT_DEGREES = (4, 6, 8)


@dataclass(frozen=True)
class TraceGenerator:
    """One independent scalar monomial in traces of ``tau=48*T``."""

    id: str
    trace_powers: tuple[int, ...]
    leading_field_degree: int


TRACE_GENERATORS = (
    TraceGenerator("tr_tau", (1,), 4),
    TraceGenerator("tr_tau2", (2,), 4),
    TraceGenerator("tr_tau3", (3,), 6),
    TraceGenerator("tr_tau4", (4,), 8),
    TraceGenerator("tr_tau^2", (1, 1), 8),
    TraceGenerator("tr_tau*tr_tau2", (1, 2), 8),
    TraceGenerator("tr_tau2^2", (2, 2), 8),
    TraceGenerator("tr_tau5", (5,), 10),
    TraceGenerator("tr_tau*tr_tau3", (1, 3), 10),
    TraceGenerator("tr_tau2*tr_tau3", (2, 3), 10),
    TraceGenerator("tr_tau6", (6,), 12),
    TraceGenerator("tr_tau*tr_tau4", (1, 4), 12),
    TraceGenerator("tr_tau2*tr_tau4", (2, 4), 12),
    TraceGenerator("tr_tau3^2", (3, 3), 12),
    TraceGenerator("tr_tau^3", (1, 1, 1), 12),
    TraceGenerator("tr_tau^2*tr_tau2", (1, 1, 2), 12),
    TraceGenerator("tr_tau*tr_tau2^2", (1, 2, 2), 12),
    TraceGenerator("tr_tau2^3", (2, 2, 2), 12),
)


def _canonical_monomial(monomial, variable_order):
    return tuple(sorted(
        monomial, key=lambda item: variable_order[item]))


def _add_matrix_term(polynomial, key, value, mod):
    value = np.asarray(value, dtype=np.int64) % mod
    if key in polynomial:
        polynomial[key] = (polynomial[key] + value) % mod
    else:
        polynomial[key] = value


def _matrix_polynomial_product(left, right, variable_order, mod,
                               max_degree=MAX_FIELD_DEGREE):
    result = {}
    for (left_degree, left_monomial), left_matrix in left.items():
        for (right_degree, right_monomial), right_matrix in right.items():
            degree = left_degree + right_degree
            if degree > max_degree:
                continue
            monomial = _canonical_monomial(
                left_monomial + right_monomial, variable_order)
            _add_matrix_term(
                result,
                (degree, monomial),
                left_matrix @ right_matrix % mod,
                mod,
            )
    return result


def _matrix_polynomial_product_directional(
    left,
    left_tangent,
    right,
    right_tangent,
    variable_order,
    mod,
    max_degree=MAX_FIELD_DEGREE,
):
    primal = _matrix_polynomial_product(
        left, right, variable_order, mod, max_degree)
    tangent = {}
    for first, second in (
        (left_tangent, right),
        (left, right_tangent),
    ):
        partial = _matrix_polynomial_product(
            first, second, variable_order, mod, max_degree)
        for key, matrix in partial.items():
            _add_matrix_term(tangent, key, matrix, mod)
    return primal, tangent


def _add_scalar_term(polynomial, key, value, mod):
    polynomial[key] = (
        polynomial.get(key, 0) + int(value)
    ) % mod


def _scalar_polynomial_product(left, right, variable_order, mod,
                               max_degree=MAX_FIELD_DEGREE):
    result = {}
    for (left_degree, left_monomial), left_value in left.items():
        for (right_degree, right_monomial), right_value in right.items():
            degree = left_degree + right_degree
            if degree > max_degree:
                continue
            monomial = _canonical_monomial(
                left_monomial + right_monomial, variable_order)
            _add_scalar_term(
                result,
                (degree, monomial),
                left_value * right_value,
                mod,
            )
    return result


def _scalar_polynomial_product_directional(
    left,
    left_tangent,
    right,
    right_tangent,
    variable_order,
    mod,
    max_degree=MAX_FIELD_DEGREE,
):
    primal = _scalar_polynomial_product(
        left, right, variable_order, mod, max_degree)
    tangent = {}
    for first, second in (
        (left_tangent, right),
        (left, right_tangent),
    ):
        partial = _scalar_polynomial_product(
            first, second, variable_order, mod, max_degree)
        for key, value in partial.items():
            _add_scalar_term(tangent, key, value, mod)
    return primal, tangent


def _bilinear_derivative_moment_mixed(left, right, mod):
    """Return ``left_mu,rho(4) right^nu,rho(4)`` as a mixed tensor."""
    left = np.asarray(left, dtype=np.int64) % mod
    right = np.asarray(right, dtype=np.int64) % mod
    signs = metric_signs(10, True) % mod
    raised_tail = right
    for axis in (1, 2, 3, 4):
        shape = [1] * 5
        shape[axis] = 10
        raised_tail = raised_tail * signs.reshape(shape) % mod
    lower = mod_einsum(
        "mabcd,nabcd->mn", [left, raised_tail], mod)
    return np.asarray(
        lower * signs.reshape(1, 10) % mod,
        dtype=np.int64,
    )


def interaction_data(five_form, registry, mod, value_degrees=None,
                     gradient_degrees=INTERACTION_GRADIENT_DEGREES):
    """Evaluate the basis data needed by traces through degree 12."""
    value_degrees = (
        tuple(registry.degrees)
        if value_degrees is None
        else tuple(int(degree) for degree in value_degrees)
    )
    gradient_degrees = tuple(int(degree) for degree in gradient_degrees)
    unknown = (set(value_degrees) | set(gradient_degrees)) - set(
        registry.degrees)
    if unknown:
        raise ValueError(f"unregistered interaction degrees: {sorted(unknown)}")

    derivative_cache = {}
    values = {}
    derivatives = {}
    for degree in gradient_degrees:
        for item in registry.basis(degree):
            value, derivative = registry.evaluate_item_with_gradient(
                item.id, five_form, mod, derivative_cache)
            values[item.id] = value
            derivatives[item.id] = derivative

    value_cache = dict(values)
    for degree in value_degrees:
        for item in registry.basis(degree):
            if item.id in values:
                continue
            values[item.id] = registry.evaluate_item(
                item.id, five_form, mod, value_cache)
    return values, derivatives


def interaction_directional_data(
    five_form,
    direction,
    registry,
    mod,
    degrees=INTERACTION_GRADIENT_DEGREES,
):
    """Evaluate interaction values/gradients and their directional changes."""
    degrees = tuple(int(degree) for degree in degrees)
    cache = {}

    def evaluate(item_id):
        if item_id in cache:
            return cache[item_id]
        item = registry.item(item_id)
        if item.kind in {"graph", "graph_record"}:
            matrix = (
                graph_from_label(item.graph)
                if item.kind == "graph"
                else graph_from_record(item.graph_record)
            )
            result = invariant_value_derivative_hvp(
                matrix, five_form, direction, 10, 5, True, mod)
        elif item.kind == "product":
            # Dual product rule for both the scalar and its form gradient.
            result = (
                1,
                0,
                np.zeros((10,) * 5, dtype=np.int64),
                np.zeros((10,) * 5, dtype=np.int64),
            )
            for factor in item.factors:
                left_value, left_tangent, left_gradient, left_hvp = result
                value, tangent, gradient, hvp = evaluate(factor)
                result = (
                    left_value * value % mod,
                    (
                        left_tangent * value
                        + left_value * tangent
                    ) % mod,
                    (
                        left_gradient * value
                        + left_value * gradient
                    ) % mod,
                    (
                        left_hvp * value
                        + left_gradient * tangent
                        + left_tangent * gradient
                        + left_value * hvp
                    ) % mod,
                )
        else:
            raise RuntimeError(
                f"{item.id} has no concrete interaction formula")
        cache[item_id] = (
            int(result[0]) % mod,
            int(result[1]) % mod,
            np.asarray(result[2], dtype=np.int64) % mod,
            np.asarray(result[3], dtype=np.int64) % mod,
        )
        return cache[item_id]

    for degree in degrees:
        for item in registry.basis(degree):
            evaluate(item.id)
    return {
        item_id: record[0] for item_id, record in cache.items()
    }, {
        item_id: record[1] for item_id, record in cache.items()
    }, {
        item_id: record[2] for item_id, record in cache.items()
    }, {
        item_id: record[3] for item_id, record in cache.items()
    }


def normalized_stress_polynomial(five_form, registry, mod,
                                 values=None, derivatives=None,
                                 max_degree=MAX_FIELD_DEGREE):
    """Return the sparse mixed-tensor polynomial for ``tau=48*T``.

    Terms with a derivative-moment degree above ten cannot enter a nonzero
    trace through total field degree twelve and are excluded.  Scalar terms
    through degree eight suffice for traces of power two and higher; the
    exact trace-one polynomial is constructed separately.
    """
    if values is None or derivatives is None:
        values, derivatives = interaction_data(five_form, registry, mod)
    variable_ids = [
        item.id
        for degree in registry.degrees
        for item in registry.basis(degree)
    ]
    variable_order = {
        item_id: index for index, item_id in enumerate(variable_ids)}
    _, m_mixed = five_form_moment(five_form, mod)
    polynomial = {(2, ()): m_mixed}
    identity = np.eye(10, dtype=np.int64)

    active_items = [
        item
        for degree in INTERACTION_GRADIENT_DEGREES
        for item in registry.basis(degree)
    ]
    for item in active_items:
        # A scalar V_8 term can first enter Tr(tau^2) at degree 12.
        _add_matrix_term(
            polynomial,
            (item.degree, (item.id,)),
            (item.degree - 2) * values[item.id] * identity,
            mod,
        )

    for left_index, left_item in enumerate(active_items):
        for right_item in active_items[left_index:]:
            derivative_degree = left_item.degree + right_item.degree - 2
            # A traceless Q_12 cannot contribute below field degree 14.
            if derivative_degree > max_degree - 2:
                continue
            left = derivatives[left_item.id]
            right = derivatives[right_item.id]
            moment = _bilinear_derivative_moment_mixed(left, right, mod)
            if left_item.id != right_item.id:
                moment = (
                    moment
                    + _bilinear_derivative_moment_mixed(
                        right, left, mod)
                ) % mod
            monomial = _canonical_monomial(
                (left_item.id, right_item.id), variable_order)
            _add_matrix_term(
                polynomial,
                (derivative_degree, monomial),
                -25 * moment,
                mod,
            )
    return polynomial, variable_order


def normalized_stress_polynomial_directional(
    five_form,
    direction,
    registry,
    mod,
    values,
    value_tangents,
    derivatives,
    derivative_tangents,
    max_degree=MAX_FIELD_DEGREE,
):
    """Return ``tau`` and its field-directional derivative polynomials."""
    variable_ids = [
        item.id
        for degree in registry.degrees
        for item in registry.basis(degree)
    ]
    variable_order = {
        item_id: index for index, item_id in enumerate(variable_ids)}
    _, m_mixed = five_form_moment(five_form, mod)
    m_tangent = (
        _bilinear_derivative_moment_mixed(
            five_form, direction, mod)
        + _bilinear_derivative_moment_mixed(
            direction, five_form, mod)
    ) % mod
    polynomial = {(2, ()): m_mixed}
    tangent_polynomial = {(2, ()): m_tangent}
    identity = np.eye(10, dtype=np.int64)
    active_items = [
        item
        for degree in INTERACTION_GRADIENT_DEGREES
        for item in registry.basis(degree)
    ]
    for item in active_items:
        key = (item.degree, (item.id,))
        _add_matrix_term(
            polynomial,
            key,
            (item.degree - 2) * values[item.id] * identity,
            mod,
        )
        _add_matrix_term(
            tangent_polynomial,
            key,
            (item.degree - 2)
            * value_tangents[item.id]
            * identity,
            mod,
        )

    for left_index, left_item in enumerate(active_items):
        for right_item in active_items[left_index:]:
            derivative_degree = left_item.degree + right_item.degree - 2
            if derivative_degree > max_degree - 2:
                continue
            left = derivatives[left_item.id]
            right = derivatives[right_item.id]
            left_tangent = derivative_tangents[left_item.id]
            right_tangent = derivative_tangents[right_item.id]
            moment = _bilinear_derivative_moment_mixed(
                left, right, mod)
            moment_tangent = (
                _bilinear_derivative_moment_mixed(
                    left_tangent, right, mod)
                + _bilinear_derivative_moment_mixed(
                    left, right_tangent, mod)
            ) % mod
            if left_item.id != right_item.id:
                moment = (
                    moment
                    + _bilinear_derivative_moment_mixed(
                        right, left, mod)
                ) % mod
                moment_tangent = (
                    moment_tangent
                    + _bilinear_derivative_moment_mixed(
                        right_tangent, left, mod)
                    + _bilinear_derivative_moment_mixed(
                        right, left_tangent, mod)
                ) % mod
            monomial = _canonical_monomial(
                (left_item.id, right_item.id), variable_order)
            key = (derivative_degree, monomial)
            _add_matrix_term(
                polynomial, key, -25 * moment, mod)
            _add_matrix_term(
                tangent_polynomial, key, -25 * moment_tangent, mod)
    return polynomial, tangent_polynomial, variable_order


def trace_polynomials(five_form, registry, mod, values=None,
                      derivatives=None, max_degree=MAX_FIELD_DEGREE):
    """Return sparse polynomials for ``Tr(tau^k)``, ``k=1,...,6``."""
    if values is None or derivatives is None:
        values, derivatives = interaction_data(five_form, registry, mod)
    tensor, variable_order = normalized_stress_polynomial(
        five_form,
        registry,
        mod,
        values,
        derivatives,
        max_degree,
    )

    # The derivative-moment term is traceless for anti-self-dual D, while
    # Tr(M)=0.  This exact formula also includes V_10 and V_12 without
    # computing their expensive derivatives.
    trace_one = {}
    for degree in registry.degrees:
        if degree > max_degree:
            continue
        for item in registry.basis(degree):
            if item.id not in values:
                continue
            _add_scalar_term(
                trace_one,
                (degree, (item.id,)),
                10 * (degree - 2) * values[item.id],
                mod,
            )

    identity = {(0, ()): np.eye(10, dtype=np.int64)}
    power = identity
    traces = {1: trace_one}
    for exponent in range(1, 7):
        power = _matrix_polynomial_product(
            power, tensor, variable_order, mod, max_degree)
        if exponent >= 2:
            traces[exponent] = {
                key: int(np.trace(matrix) % mod)
                for key, matrix in power.items()
            }
    return traces, variable_order


def trace_polynomials_directional(
    five_form,
    direction,
    registry,
    mod,
    values,
    value_tangents,
    derivatives,
    derivative_tangents,
    max_degree=MAX_FIELD_DEGREE,
):
    """Return trace polynomials and exact field-directional derivatives."""
    tensor, tensor_tangent, variable_order = (
        normalized_stress_polynomial_directional(
            five_form,
            direction,
            registry,
            mod,
            values,
            value_tangents,
            derivatives,
            derivative_tangents,
            max_degree,
        )
    )
    trace_one = {}
    trace_one_tangent = {}
    for degree in registry.degrees:
        if degree > max_degree:
            continue
        for item in registry.basis(degree):
            if item.id not in values:
                continue
            key = (degree, (item.id,))
            _add_scalar_term(
                trace_one,
                key,
                10 * (degree - 2) * values[item.id],
                mod,
            )
            _add_scalar_term(
                trace_one_tangent,
                key,
                10 * (degree - 2) * value_tangents[item.id],
                mod,
            )

    power = {(0, ()): np.eye(10, dtype=np.int64)}
    power_tangent = {
        (0, ()): np.zeros((10, 10), dtype=np.int64)}
    traces = {1: trace_one}
    trace_tangents = {1: trace_one_tangent}
    for exponent in range(1, 7):
        power, power_tangent = (
            _matrix_polynomial_product_directional(
                power,
                power_tangent,
                tensor,
                tensor_tangent,
                variable_order,
                mod,
                max_degree,
            )
        )
        if exponent >= 2:
            traces[exponent] = {
                key: int(np.trace(matrix) % mod)
                for key, matrix in power.items()
            }
            trace_tangents[exponent] = {
                key: int(np.trace(matrix) % mod)
                for key, matrix in power_tangent.items()
            }
    return traces, trace_tangents, variable_order


def flow_generator_polynomials(five_form, registry, mod, values=None,
                               derivatives=None,
                               max_degree=MAX_FIELD_DEGREE):
    """Expand every independent analytic scalar generator through degree 12."""
    traces, variable_order = trace_polynomials(
        five_form,
        registry,
        mod,
        values,
        derivatives,
        max_degree,
    )
    result = {}
    for generator in TRACE_GENERATORS:
        polynomial = {(0, ()): 1}
        for power in generator.trace_powers:
            polynomial = _scalar_polynomial_product(
                polynomial,
                traces[power],
                variable_order,
                mod,
                max_degree,
            )
        result[generator.id] = polynomial
    return result


def flow_generator_directional_polynomials(
    five_form,
    direction,
    registry,
    mod,
    max_degree=MAX_FIELD_DEGREE,
    return_interaction_data=False,
):
    """Expand all generators and their exact directional derivatives."""
    data = interaction_directional_data(
        five_form, direction, registry, mod)
    traces, trace_tangents, variable_order = (
        trace_polynomials_directional(
            five_form,
            direction,
            registry,
            mod,
            *data,
            max_degree,
        )
    )
    result = {}
    tangents = {}
    for generator in TRACE_GENERATORS:
        polynomial = {(0, ()): 1}
        tangent = {(0, ()): 0}
        for power in generator.trace_powers:
            polynomial, tangent = (
                _scalar_polynomial_product_directional(
                    polynomial,
                    tangent,
                    traces[power],
                    trace_tangents[power],
                    variable_order,
                    mod,
                    max_degree,
                )
            )
        result[generator.id] = polynomial
        tangents[generator.id] = tangent
    if return_interaction_data:
        return result, tangents, data
    return result, tangents


def evaluate_matrix_polynomial(polynomial, coefficients, mod):
    """Evaluate a sparse matrix polynomial at interaction coefficients."""
    result = np.zeros((10, 10), dtype=np.int64)
    for (_, monomial), matrix in polynomial.items():
        scalar = 1
        for item_id in monomial:
            scalar = scalar * int(coefficients.get(item_id, 0)) % mod
        result = (result + scalar * matrix) % mod
    return result
