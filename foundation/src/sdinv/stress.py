"""Paper-normalized stress tensors for a self-dual five-form in D=10.

The conventions follow Hutomo--Lechner--Sorokin, arXiv:2509.14351v2:

* metric ``eta = diag(-1,+1,...,+1)``;
* ``Lambda = *Lambda`` with ``epsilon^{01...9} = -1``;
* ``M_mu^nu = Lambda_{mu rho(4)} Lambda^{nu rho(4)}``;
* the free stress tensor is ``T_{mu nu}=M_{mu nu}/(2*4!)``.

All public numerical routines operate exactly over an odd prime field.  The
``reference`` backends deliberately use direct two-operand ``numpy.einsum``
calls.  Their reduction bounds are small enough to avoid int64 overflow and
they provide a simple oracle for the optimized pairwise implementation.

The full 1050 and 4125 components are available as convention tests.  The
production stress path materializes only the contraction of the 4125 irrep
which enters the stress tensor:

    (N^(4125) M M)_{mu nu}.

Equation (2.17), inserted into equations (3.3)--(3.4), gives it exactly in
terms of the unprojected N tensor.  The identity is tested independently.
"""

import itertools
from math import factorial

import numpy as np

from .forms import metric_signs, perm_sign
from .modp import P, inv, mod_einsum


D = 10
PDEG = 5
FOUR_FACTORIAL = factorial(4)


def _require_five_form(five_form):
    tensor = np.asarray(five_form, dtype=np.int64)
    if tensor.shape != (D,) * PDEG:
        raise ValueError(
            f"five-form must have shape {(D,) * PDEG}, got {tensor.shape}")
    return tensor


def _metric(mod):
    signs = metric_signs(D, True) % mod
    return signs, np.diag(signs).astype(np.int64) % mod


def _raise_axes(tensor, axes, mod):
    """Raise selected covariant axes with the diagonal Lorentz metric."""
    result = np.asarray(tensor, dtype=np.int64) % mod
    signs = metric_signs(D, True) % mod
    for axis in axes:
        shape = [1] * result.ndim
        shape[int(axis)] = D
        result = (result * signs.reshape(shape)) % mod
    return result


def _scale(tensor, factors, mod):
    """Multiply an array by modular scalars with reduction at every step."""
    result = np.asarray(tensor, dtype=np.int64) % mod
    for factor in factors:
        result = result * (int(factor) % mod) % mod
    return result


def matrix_trace_power(mixed_matrix, power, mod=P):
    """Return ``Tr(A**power)`` over ``F_mod`` for a mixed tensor A_mu^nu."""
    if power < 1:
        raise ValueError("power must be positive")
    matrix = np.asarray(mixed_matrix, dtype=np.int64) % mod
    if matrix.shape != (D, D):
        raise ValueError(f"matrix must have shape {(D, D)}")
    product = np.eye(D, dtype=np.int64)
    for _ in range(int(power)):
        product = (product @ matrix) % mod
    return int(np.trace(product) % mod)


def symmetric_inner(left_lower, right_lower, mod=P):
    """Lorentz contraction ``A_{mu nu} B^{mu nu}`` over ``F_mod``."""
    left = np.asarray(left_lower, dtype=np.int64) % mod
    right = np.asarray(right_lower, dtype=np.int64) % mod
    if left.shape != (D, D) or right.shape != (D, D):
        raise ValueError("both tensors must have shape (10, 10)")
    signs = metric_signs(D, True) % mod
    right_upper = (
        right * signs.reshape(D, 1) * signs.reshape(1, D)
    ) % mod
    return int(np.sum(left * right_upper, dtype=np.int64) % mod)


def five_form_moment(five_form, mod=P, backend="optimized"):
    """Return ``(M_lower, M_mixed)`` for the paper's quadratic moment.

    ``M_lower[mu,nu] = Lambda_{mu rho(4)} Lambda_nu^{ rho(4)}`` and
    ``M_mixed = M_lower * eta^{-1}`` on its second index.
    """
    five_form = _require_five_form(five_form) % mod
    raised_tail = _raise_axes(five_form, (1, 2, 3, 4), mod)
    if backend == "optimized":
        lower = mod_einsum(
            "mabcd,nabcd->mn", [five_form, raised_tail], mod)
    elif backend == "reference":
        # At most 10^4 terms of size < mod^2 enter each output, safely below
        # int64.  Keeping this direct expression makes it a useful oracle.
        lower = np.einsum(
            "mabcd,nabcd->mn", five_form, raised_tail) % mod
    else:
        raise ValueError(f"unknown backend: {backend!r}")
    signs = metric_signs(D, True) % mod
    mixed = (lower * signs.reshape(1, D)) % mod
    return lower.astype(np.int64), mixed.astype(np.int64)


def composite_n(five_form, mod=P, backend="optimized"):
    """Return ``N_{abc,def}=Lambda_{abc rs} Lambda_{def}{}^{rs}``.

    This is the unprojected composite tensor of paper equation (2.16).
    """
    five_form = _require_five_form(five_form) % mod
    raised_pair = _raise_axes(five_form, (3, 4), mod)
    if backend == "optimized":
        result = mod_einsum(
            "abcrs,defrs->abcdef", [five_form, raised_pair], mod)
    elif backend == "reference":
        # Each output sums only 10^2 products, so int64 is exact here.
        result = np.einsum(
            "abcrs,defrs->abcdef", five_form, raised_pair) % mod
    else:
        raise ValueError(f"unknown backend: {backend!r}")
    return np.asarray(result, dtype=np.int64)


def _antisymmetrize_axes(tensor, axes, mod):
    """Normalized antisymmetrization over selected tensor axes."""
    axes = tuple(int(axis) for axis in axes)
    result = np.zeros_like(tensor, dtype=np.int64)
    for permutation in itertools.permutations(range(len(axes))):
        order = list(range(tensor.ndim))
        for destination, source in enumerate(permutation):
            order[axes[destination]] = axes[source]
        sign = perm_sign(permutation)
        result = (result + sign * tensor.transpose(order)) % mod
    return result * inv(factorial(len(axes)), mod) % mod


def interaction_gradient_i4(five_form, v_i, mod=P, backend="optimized"):
    """Return ``dV/dLambda^{mu(5)}`` for ``V=V(I4)``.

    Equation (3.2) of arXiv:2509.14351v2 states

    ``F^- = 20 M_[mu1]^rho Lambda_mu2...mu5 rho V_I``,

    while equation (2.20) states ``F^-=5 dV/dLambda^upper``.  Consequently
    the covariant, anti-self-dual derivative needed in equation (2.33) is

    ``dV/dLambda^{mu(5)}
      = 4 M_[mu1]^rho Lambda_mu2...mu5 rho V_I``.

    Square brackets are normalized throughout this module.
    """
    five_form = _require_five_form(five_form) % mod
    _, m_mixed = five_form_moment(five_form, mod, backend)
    if backend == "optimized":
        raw = mod_einsum(
            "mr,abcdr->mabcd", [m_mixed, five_form], mod)
    elif backend == "reference":
        # Only ten products enter each output component.
        raw = np.einsum(
            "mr,abcdr->mabcd", m_mixed, five_form) % mod
    else:
        raise ValueError(f"unknown backend: {backend!r}")
    derivative = _antisymmetrize_axes(raw, range(PDEG), mod)
    return _scale(derivative, (4, int(v_i)), mod)


def interacting_stress(five_form, derivative_lower, v, euler, mod=P,
                       backend="optimized"):
    """Return the fully interacting INZ stress tensor, equation (2.33).

    Parameters are exact field values:

    * ``derivative_lower`` is the covariant anti-self-dual tensor
      ``dV/dLambda^{mu(5)}``;
    * ``v`` is ``V(Lambda)``;
    * ``euler`` is ``Lambda_(5) dV/dLambda_(5)``.

    Passing the Euler contraction explicitly avoids silently choosing an
    ambient extension of a function defined only on the self-dual subspace.
    For a homogeneous field-degree-``d`` interaction it is simply ``d*v``.
    """
    five_form = _require_five_form(five_form) % mod
    derivative = _require_five_form(derivative_lower) % mod
    m_lower, _ = five_form_moment(five_form, mod, backend)
    derivative_moment, _ = five_form_moment(derivative, mod, backend)
    _, eta = _metric(mod)
    tensor_part = _scale(
        (m_lower - 25 * derivative_moment) % mod,
        (inv(2 * FOUR_FACTORIAL, mod),),
        mod,
    )
    scalar_part = _scale(
        eta,
        (
            -inv(FOUR_FACTORIAL, mod),
            (int(v) - inv(2, mod) * int(euler)) % mod,
        ),
        mod,
    )
    return np.asarray((tensor_part + scalar_part) % mod, dtype=np.int64)


def interacting_trace_formula(v, euler, mod=P):
    """Trace implied algebraically by the displayed equation (2.33).

    With ``eta=(-,+,...,+)``, its first bracket is traceless and the
    ``-g_mu_nu/4!`` term gives

    ``Tr(T) = -5/12 * (V - (Lambda.dV)/2)``.

    Equation (2.36) of arXiv:2509.14351v2 displays the opposite overall sign;
    the repository records this internal source discrepancy explicitly
    instead of silently changing equation (2.33).
    """
    return int(
        -5 * inv(12, mod)
        * (int(v) - inv(2, mod) * int(euler))
        % mod
    )


def composite_n1050(five_form, mod=P, backend="optimized"):
    """Return the 1050 component defined in paper equation (2.15).

    The first five indices are normalized-antisymmetrized:

    ``N1050_[abc,de]f = Lambda^{mn}_[abc Lambda_{de]fmn}``.
    """
    five_form = _require_five_form(five_form) % mod
    raised_pair = _raise_axes(five_form, (0, 1), mod)
    if backend == "optimized":
        raw = mod_einsum(
            "mnabc,defmn->abcdef", [raised_pair, five_form], mod)
    elif backend == "reference":
        raw = np.einsum(
            "mnabc,defmn->abcdef", raised_pair, five_form) % mod
    else:
        raise ValueError(f"unknown backend: {backend!r}")
    return _antisymmetrize_axes(raw, (0, 1, 2, 3, 4), mod)


def _n_trace_54(m_mixed, mod):
    """The 54 trace term in the N decomposition, equation (2.17)."""
    base = np.zeros((D,) * 6, dtype=np.int64)
    for first in range(D):
        for second in range(D):
            base[first, second, :, first, second, :] = m_mixed
    projected = _antisymmetrize_axes(base, (0, 1, 2), mod)
    projected = _antisymmetrize_axes(projected, (3, 4, 5), mod)
    # The alpha indices in equation (2.17) are contracted with three lower
    # metrics to become the second covariant 3-block.
    projected = _raise_axes(projected, (3, 4, 5), mod)
    return _scale(projected, (9, inv(28, mod)), mod)


def composite_n4125(five_form, mod=P, backend="optimized"):
    """Materialize the full 4125 irrep in paper equation (2.17).

    This reference-facing routine is intentionally separate from the faster
    stress path, which needs only ``(N^(4125)MM)_{mu nu}``.
    """
    n_lower = composite_n(five_form, mod, backend)
    n1050 = composite_n1050(five_form, mod, backend)
    # The red antisymmetrization in equation (2.17) is performed after the
    # five-index (black-bracket) antisymmetrization in N1050.
    n1050_term = _scale(
        _antisymmetrize_axes(n1050, (3, 4, 5), mod),
        (5,),
        mod,
    )
    _, m_mixed = five_form_moment(five_form, mod, backend)
    trace_term = _n_trace_54(m_mixed, mod)
    return (n_lower - n1050_term - trace_term) % mod


def raw_n_mm(n_lower, m_mixed, mod=P):
    """Contract the unprojected N tensor with two copies of M.

    Returns
    ``N_{mu a b,nu}{}^{c d} M_c{}^a M_d{}^b`` with both free indices
    lowered, as it appears in paper equation (3.3).
    """
    tensor = np.asarray(n_lower, dtype=np.int64) % mod
    matrix = np.asarray(m_mixed, dtype=np.int64) % mod
    if tensor.shape != (D,) * 6:
        raise ValueError("N must have shape (10,)*6")
    if matrix.shape != (D, D):
        raise ValueError("M must have shape (10,10)")
    # Raise the last two indices of N; the first index of its second 3-block
    # (nu) remains lowered.
    tensor = _raise_axes(tensor, (4, 5), mod)
    return np.asarray(mod_einsum(
        "mabncd,ca,db->mn", [tensor, matrix, matrix], mod),
        dtype=np.int64,
    )


def _m_power_lower(m_mixed, power, mod):
    mixed = np.eye(D, dtype=np.int64)
    for _ in range(int(power)):
        mixed = (mixed @ m_mixed) % mod
    signs = metric_signs(D, True) % mod
    return (mixed * signs.reshape(1, D)) % mod


def stress_correction_r(five_form, mod=P, backend="optimized",
                        return_intermediates=False):
    """Return the paper's symmetric traceless ModMax correction R.

    We define

    ``R = 5/7 (M^3 - eta Tr(M^3)/10) - 12 (N^(4125) M M)``.

    Rather than implementing a convention-sensitive six-index Young
    projector, equation (2.17) gives the exactly equivalent expression

    ``R = 2 M^3 - eta Tr(M^3)/2 - 12 (N M M) - 9 I4 M/14``.
    """
    m_lower, m_mixed = five_form_moment(five_form, mod, backend)
    n_lower = composite_n(five_form, mod, backend)
    nmm = raw_n_mm(n_lower, m_mixed, mod)
    m3_lower = _m_power_lower(m_mixed, 3, mod)
    tr_m3 = matrix_trace_power(m_mixed, 3, mod)
    i4 = matrix_trace_power(m_mixed, 2, mod)
    _, eta = _metric(mod)
    half = inv(2, mod)
    nine_fourteenths = 9 * inv(14, mod) % mod
    result = (
        2 * m3_lower
        - half * tr_m3 * eta
        - 12 * nmm
        - nine_fourteenths * i4 * m_lower
    ) % mod
    if return_intermediates:
        return result, {
            "M_lower": m_lower,
            "M_mixed": m_mixed,
            "N_lower": n_lower,
            "NMM_lower": nmm,
            "M3_lower": m3_lower,
            "tr_M3": tr_m3,
            "I4": i4,
        }
    return result


def n4125_mm(five_form, mod=P, backend="optimized"):
    """Return the relevant ``(N^(4125) M M)_{mu nu}`` contraction.

    This follows by solving the defining equation for R.  It is the complete
    4125-irrep information needed by equations (3.4), (3.12)--(3.15).
    """
    r_tensor, data = stress_correction_r(
        five_form, mod, backend, return_intermediates=True)
    m3_traceless = (
        data["M3_lower"]
        - inv(10, mod) * data["tr_M3"] * _metric(mod)[1]
    ) % mod
    return (
        (5 * inv(7, mod) % mod) * m3_traceless - r_tensor
    ) * inv(12, mod) % mod


def paper_i4_i8_i12(five_form, mod=P, backend="optimized"):
    """Return the paper's ``(I4, Tr(M^3), I8, I12)`` exactly.

    Equations (3.14)--(3.15) are equivalently
    ``I8=M^{mu nu}R_{mu nu}`` and ``I12=R^{mu nu}R_{mu nu}``.
    """
    r_tensor, data = stress_correction_r(
        five_form, mod, backend, return_intermediates=True)
    i8 = symmetric_inner(data["M_lower"], r_tensor, mod)
    i12 = symmetric_inner(r_tensor, r_tensor, mod)
    return {
        "I4": int(data["I4"]),
        "tr_M3": int(data["tr_M3"]),
        "I8": int(i8),
        "I12": int(i12),
    }


def paper_i8_i12_expanded(five_form, mod=P, backend="optimized"):
    """Evaluate equations (3.14)--(3.15) in their displayed expanded form.

    This is intentionally independent of the compact ``M.R`` / ``R.R``
    evaluation in :func:`paper_i4_i8_i12` and is used as a regression oracle.
    """
    r_tensor, data = stress_correction_r(
        five_form, mod, backend, return_intermediates=True)
    projected_nmm = n4125_mm(five_form, mod, backend)
    m_mixed = data["M_mixed"]
    m_lower = data["M_lower"]
    m3_lower = data["M3_lower"]
    tr_m3 = data["tr_M3"]
    tr_m4 = matrix_trace_power(m_mixed, 4, mod)
    tr_m6 = matrix_trace_power(m_mixed, 6, mod)
    nmm_m = symmetric_inner(projected_nmm, m_lower, mod)
    nmm_m3 = symmetric_inner(projected_nmm, m3_lower, mod)
    nmm_squared = symmetric_inner(projected_nmm, projected_nmm, mod)
    i8 = (
        5 * inv(7, mod) * tr_m4 - 12 * nmm_m
    ) % mod
    i12 = (
        25 * inv(49, mod) * tr_m6
        - 5 * inv(98, mod) * tr_m3 * tr_m3
        - 120 * inv(7, mod) * nmm_m3
        + 144 * nmm_squared
    ) % mod
    return {
        "I8": int(i8),
        "I12": int(i12),
        "N4125MM_lower": projected_nmm,
        "R_lower": r_tensor,
    }


def stress_v_i4(five_form, v, v_i, mod=P, backend="optimized"):
    """Stress tensor for an interaction ``V(I4)`` (paper equation 3.3).

    ``v`` and ``v_i=dV/dI4`` are field values in ``F_mod``.  Supplying them
    separately also supports algebraic extensions such as ``sqrt(I4)``
    without choosing a finite-field square-root branch.
    """
    r_tensor, data = stress_correction_r(
        five_form, mod, backend, return_intermediates=True)
    m_lower = data["M_lower"]
    i4 = data["I4"]
    _, eta = _metric(mod)
    v = int(v) % mod
    v_i = int(v_i) % mod
    v_i_squared = v_i * v_i % mod
    first_coefficient = (
        1 - 96 * inv(7, mod) * v_i_squared * i4
    ) % mod
    first = _scale(
        m_lower, (inv(2 * FOUR_FACTORIAL, mod), first_coefficient), mod)
    scalar_coefficient = (
        -inv(FOUR_FACTORIAL, mod) * (v - 2 * i4 * v_i)
    ) % mod
    scalar = _scale(eta, (scalar_coefficient,), mod)
    correction = _scale(
        r_tensor, (inv(3, mod), v_i_squared), mod)
    return np.asarray((first + scalar + correction) % mod, dtype=np.int64)


def stress_v_i4_raw(five_form, v, v_i, mod=P, backend="optimized"):
    """Reference equation (3.3), before the N-irrep decomposition."""
    m_lower, m_mixed = five_form_moment(five_form, mod, backend)
    n_lower = composite_n(five_form, mod, backend)
    nmm = raw_n_mm(n_lower, m_mixed, mod)
    m3_lower = _m_power_lower(m_mixed, 3, mod)
    tr_m3 = matrix_trace_power(m_mixed, 3, mod)
    i4 = matrix_trace_power(m_mixed, 2, mod)
    _, eta = _metric(mod)
    v = int(v) % mod
    v_i = int(v_i) % mod
    v_i_squared = v_i * v_i % mod
    first_coefficient = (1 - 24 * v_i_squared * i4) % mod
    first = _scale(
        m_lower, (inv(2 * FOUR_FACTORIAL, mod), first_coefficient), mod)
    scalar_coefficient = (
        -inv(FOUR_FACTORIAL, mod) * (v - 2 * i4 * v_i)
    ) % mod
    scalar = _scale(eta, (scalar_coefficient,), mod)
    bracket = (
        2 * m3_lower
        - inv(2, mod) * tr_m3 * eta
        - 12 * nmm
    ) % mod
    correction = _scale(bracket, (inv(3, mod), v_i_squared), mod)
    return np.asarray((first + scalar + correction) % mod, dtype=np.int64)


def modmax_stress(five_form, b, mod=P, backend="optimized"):
    """D=10 ModMax stress tensor for ``V=b*sqrt(I4)``.

    The square root cancels from T.  This is paper equation (3.12):

    ``T = (1-24 b^2/7) M/(2*4!) + b^2 R/(12 I4)``.
    """
    r_tensor, data = stress_correction_r(
        five_form, mod, backend, return_intermediates=True)
    i4 = int(data["I4"]) % mod
    if i4 == 0:
        raise ZeroDivisionError("ModMax stress is singular on I4=0")
    b_squared = int(b) * int(b) % mod
    first_coefficient = (
        1 - 24 * inv(7, mod) * b_squared
    ) % mod
    first = _scale(
        data["M_lower"],
        (inv(2 * FOUR_FACTORIAL, mod), first_coefficient),
        mod,
    )
    correction = _scale(
        r_tensor, (b_squared, inv(12, mod), inv(i4, mod)), mod)
    return np.asarray((first + correction) % mod, dtype=np.int64)


def free_stress(five_form, mod=P, backend="optimized"):
    """Free chiral four-form stress tensor, ``M/(2*4!)``."""
    m_lower, _ = five_form_moment(five_form, mod, backend)
    return _scale(m_lower, (inv(2 * FOUR_FACTORIAL, mod),), mod)


def stress_mixed(stress_lower, mod=P):
    """Raise the second index of a covariant symmetric stress tensor."""
    tensor = np.asarray(stress_lower, dtype=np.int64) % mod
    if tensor.shape != (D, D):
        raise ValueError("stress tensor must have shape (10,10)")
    signs = metric_signs(D, True) % mod
    return tensor * signs.reshape(1, D) % mod


def stress_traces(stress_lower, max_power=10, mod=P):
    """Return ``Tr(T^k)`` for ``k=1,...,max_power``."""
    if not 1 <= max_power <= D:
        raise ValueError("require 1 <= max_power <= the D=10 CH limit")
    mixed = stress_mixed(stress_lower, mod)
    return {
        power: matrix_trace_power(mixed, power, mod)
        for power in range(1, max_power + 1)
    }


def modmax_stress_square_formula(five_form, b, mod=P,
                                 backend="optimized"):
    """Evaluate the exact right side of paper equations (3.13)--(3.15)."""
    invariants = paper_i4_i8_i12(five_form, mod, backend)
    i4, i8, i12 = (
        invariants["I4"], invariants["I8"], invariants["I12"])
    if i4 == 0:
        raise ZeroDivisionError("ModMax formula is singular on I4=0")
    b_squared = int(b) * int(b) % mod
    a = (1 - 24 * inv(7, mod) * b_squared) % mod
    base = inv(4 * FOUR_FACTORIAL ** 2, mod) * i4 * a * a
    delta8 = (
        inv(4 * FOUR_FACTORIAL ** 2, mod)
        * 8 * b_squared * inv(i4, mod) * a * i8
    )
    delta12 = (
        inv(4 * FOUR_FACTORIAL ** 2, mod)
        * 16 * b_squared * b_squared
        * inv(i4, mod) * inv(i4, mod) * i12
    )
    return {
        "stress_square": int((base + delta8 + delta12) % mod),
        "root_term": int(base % mod),
        "delta_I8": int(delta8 % mod),
        "delta_I12": int(delta12 % mod),
        **invariants,
    }
