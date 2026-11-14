#!/usr/bin/env python3
"""Independent exact algebra replay. Imports no chiral4form/foundation modules.

Uses SymPy QQ arithmetic on committed polynomial data. This checks reduced
models, not the physics-to-model derivation or global polynomial-ring claims.
"""
from __future__ import annotations
import argparse
from copy import deepcopy
from functools import lru_cache
import hashlib
import json
from pathlib import Path
import sys
import sympy as s

ROOT = Path(__file__).resolve().parents[1]


def require(ok, message):
    if not ok:
        raise ValueError(message)


def read(root, name):
    return json.loads((root / name).read_text())


def polynomial(data, variables):
    require(data['nvars'] == len(variables), 'polynomial dimension mismatch')
    require(data.get('prime') is None, 'exact QQ polynomial required')
    out = s.S.Zero
    for term in data['terms']:
        powers = term['powers']
        require(len(powers) == len(variables) and all(type(e) is int and e >= 0 for e in powers), 'invalid exponent')
        out += s.Rational(term['coefficient']) * s.prod(v**e for v, e in zip(variables, powers))
    return s.expand(out)


def deriv(q, field, variables):
    return s.expand(sum(s.diff(q, v) * f for v, f in zip(variables, field) if f))


def degree10(root):
    model = read(root, 'verification/degree10/rational_model.json')
    obs = read(root, 'verification/degree10/obstruction_certificate.json')
    require(model['coordinate_ids'] == obs['coordinate_ids'], 'degree-ten coordinate mismatch')
    require(len(model['coordinate_ids']) == 24, 'degree-ten ambient dimension mismatch')
    x = s.symbols('x0:24')
    fields = {name: [polynomial(q, x) for q in field] for name, field in model['fields'].items()}
    A, B, C, D, U, V = base = s.symbols('A B C D U V')
    lower = [A**4, A*A*B, A*C, A*D, B*B]
    rows = s.Matrix(obs['obstructions_integer_coordinates'])
    require(rows.shape == (12, 19), 'degree-ten obstruction shape mismatch')
    dep = [j for j in range(14) if j not in (0, 12)]
    rhs = -rows[:, 0] * U - rows[:, 12] * V - rows[:, 14:] * s.Matrix(lower)
    dependent = rows[:, dep].inv() * rhs
    z = [None] * 14
    z[0], z[12] = U, V
    for j, value in zip(dep, dependent):
        z[j] = s.expand(value)
    graph = [A, B, 0, C, -16*A**3, 0, 0, 0, 0, D] + z
    subs = dict(zip(x, graph))
    lower_constraints = [x[2], x[5], x[6], x[7], x[8], x[4] + 16*x[0]**3]
    candidates = list(x[10:]) + [x[0]**4, x[0]**2*x[1], x[0]*x[3], x[0]*x[9], x[1]**2]
    constraints = lower_constraints + list(rows * s.Matrix(candidates))
    for q in constraints:
        require(s.expand(q.subs(subs, simultaneous=True)) == 0, 'degree-ten constraint does not vanish on graph')
        for name, field in fields.items():
            require(s.expand(deriv(q, field, x).subs(subs, simultaneous=True)) == 0, f'degree-ten tangency failure: {name}')
    reduced = {name: [s.expand(field[j].subs(subs, simultaneous=True)) for j in (0, 1, 3, 9, 10, 22)] for name, field in fields.items()}
    return x, fields, base, graph, reduced


def check8(fields, x):
    v = x[:10]
    fs = {name: field[:10] for name, field in fields.items() if any(field[:10])}
    require(len(fs) == 7, 'degree-eight generator catalogue differs')
    a, j, k, b1, b2, b3, b4, b5, b6, p = v
    omega = b2 + 16*a**3
    qs = [k, b3, b4, b5, b6, omega]
    graph = {k: 0, b3: 0, b4: 0, b5: 0, b6: 0, b2: -16*a**3}
    for name, field in fs.items():
        require(s.expand(deriv(omega, field, v) - (60*omega if name == 'tr1' else 0)) == 0, 'Omega8 transport failure')
        for q in qs:
            require(s.expand(deriv(q, field, v).subs(graph, simultaneous=True)) == 0, 'degree-eight ideal not tangent')
    extras = (2, 4, 5, 6, 7, 8)
    chosen = ['tr2', 'tr3', 'tr4', 'tr2*tr2']
    controls = [fs[name] for name in chosen] + [[s.Integer(i == j) for i in range(10)] for j in extras]
    determinant = s.expand(s.Matrix.hstack(*[s.Matrix(f) for f in controls]).det())
    require(determinant != 0 and not determinant.free_symbols, 'degree-eight completion is not everywhere full rank')
    # Each omitted extra leaves its corresponding relation invariant.
    for omitted, q in zip(extras, (k, omega, b3, b4, b5, b6)):
        for added in extras:
            if added != omitted:
                require(s.diff(q, v[added]) == 0, 'catalogue-removal obstruction failed')
    return {'orbit_dimension': 4, 'ambient_dimension': 10, 'completion_minor': str(determinant),
            'six_catalogue_removal_obstructions': True,
            'excluded_interaction': 'a=0, b2=1 violates Omega8=0',
            'extra_generator_escape': 'adding I8_2 translates b2 and leaves the pure-stress orbit'}


def degree12(root, base, reduced, source_override=None):
    source = deepcopy(source_override) if source_override is not None else read(root, 'verification/all_orders/degree12_exact/source_equation_space.json')
    target = read(root, 'verification/all_orders/degree12_exact/exact_source_targets.json')
    kernel = read(root, 'verification/all_orders/degree12_exact/kernel_lift.json')
    ideal = read(root, 'verification/all_orders/degree12_exact/ideal.json')
    rows = source['all_source_rows']
    require(len(rows) == 36 and [r['source_row'] for r in rows] == list(range(36)), '36 ordered source rows required')
    targets = {r['source_row']: r for r in target['rows']}
    A, B, C, D, U, V = base
    lower = [A**5, A**3*B, A*A*C, A*A*D, A*B*B, A*U, A*V, B*C, B*D]
    equations = []
    Z = s.symbols('Z0:72')
    variables = tuple(base) + Z
    fields = {name: list(field) + [s.S.Zero]*72 for name, field in reduced.items()}
    fields['tr1'] = list(reduced['tr1']) + [100*z for z in Z]
    for meta in rows:
        idx = meta['source_row']
        require(len(meta['row']) == 81, 'source row width mismatch')
        row = list(map(s.Rational, meta['row']))
        require(row[:72] == list(map(s.Rational, targets[idx]['coordinates'])), 'source/physical-target coordinates differ')
        exp = [0]*6
        for factor in meta['output_monomial']:
            var = str(factor['variable'])
            require(var in ('A', 'B', 'C', 'D', 'U', 'V'), 'invalid source variable')
            exp[('A', 'B', 'C', 'D', 'U', 'V').index(var)] += int(factor['power'])
        monomial = s.prod(v**e for v, e in zip(base, exp))
        field = reduced.get(meta['generator'], [s.S.Zero]*6)
        tail = [s.Poly(deriv(q, field, base), *base).coeff_monomial(tuple(exp)) for q in lower]
        require(row[72:] == tail, 'independent degree-ten tail derivation differs')
        equations.append(row)
        if meta['generator'] not in fields:
            fields[meta['generator']] = [s.S.Zero]*78
        for j, coefficient in enumerate(row[:72]):
            fields[meta['generator']][6+j] += coefficient*monomial
    E = s.Matrix(equations)
    K = s.Matrix([[s.Rational(q) for q in row] for row in kernel['lift']['rational_rows']])
    require(E.shape == (36, 81) and E.rank() == 13, 'degree-twelve source rank differs')
    require(K.shape == (68, 81) and K.rank() == 68, 'degree-twelve kernel rank differs')
    require(E*K.T == s.zeros(36, 68), 'degree-twelve tangency residual is nonzero')
    require(K[:, :72].rank() == 68, 'leading obstruction rank differs')
    rref, pivots = K.rref()
    require(rref == K and list(pivots) == kernel['pivots'], 'kernel RREF metadata differs')
    free = [j for j in range(72) if j not in pivots]
    require(free == kernel['free_degree12_columns'] and len(free) == 4, 'free fiber coordinates differ')
    basis = s.Matrix(list(Z)+lower)
    constraints = list(K*basis)
    require(ideal['ansatz_labels'] == kernel['ansatz_labels'], 'ideal coordinate ordering differs')
    require(len(ideal['constraints']) == 68, 'ideal constraint count differs')
    for q, record in zip(constraints, ideal['constraints']):
        require(s.expand(q-polynomial(record, variables)) == 0, 'ideal polynomial differs from exact kernel')
    vertical = ['tr2*tr2*tr2', 'tr2*tr4', 'tr3*tr3', 'tr6']
    for name in vertical:
        require(all(q == 0 for q in fields[name][:6]), 'vertical generator changes base')
        require(all(not s.sympify(q).free_symbols for q in fields[name][6:]), 'vertical action is not constant')
    translation = s.Matrix([[fields[name][6+j] for name in vertical] for j in free])
    det = translation.det()
    require(det != 0, 'vertical translation matrix singular')
    # Every constant vertical direction lies in the homogeneous fiber kernel.
    require(K[:, :72] * s.Matrix.hstack(*[s.Matrix(fields[n][6:]) for n in vertical]) == s.zeros(68, 4), 'vertical direction escapes ideal')
    return {'source_rank': 13, 'relative_constraints': 68, 'orbit_dimension': 10,
            'ambient_dimension': 96, 'vertical_free_columns': free,
            'vertical_matrix_QQ': [[str(q) for q in row] for row in translation.tolist()],
            'vertical_determinant_QQ': str(det), 'translation_matrix_scope': 'exact reduced model loaded from committed source-target coordinates'}, variables, fields, constraints, free, translation, graph_from_kernel(K, pivots, basis, variables)


def graph_from_kernel(K, pivots, basis, variables):
    graph = list(variables)
    for row, pivot in enumerate(pivots):
        graph[6+pivot] = s.expand(-sum(K[row, j]*basis[j] for j in range(81) if j != pivot))
    return graph


def constructive_structure(base, variables, fields):
    """Symbolic conditions establishing the endpoint algorithm for all targets."""
    first = ('tr2', 'tr3', 'tr4', 'tr2*tr2')
    for i, name in enumerate(first):
        require(all(q == 0 for q in fields[name][:i]), 'constructive base control changes an earlier coordinate')
        require(fields[name][i] != 0 and not s.sympify(fields[name][i]).free_symbols, 'constructive base diagonal is not constant')
    pair = ('tr2*tr3', 'tr5')
    require(all(q == 0 for name in pair for q in fields[name][:4]), 'last base controls disturb lower coordinates')
    matrix = s.Matrix([[fields[name][i] for name in pair] for i in (4, 5)])
    require(all(not s.sympify(q).free_symbols for q in matrix) and matrix.det() != 0, 'last base-control block is not constant invertible')
    weights = (20, 40, 60, 60, 80, 80) + (100,)*72
    selected = first + pair + ('tr2*tr2*tr2', 'tr2*tr4', 'tr3*tr3', 'tr6')
    for name in selected:
        for i, q in enumerate(fields[name]):
            if q == 0:
                continue
            poly = s.Poly(q, *variables)
            require(all(sum(w*e for w,e in zip(weights,m)) < weights[i] for m in poly.monoms()), 'constructive field does not strictly lower coupling weight')
    return {'base_triangular_conditions_exact': True,
            'last_base_control_determinant_QQ': str(matrix.det()),
            'all_ten_selected_fields_strictly_lower_weight': True,
            'consequence': 'global polynomial flows and constructive reachability for arbitrary real targets; exact rational CLI endpoints for rational targets'}


def step(variables, field, point, amount):
    """Exact terminating exponential of a degree-lowering polynomial field."""
    initial = dict(zip(variables, point))
    out = []
    for coordinate in variables:
        q, value = coordinate, s.S.Zero
        for order in range(7):
            value += q.subs(initial, simultaneous=True) * amount**order / s.factorial(order)
            q = deriv(q, field, variables)
            if q == 0:
                break
        require(q == 0, 'Lie series did not terminate within field-degree bound')
        out.append(s.cancel(value))
    return out


def endpoint(variables, fields, constraints, free, translation, target):
    require(len(target) == 10, 'endpoint requires six base and four fiber targets')
    target = list(map(s.Rational, target))
    point = [s.S.Zero]*78
    segments = []
    for index, name in enumerate(('tr2', 'tr3', 'tr4', 'tr2*tr2')):
        speed = fields[name][index]
        require(not speed.free_symbols and speed != 0, 'base control speed not constant')
        amount = (target[index]-point[index])/speed
        point = step(variables, fields[name], point, amount)
        segments.append({'generator': name, 'amount': str(amount)})
    pair = ('tr2*tr3', 'tr5')
    speeds = s.Matrix([[fields[name][j] for name in pair] for j in (4, 5)])
    require(all(not q.free_symbols for q in speeds), 'last base controls not constant')
    amounts = speeds.inv() * s.Matrix([target[4]-point[4], target[5]-point[5]])
    for name, amount in zip(pair, amounts):
        point = step(variables, fields[name], point, amount)
        segments.append({'generator': name, 'amount': str(amount)})
    require(point[:6] == target[:6], 'six-stage base endpoint failed')
    shifts = translation.inv()*s.Matrix([target[6+i]-point[6+j] for i, j in enumerate(free)])
    for name, amount in zip(('tr2*tr2*tr2', 'tr2*tr4', 'tr3*tr3', 'tr6'), shifts):
        point = step(variables, fields[name], point, amount)
        segments.append({'generator': name, 'amount': str(amount)})
    require(point[:6] == target[:6] and [point[6+j] for j in free] == target[6:], 'final endpoint mismatch')
    substitution = dict(zip(variables, point))
    require(all(q.subs(substitution, simultaneous=True) == 0 for q in constraints), 'endpoint violates orbit ideal')
    return {'scope': 'exact rational endpoint in the declared finite degree-12 reduced model',
            'target': list(map(str, target)), 'segments': segments,
            'endpoint_reduced_78': list(map(str, point)), 'all_68_constraints_zero': True}


def verify(root):
    lock = read(root, 'data/combined/source-lock.json')
    for path, digest in lock['foundation_files'].items():
        require(hashlib.sha256((root/'foundation'/path).read_bytes()).hexdigest() == digest, f'foundation hash mismatch: {path}')
    require(hashlib.sha256((root/'data/combined/registry.json').read_bytes()).hexdigest() == lock['registry_sha256'], 'registry hash mismatch')
    archive = lock.get('complete_source_archive')
    if archive:
        require(hashlib.sha256((root/archive['path']).read_bytes()).hexdigest() == archive['sha256'], 'complete foundation archive changed')
    for path, target in lock.get('foundation_symlinks', {}).items():
        require((root/'foundation'/path).is_symlink() and str((root/'foundation'/path).readlink()) == target, f'foundation symlink changed: {path}')
    registry = read(root, 'data/combined/registry.json')
    dimensions = {int(d): len(ids) for d, ids in registry['degree_bases'].items()}
    character = read(root, 'foundation/science/character/independent-d5.json')
    counts = {r['degree']: int(r['invariants']) for r in character['degrees']}
    require(dimensions == {d: counts[d] for d in dimensions}, 'registry/character dimension mismatch')
    x, fields10, base, graph10, reduced = degree10(root)
    eight = check8(fields10, x)
    twelve, variables, fields, qs, free, matrix, graph12 = degree12(root, base, reduced)
    structure = constructive_structure(base, variables, fields)
    demo = endpoint(variables, fields, qs, free, matrix, ['2', '-3/7', '5', '-11/13', '1/3', '-2/5', '7/11', '-1', '4/9', '2/3'])
    report = {'schema': 1, 'passed': True, 'verification_level': 'independent exact reduced-model algebra and saved source evidence; no tensor-map recomputation',
              'homogeneous_dimensions': dimensions, 'degree8': eight,
              'degree10': {'orbit_dimension': 6, 'ambient_dimension': 24, 'constraints': 18, 'all_generators_tangent': True},
              'degree12': twelve, 'constructive_structure': structure, 'endpoint_demo': demo,
              'not_established_by_this_checker': ['correctness of tensor-to-polynomial maps', 'fresh degree-12 physical interpolation', 'physical causality/unitarity', 'global polynomial invariant-ring presentation', 'novelty']}
    return report


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--root', type=Path, default=ROOT)
    parser.add_argument('--output', type=Path)
    parser.add_argument('--target', help='six rational base targets then four rational fiber targets, comma separated')
    args = parser.parse_args()
    if args.target:
        _, _, base, _, reduced = degree10(args.root)
        _, variables, fields, qs, free, matrix, _ = degree12(args.root, base, reduced)
        report = endpoint(variables, fields, qs, free, matrix, args.target.split(','))
    else:
        report = verify(args.root)
    data = json.dumps(report, indent=2)+'\n'
    if args.output:
        args.output.parent.mkdir(parents=True, exist_ok=True)
        args.output.write_text(data)
    print(json.dumps({'passed': report.get('passed', report.get('all_68_constraints_zero')), 'output': str(args.output), 'scope': report.get('verification_level', report.get('scope'))}))


if __name__ == '__main__':
    main()
