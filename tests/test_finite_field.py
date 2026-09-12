from chiral4form.finite_field import annihilator_basis, matvec, matrix_rank, nullspace, span_contains


def test_rank_and_membership():
    p = 101
    rows = [[1, 0, 1], [0, 1, 1], [1, 1, 2]]
    assert matrix_rank(rows, p) == 2
    assert span_contains(rows[:2], rows[2], p)
    assert not span_contains(rows[:1], rows[1], p)


def test_nullspace():
    p = 101
    rows = [[1, 2, 3], [0, 1, 1]]
    ns = nullspace(rows, p)
    assert len(ns) == 1
    assert matvec(rows, ns[0], p) == [0, 0]


def test_annihilator_is_right_nullspace():
    p = 103
    reachable = [[1, 0, 1, 0], [0, 1, 0, 1]]
    anns = annihilator_basis(reachable, p)
    assert len(anns) == 2
    assert all(matvec(reachable, ell, p) == [0, 0] for ell in anns)
