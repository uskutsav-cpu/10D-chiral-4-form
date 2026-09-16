from chiral4form.degree12_explicit_ideal import ansatz_polynomials, LOWER_LABELS


def test_degree12_ansatz_is_81_dimensional():
    basis, labels = ansatz_polynomials()
    assert len(basis) == 81
    assert len(labels) == 81
    assert len(LOWER_LABELS) == 9


def test_degree12_ansatz_has_72_new_coordinates():
    _, labels = ansatz_polynomials()
    assert labels[0] == "Z1"
    assert labels[71] == "Z72"
    assert tuple(labels[-9:]) == tuple(LOWER_LABELS)
