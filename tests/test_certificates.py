from chiral4form.certificates import build_quotient_certificate


def test_quotient_certificate_constructs_annihilators():
    p = 101
    rows = [[1, 0, 1, 0], [0, 1, 0, 1]]
    cert = build_quotient_certificate(rows, ambient_dim=4, prime=p)
    assert cert.reachable_rank == 2
    assert cert.quotient_dim == 2
    assert len(cert.annihilators) == 2
