from chiral4form.degree12_all_orders import degree12_lift_readiness

def test_degree12_gate_never_confuses_rank_certificate_with_ideal_tangency():
    r=degree12_lift_readiness()
    assert r["finite_characteristic_zero_local_theorem"] is True
    assert r["finite_global_reachability_theorem"] is True
    if not r["explicit_exact_QQ_orbit_ideal_present"]:
        assert r["all_orders_lift_ready"] is False
