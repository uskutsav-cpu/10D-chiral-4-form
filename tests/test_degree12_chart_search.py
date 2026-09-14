from chiral4form.degree12_chart_search import candidate_charts,CERTIFIED_FREE

def test_certified_chart_first():
    it=iter(candidate_charts())
    assert next(it)==CERTIFIED_FREE

def test_candidate_charts_have_four_distinct_leading_columns():
    for i,c in zip(range(100),candidate_charts()):
        assert len(c)==4
        assert len(set(c))==4
        assert all(0<=x<72 for x in c)
