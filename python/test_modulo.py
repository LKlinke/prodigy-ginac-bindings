# Sections modulo q of rational gfs (arithmetic progressions, issue #10) are rational functions with rational
# coefficients, without roots of unity. Their coefficients are compared with filterEq, which truncates the series of
# the original gf and does not use the multisection.
import pytest

import prodigy_ginac as pg

GFS = [
    ("1/2/(1-1/2*x)", []),
    ("(1+x^2)/((2-x)*(3-x))", []),
    ("x^4/(2-x)^3", []),
    ("(1/4+x/4)^3*(1+x^7)", []),
    ("1/(3-x-y)", []),
    ("y*(1+x*y)/((2-x*y)*(4-x-y^2))", []),
    ("p/(1-(1-p)*x)", ["p"]),
    ("(1-p)*(1+a*x)/(1-p*x^2*y)", ["p", "a"]),
    ("exp(-l)/(2-x)", ["l"]),
]


@pytest.mark.parametrize("gf,params", GFS)
@pytest.mark.parametrize("q", [2, 3, 5, 6])
def test_sections(gf, params, q):
    d = pg.Dist(gf, params)
    sections = d.arithmetic_progression("x", str(q))
    assert len(sections) == q
    total = pg.Dist("0", params)
    for r, section in enumerate(sections):
        s = str(section)
        assert "Pi" not in s and "I" not in s, f"{gf}, q = {q}, r = {r}: {s}"
        for k in range(2 * q + 2):
            expected = d.filterEq("x", str(k)) if k % q == r else pg.Dist("0", params)
            assert section.filterEq("x", str(k)) == expected, f"{gf}, q = {q}, r = {r}, k = {k}: {s}"
        total = total + section
    assert total == d


@pytest.mark.parametrize("q", [2, 3, 5, 6])
def test_geometric(q):
    # Prodigy#35: 1/2/(1-1/2*x) filtered with x % q = r
    for r, section in enumerate(pg.Dist("1/2/(1-1/2*x)").arithmetic_progression("x", str(q))):
        assert section == pg.Dist(f"x^{r}/2^{r + 1}/(1-x^{q}/2^{q})")


def test_update_modulo():
    d = pg.Dist("1/2/(1-1/2*x)")
    res = d.update_modulo("y", "x", "3", {"x", "y"}, False, None)
    assert res == pg.Dist("(1+x*y/2+x^2*y^2/4)/(2-x^3/4)")
    assert "Pi" not in str(res)
