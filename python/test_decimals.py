# Decimal literals are parsed as exact rationals; GiNaC's parser alone makes them floating-point numbers (issue #9).
import pytest

import prodigy_ginac as pg


def test_issue_9_examples():
    assert not pg.Dist("1.00000000000000000001*x") == pg.Dist("x")
    assert not pg.Dist("(1e20+1)*x") == pg.Dist("1e20*x")
    assert pg.Dist("(1e20+1)*x") == pg.Dist("100000000000000000001*x")


@pytest.mark.parametrize("decimal, exact", [
    ("0.1*x", "1/10*x"),
    (".5*x", "1/2*x"),
    ("2.5e-3*x", "1/400*x"),
    ("1e20*x", "100000000000000000000*x"),
    ("3e0*x", "3*x"),
    ("0.1*x + 0.2*x", "3/10*x"),
    ("0.5*x1 + 0.5*e", "1/2*x1 + 1/2*e"),
])
def test_exact(decimal, exact):
    assert pg.Dist(decimal) == pg.Dist(exact)
    assert (pg.Dist(decimal) - pg.Dist(exact)).is_zero() == pg.troolean.true


def test_geometric_invariant():
    # the wrong EVT invariant of Prodigy#117 must not equal the correct one
    wrong = pg.Dist("(x + 1/2)/(1 - 0.50000000000000000001*c)")
    right = pg.Dist("(x + 1/2)/(1 - 1/2*c)")
    assert not wrong == right
    assert (wrong - right).is_zero() == pg.troolean.false


def test_malformed():
    with pytest.raises(ValueError):
        pg.Dist("1.2.3*x")
