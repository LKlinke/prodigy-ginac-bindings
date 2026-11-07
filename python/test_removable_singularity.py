# Updates evaluate the gf at x = 1, where it may have a removable singularity (issue #20).
import pytest

import prodigy_ginac as pg

# x ~ uniform{1, 2, 3}, y = 1, written with the factor (x - 1) not cancelled
GF = "x*(x^3-1)/(3*(x-1))*y"


def test_assign_variable():
    assert pg.Dist(GF).update_var("x", "y") == pg.Dist("x*y")


def test_assign_literal():
    assert pg.Dist(GF).update_var("x", "5") == pg.Dist("x^5*y")


def test_assign_other_variable():
    assert pg.Dist(GF).update_var("y", "x") == pg.Dist("(x*y + x^2*y^2 + x^3*y^3)/3")


def test_mass():
    assert pg.Dist(GF).mass() == "1"


def test_genuine_pole_still_raises():
    with pytest.raises(ValueError):
        pg.Dist("x/(1-x)").update_var("x", "5")
