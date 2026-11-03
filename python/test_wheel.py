# Checks of an installed prodigy_ginac (cibuildwheel runs them on every wheel).
import prodigy_ginac as pg


def test_version():
    assert pg.__version__ != "dev"


def test_dist_from_string():
    dist = pg.Dist("1/(2-x)")
    assert pg.find_symbols(str(dist)) == ["x"]
    assert dist.is_polynomial({"x"}) == pg.troolean.false


def test_geometric_expected_value():
    dist = pg.geometric("x", "1/2")
    assert dist.E("x") == "1"
    assert dist.is_zero() == pg.troolean.false
