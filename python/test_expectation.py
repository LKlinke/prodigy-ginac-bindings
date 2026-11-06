# Expected values of polynomial expressions in the variables (issue #16), compared with hand-computed values.
import pytest

import prodigy_ginac as pg


def test_sums_and_linear_combinations():
    # x ~ uniform{0, 1, 2}, y ~ Bernoulli(1/4), independent
    dist = pg.Dist("(1+x+x^2)/3 * (3/4+y/4)")
    assert dist.E("x + y") == "5/4"
    assert dist.E("2*x - 3*y") == "5/4"
    assert dist.E("x/2 + 1/3*y + 7") == "91/12"
    assert dist.E("0.5*x + y") == "3/4"
    assert dist.E("5") == "5"


def test_dependent_variables():
    # (x, y) in {(0, 0), (1, 2), (2, 1)} with probabilities 1/2, 1/4, 1/4
    dist = pg.Dist("1/2 + 1/4*x*y^2 + 1/4*x^2*y")
    assert dist.E("x + y") == "3/2"
    assert dist.E("x*y") == "1"
    assert dist.E("(x - y)^2") == "1/2"


def test_infinite_support():
    # x ~ geometric(1/2), y ~ geometric(1/3): E[x] = 1, E[y] = 2, E[x^2] = 3, E[y^2] = 10
    dist = pg.geometric("x", "1/2") * pg.geometric("y", "1/3")
    assert dist.E("x + y") == "3"
    assert dist.E("3*x - y/2") == "2"
    assert dist.E("(x + y)^2") == "17"


def test_parameters():
    # x ~ geometric(p): E[x] = (1-p)/p
    dist = pg.geometric("x", "p") * pg.Dist("1/2 + y/2")
    assert pg.Dist(dist.E("x + y"), "p") == pg.Dist("(1-p)/p + 1/2", "p")
    assert pg.Dist(dist.E("p*x + 2*y"), "p") == pg.Dist("2-p", "p")


@pytest.mark.parametrize("expr", ["1/(x+1)", "x^(1/2)", "exp(x)", "x/y"])
def test_unsupported_expressions(expr):
    dist = pg.geometric("x", "1/2") * pg.geometric("y", "1/2")
    with pytest.raises(NotImplementedError, match="polynomial"):
        dist.E(expr)


# issue #18: symbols that are neither variables nor parameters must raise instead of being treated as 0
def test_constants_and_variables_still_work():
    dist = pg.Dist("1/8*(x^3-x^11)/(1-x)")
    assert dist.E("x") == "13/2"
    assert dist.E("3") == "3"


@pytest.mark.parametrize("expr", ["z", "x+z", "x*z", "2*z*x + x + 1"])
def test_unknown_symbols(expr):
    dist = pg.Dist("1/8*(x^3-x^11)/(1-x)")
    with pytest.raises(ValueError, match="neither a variable nor a parameter"):
        dist.E(expr)
