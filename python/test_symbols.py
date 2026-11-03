# Every operation that returns a Dist must link its gf to the symbol table: a variable or parameter is the same GiNaC
# symbol wherever its name occurs. Otherwise equal terms do not cancel and r - Dist(str(r)) is not zero (issue #3).
import pytest

import prodigy_ginac as pg


def assert_linked(r):
    params = [s for s in r.get_symbols() if r.is_known_as_param(s)]
    assert (r - pg.Dist(str(r), params)).is_zero() == pg.troolean.true, str(r)


def coefficients():
    it = pg.geometric("x", "p").coefficient_iterator("x")
    return [it.next(), it.next(), it.rest()]


INFINITE = pg.Dist("1/(2-x)")
FINITE = pg.Dist("1/2*x^3+1/2*x^4")
FINITE_XY = pg.Dist("1/4*x^2*y + 1/4*x^3*y^2 + 1/2*x^4*y^2")
PARAMETRIC = pg.Dist("1/2*x+1/2", ["p"])

CASES = {
    "Dist": lambda: [pg.Dist("x*y^2"), pg.Dist("p*x + (1-p)", ["p"]), pg.Dist("p*x + (1-p)", "p")],
    "geometric": lambda: [pg.geometric("x", "1/2"), pg.geometric("x", "p")],
    "uniform": lambda: [pg.uniform("x", "0", "3"), pg.uniform("x", "a", "b")],
    "dirac": lambda: [pg.dirac("x", "3"), pg.dirac("x", "n")],
    "catalan": lambda: [pg.catalan("x", "1/2"), pg.catalan("x", "p")],
    "normalize": lambda: [pg.Dist("1/4*x + 1/4").normalize()],
    "marginal": lambda: [FINITE_XY.marginal("x"), FINITE_XY.marginal(["y"])],
    "arithmetic_progression": lambda: INFINITE.arithmetic_progression("x", "2")
        + FINITE.arithmetic_progression("x", "3"),
    "updateIid": lambda: [
        PARAMETRIC.updateIid("y", pg.Dist("p*n+(1-p)", ["p"]), "x"),
        PARAMETRIC.updateIid("x", pg.Dist("q*n+(1-q)", ["q"]), "x"),
        pg.Dist("x").updateIid("y", pg.geometric("T", "p"), "x"),
        pg.Dist("x").updateIid("y", pg.geometric("x", "1/2"), "2"),
    ],
    "filters": lambda: [
        INFINITE.filterLess("x", "2"), INFINITE.filterLeq("x", "2"), INFINITE.filterEq("x", "2"),
        INFINITE.filterNeq("x", "2"), INFINITE.filterGreater("x", "2"), INFINITE.filterGeq("x", "2"),
        INFINITE.filterLess("x", "0"),
    ],
    "update_var": lambda: [FINITE_XY.update_var("x", "y"), FINITE_XY.update_var("z", "x"),
                           FINITE_XY.update_var("x", "5")],
    "update_sum": lambda: [FINITE_XY.update_sum("z", "x", "y"), FINITE_XY.update_sum("x", "x", "2"),
                           PARAMETRIC.update_sum("y", "x", "p"), FINITE_XY.update_sum("z", "1", "2")],
    "update_product": lambda: [
        FINITE_XY.update_product("z", "x", "y", {"x", "y"}, True, None),
        INFINITE.update_product("y", "x", "x", {"x"}, False, "0.9"),
        FINITE_XY.update_product("z", "x", "3", {"x", "y"}, True, None),
        FINITE_XY.update_product("z", "2", "3", {"x", "y"}, True, None),
    ],
    "update_subtraction": lambda: [FINITE_XY.update_subtraction("z", "x", "y"),
                                   FINITE_XY.update_subtraction("x", "x", "1")],
    "update_modulo": lambda: [
        FINITE.update_modulo("x", "x", "2", {"x"}, True, None),
        FINITE_XY.update_modulo("z", "x", "y", {"x", "y"}, None, None),
        INFINITE.update_modulo("y", "x", "2", {"x"}, False, None),
        (INFINITE * pg.Dist("1/2*y + 1/2*y^2")).update_modulo("z", "x", "y", {"x", "y"}, False, None),
        FINITE.update_modulo("y", "7", "3", {"x"}, None, None),
    ],
    "update_division": lambda: [pg.Dist("1/2*x^4*y^2 + 1/2*x^3*y").update_division("z", "x", "y", None),
                                pg.Dist("1/2*x^2 + 1/2*x^4").update_division("y", "x", "2", None)],
    "update_power": lambda: [FINITE_XY.update_power("z", "x", "y", None), FINITE_XY.update_power("z", "x", "2", None),
                             FINITE_XY.update_power("z", "2", "y", None)],
    "update_root": lambda: [pg.Dist("1/2*x^4 + 1/2*x^9").update_root("y", "x", "2", None)],
    "approximate_unilaterally": lambda: [INFINITE.approximate_unilaterally("x", "0.9")],
    "set_assumptions": lambda: [FINITE.set_assumptions([])],
    "operators": lambda: [
        INFINITE * pg.geometric("x", "1/2"), INFINITE + pg.Dist("x"), INFINITE - pg.Dist("x"),
        PARAMETRIC * pg.Dist("p*y + (1-p)", ["p"]), PARAMETRIC - pg.Dist("1/2*x+1/2", ["p"]),
    ],
    "coefficient_iterator": coefficients,
}


@pytest.mark.parametrize("operation", CASES)
def test_results_use_registered_symbols(operation):
    for r in CASES[operation]():
        assert_linked(r)


def test_iid_shares_parameter_symbol():
    r = PARAMETRIC.updateIid("y", pg.Dist("p*n+(1-p)", ["p"]), "x")
    # with two different symbols p, this would not cancel
    assert (r - pg.Dist("1/2*x*(p*y+(1-p))+1/2", ["p"])).is_zero() == pg.troolean.true
    assert r.get_symbols() == ["p", "x", "y"]


def test_iid_parameter_known_as_variable():
    with pytest.raises(Exception, match="known as variable"):
        pg.Dist("p*x").updateIid("y", pg.Dist("p*n+(1-p)", ["p"]), "x")


def test_variables_are_positive():
    # issue #5: (x^2)^(9/2) = x^9 and sqrt(x^2) = x hold for the power series in a variable x
    assert pg.Dist("(x^2)^(9/2)") == pg.Dist("x^9")
    assert pg.Dist("sqrt(x^2)") == pg.Dist("x")
    assert pg.dirac("x", "n") * pg.Dist("sqrt(y^2)") == pg.Dist("x^n*y", ["n"])
    # the simplification happens before x -> -x in the arithmetic progression: odd part x, even part 0
    even, odd = pg.Dist("sqrt(x^2)").arithmetic_progression("x", "2")
    assert even.is_zero() == pg.troolean.true
    assert odd == pg.Dist("x")


def test_parameters_may_be_negative():
    assert (pg.Dist("sqrt(a^2)", ["a"]) - pg.Dist("a", ["a"])).is_zero() != pg.troolean.true
    assert (pg.geometric("x", "a") * pg.Dist("sqrt(a^2)", ["a"]) - pg.geometric("x", "a") * pg.Dist("a", ["a"])) \
        .is_zero() != pg.troolean.true


def test_fresh_variable_is_a_variable():
    # a fresh variable (a possymbol) was sliced into a plain symbol by getSymbols, so isFiniteSupport took this
    # infinite gf for finite and update_modulo failed with "can only iterate over terms of a finite polynomial"
    d = pg.geometric("x", "1/2")
    v = d.get_fresh_variable()
    d = d.update_var(v, "x").update_var("x", "0")
    assert d.update_modulo("z", v, "2", {"x", v}, None, None) == \
        pg.Dist(f"-(2+{v}*z)*(-2+{v})^(-1)*(2+{v})^(-1)")
