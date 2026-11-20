# Polynomials whose coefficients contain parameters (Prodigy #444): get_terms lists each state once, with its full
# coefficient, and the term-wise updates keep a coefficient that is a sum together.
import pytest

import prodigy_ginac as pg

PARAMS = ["p", "q"]

# gf -> variables
POLYNOMIALS = {
    "p*x + (1-p)": {"x"},
    "p*x*y + (1-p)/2*x^2*y^4 + (1-p)/2": {"x", "y"},
    "(p*x + 1 - p)*(q*y + 1 - q)": {"x", "y"},
    "p*q*x + p*(1-q) + (1-p)*x": {"x"},
    "p*x + (1-p)*x^3": {"x"},
    "p*x*w + (1-p)*w^2": {"x", "w"},
    "1/3*x + 1/3*x^2 + 1/3": {"x"},
    "p + 1 - p": {"x"},
}


def dist(gf):
    return pg.Dist(gf, PARAMS)


def equal(a, b):
    return dist(f"({a}) - ({b})").is_zero() == pg.troolean.true


@pytest.mark.parametrize("gf, variables", POLYNOMIALS.items())
def test_get_terms_one_item_per_state(gf, variables):
    terms = dist(gf).get_terms(variables)
    states = [tuple(sorted(vals.items())) for _, vals in terms]
    assert len(states) == len(set(states)), terms
    assert all(not equal(prob, "0") for prob, _ in terms), terms
    rebuilt = " + ".join(f"({prob})" + "".join(f"*{var}^{val}" for var, val in vals.items())
                         for prob, vals in terms) or "0"
    assert equal(rebuilt, gf), terms


CHOICE = "p*x + (1-p)/2*x^2 + (1-p)/2"


@pytest.mark.parametrize("update, expected", [
    (lambda d: d.update_product("y", "x", "x", {"x", "y"}, True, None), "p*x*y + (1-p)/2*x^2*y^4 + (1-p)/2"),
    (lambda d: d.update_modulo("y", "x", "2", {"x", "y"}, True, None), "p*x*y + (1-p)/2*x^2 + (1-p)/2"),
    (lambda d: d.update_modulo("y", "2", "3", {"x", "y"}, True, None), "(p*x + (1-p)/2*x^2 + (1-p)/2)*y^2"),
    (lambda d: d.update_power("y", "x", "2", None), "p*x*y + (1-p)/2*x^2*y^4 + (1-p)/2"),
])
def test_term_wise_updates_with_parameters(update, expected):
    assert equal(update(dist(CHOICE)), expected)


def test_approximated_marginal_counts_zero_once():
    # the approximation 1 - p of (1-p)/(1-p*x) is one state x = 0, not the two terms 1 and -p
    geometric = dist("(1-p)/(1-p*x)")
    assert equal(geometric.update_product("y", "x", "x", {"x", "y"}, False, "0.5"), "1 - p")
