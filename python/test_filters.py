# The filters truncate the gf with GiNaC's series, which does not truncate products of polynomials (issue #6).
import pytest

import prodigy_ginac as pg

B = "(1/3+2/3*b)"

# gf in factored form -> its coefficients of c^k, as a dict for finitely many terms or as a function of k
GFS = {
    "(1+c)*(1+c^5)": {0: "1", 1: "1", 5: "1", 6: "1"},
    "(1+c*a)*(1+c^5)": {0: "1", 1: "a", 5: "1", 6: "a"},
    "1/4*(1+c)*(1+c^3)*(1+c^5)": {k: "1/4" for k in [0, 1, 3, 4, 5, 6, 8, 9]},
    f"(1/2+1/2*c*a)*(1/2+1/2*c^5)*{B}": {0: f"1/4*{B}", 1: f"1/4*a*{B}", 5: f"1/4*{B}", 6: f"1/4*a*{B}"},
    "(p+(1-p)*c^2)*(1+c^5)": {0: "p", 2: "(1-p)", 5: "p", 7: "(1-p)"},
    "(1+c*a)*(1+c^5)/(2-a)": {0: "1/(2-a)", 1: "a/(2-a)", 5: "1/(2-a)", 6: "a/(2-a)"},
    # infinite in c
    "(1+c^5)/(2-c)": lambda k: f"1/2^{k + 1}" + (f"+1/2^{k - 4}" if k >= 5 else ""),
    "(1+c*a)*(1+c^5)/(2-c)": lambda k: f"1/2^{k + 1}" + (f"+a/2^{k}" if k >= 1 else "")
                                       + (f"+1/2^{k - 4}" if k >= 5 else "") + (f"+a/2^{k - 5}" if k >= 6 else ""),
}

KEEP = {
    "filterLess": lambda k, n: k < n,
    "filterLeq": lambda k, n: k <= n,
    "filterGreater": lambda k, n: k > n,
    "filterGeq": lambda k, n: k >= n,
    "filterEq": lambda k, n: k == n,
    "filterNeq": lambda k, n: k != n,
}


def dist(gf):
    return pg.Dist(gf, ["p"]) if "p" in gf else pg.Dist(gf)


def term(gf, k):
    coefficients = GFS[gf]
    coefficient = coefficients(k) if callable(coefficients) else coefficients.get(k, "0")
    return dist(f"({coefficient})*c^{k}")


def expected(gf, filter_name, n):
    keep = KEEP[filter_name]
    if keep(n + 1, n):
        # infinitely many terms are kept: the gf minus the dropped terms, which all have degree <= n
        result = dist(gf)
        for k in range(n + 1):
            if not keep(k, n):
                result = result - term(gf, k)
        return result
    result = pg.Dist("0")
    for k in range(n + 1):
        if keep(k, n):
            result = result + term(gf, k)
    return result


@pytest.mark.parametrize("gf", GFS)
@pytest.mark.parametrize("n", [0, 1, 3, 5, 6, 7, 12])
@pytest.mark.parametrize("filter_name", KEEP)
def test_filter_of_product(gf, n, filter_name):
    result = getattr(dist(gf), filter_name)("c", str(n))
    assert result == expected(gf, filter_name, n), f"{filter_name}({gf}, c, {n}) = {result}"


def test_issue_examples():
    d = pg.Dist("(1+c*a)*(1+c^5)")
    assert d.filterLess("c", "6") == pg.Dist("1+a*c+c^5")
    assert d.filterGreater("c", "5") == pg.Dist("a*c^6")
    assert d.filterEq("c", "6") == pg.Dist("a*c^6")
    assert pg.Dist("(1+c)*(1+c^5)").filterGreater("c", "5") == pg.Dist("c^6")


@pytest.mark.parametrize("gf", GFS)
def test_coefficient_iterator_of_product(gf):
    it = dist(gf).coefficient_iterator("c")
    for k in range(10):
        coefficient = it.next()
        assert coefficient * pg.Dist(f"c^{k}") == term(gf, k), f"coefficient {k} of {gf} = {coefficient}"


@pytest.mark.parametrize("gf", ["1/8*(1+x)*(1+x^5)*(y+y^2)", "1/8*(1+x*a)*(1+x^5)*(1+y^3)*y"])
def test_infinite_paths_of_updates_on_products(gf):
    # with finite=False, update_modulo and update_product split the gf with filterEq; the finite path does not
    d = pg.Dist(gf)
    variables = set(d.get_symbols())
    for update in [d.update_modulo, d.update_product]:
        assert update("z", "x", "y", variables, False, None) == update("z", "x", "y", variables, True, None)
