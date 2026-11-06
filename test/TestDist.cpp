//
// Created by Tobias Winkler on 22.10.21.
//

#include <DistributionFactory.h>
#include "gtest/gtest.h"

using namespace std;
using namespace prodigy;
using Fac = DistributionFactory;

TEST(Dist, mass) {
    Dist dist = Fac::geometric("x", "1/2");
    ASSERT_EQ(dist.mass(), 1);
    dist = Fac::geometric("x", "p");
    ASSERT_EQ(dist.mass(), 1);
}

TEST(Dist, marginal) {
    Dist dist1 = Fac::geometric("x", "1/2");
    Dist dist2 = Fac::geometric("y", "1/2");
    Dist dist = dist1 * dist2;
    dist = dist.marginal("x");
    ASSERT_EQ(dist, Fac::geometric("x", "1/2"));
}

TEST(Dist, E) {
    Dist dist = Fac::geometric("x", "1/2");
    ASSERT_EQ(dist.E("x"), 1);
}

// issue #16: expected values of polynomial expressions in the variables
TEST(Dist, ESumsAndLinearCombinations) {
    // x ~ uniform{0, 1, 2}, y ~ Bernoulli(1/4), independent: E[x] = 1, E[y] = 1/4
    Dist dist("(1+x+x^2)/3 * (3/4+y/4)");
    ASSERT_EQ(dist.E("x"), 1);
    ASSERT_EQ(dist.E("y"), GiNaC::numeric(1, 4));
    ASSERT_EQ(dist.E("x+y"), GiNaC::numeric(5, 4));
    ASSERT_EQ(dist.E("2*x - 3*y"), GiNaC::numeric(5, 4));
    ASSERT_EQ(dist.E("x/2 + 1/3*y + 7"), GiNaC::numeric(91, 12));
    ASSERT_EQ(dist.E("0.5*x + y"), GiNaC::numeric(3, 4));
    ASSERT_EQ(dist.E("x - x"), 0);
    ASSERT_EQ(dist.E("5"), 5);
    ASSERT_EQ(dist.E("-1/3"), GiNaC::numeric(-1, 3));
}

TEST(Dist, EDependentVariables) {
    // (x, y) in {(0, 0), (1, 2), (2, 1)} with probabilities 1/2, 1/4, 1/4
    Dist dist("1/2 + 1/4*x*y^2 + 1/4*x^2*y");
    ASSERT_EQ(dist.E("x + y"), GiNaC::numeric(3, 2));
    ASSERT_EQ(dist.E("x*y"), 1);
    ASSERT_EQ(dist.E("x^2"), GiNaC::numeric(5, 4));
    ASSERT_EQ(dist.E("(x - y)^2"), GiNaC::numeric(1, 2));
    ASSERT_EQ(dist.E("x^2*y - 3*x + 1"), GiNaC::numeric(1, 4));
}

TEST(Dist, EInfiniteSupport) {
    // x ~ geometric(1/2) and y ~ geometric(1/3) on {0, 1, ...}: E[x] = 1, E[y] = 2, E[x^2] = 3, E[y^2] = 10
    Dist dist = Fac::geometric("x", "1/2") * Fac::geometric("y", "1/3");
    ASSERT_EQ(dist.E("x + y"), 3);
    ASSERT_EQ(dist.E("3*x - y/2"), 2);
    ASSERT_EQ(dist.E("x*y"), 2);
    ASSERT_EQ(dist.E("x^2 + y^2"), 13);
    ASSERT_EQ(dist.E("(x + y)^2"), 17);
    // the variables are counted in the marginal only: a third variable z does not matter
    Dist withZ = dist * Dist("1/2 + z/2");
    ASSERT_EQ(withZ.E("x + y"), 3);
    ASSERT_EQ(withZ.E("x + y + 2*z"), 4);
}

TEST(Dist, EParameters) {
    // x ~ geometric(p): E[x] = (1-p)/p, y ~ Bernoulli(q): E[y] = q
    Dist dist = Fac::geometric("x", "p") * Dist("1-q+q*y", "q");
    GiNaC::ex p, q;
    for (const auto &sym : dist.getSymbols()) {
        const std::string &name = GiNaC::ex_to<GiNaC::symbol>(sym).get_name();
        if (name == "p") p = sym;
        if (name == "q") q = sym;
    }
    ASSERT_TRUE(GiNaC::is_a<GiNaC::symbol>(p) && GiNaC::is_a<GiNaC::symbol>(q));
    EXPECT_TRUE((dist.E("x + y") - ((1 - p) / p + q)).normal().is_zero());
    EXPECT_TRUE((dist.E("p*x + 2*y") - (1 - p + 2 * q)).normal().is_zero());
    EXPECT_TRUE((dist.E("x*y") - q * (1 - p) / p).normal().is_zero());
    EXPECT_TRUE((dist.E("q + 1")  - (q + 1)).normal().is_zero());
}

// issue #18: symbols that are neither variables nor parameters have no known value, E must not treat them as 0
TEST(Dist, EUnknownSymbols) {
    // x ~ uniform{3, ..., 10}: E[x] = 13/2
    Dist dist("1/8*(x^3-x^11)/(1-x)");
    ASSERT_EQ(dist.E("x"), GiNaC::numeric(13, 2));
    ASSERT_EQ(dist.E("3"), 3);
    EXPECT_THROW(dist.E("z"), std::invalid_argument);
    EXPECT_THROW(dist.E("x + z"), std::invalid_argument);
    EXPECT_THROW(dist.E("x*z"), std::invalid_argument);
    EXPECT_THROW(dist.E("x + 2*z*x + 1"), std::invalid_argument);
    EXPECT_THROW(dist.E("1/z"), std::invalid_argument);
    // parameters are still allowed as coefficients
    Dist withParam = Fac::geometric("x", "p");
    EXPECT_NO_THROW(withParam.E("p*x + 1"));
    EXPECT_THROW(withParam.E("p*x + z"), std::invalid_argument);
}

TEST(Dist, EUnsupportedExpressions) {
    Dist dist = Fac::geometric("x", "1/2") * Fac::geometric("y", "1/2");
    EXPECT_THROW(dist.E("1/(x+1)"), not_implemented_error);
    EXPECT_THROW(dist.E("x^(1/2)"), not_implemented_error);
    EXPECT_THROW(dist.E("exp(x)"), not_implemented_error);
    EXPECT_THROW(dist.E("x/y"), not_implemented_error);
    Dist withParam = Fac::geometric("x", "p");
    EXPECT_THROW(withParam.E("x^p"), not_implemented_error);
    EXPECT_THROW(withParam.E("p"), std::invalid_argument);
}

TEST(Dist, update) {
    Dist dist = Fac::geometric("x", "1/2");
    dist = dist.updateSum("x", "x", "5");
    ASSERT_EQ(dist.E("x"), 6);
}

TEST(Dist, updateVar) {

    Dist dist = Fac::uniform("x", "0", "4");
    

    EXPECT_ANY_THROW(dist.updateVar("x", "x+123"));

    dist = dist.updateVar("y", "x");

    dist = dist.updateVar("y", "w");
    ASSERT_EQ(dist, Dist {"1/5 * (1+w*x+(w*x)^2+(w*x)^3+(w*x)^4)"});
}

TEST(Dist, updateMonus) {
    Dist dist {"1/2 + 1/2 * x^2"};
    EXPECT_ANY_THROW(dist.updateSubtraction("x", "x", "1"));
    //ASSERT_EQ(dist.updateSubtraction("x", "x", "1"), Dist {"1/2 + 1/2 * x"});
}

TEST(Dist, equality) {
    Dist dist {"(x^2)^(9/2)"}, dist2 {"x^9"};
    ASSERT_EQ(dist, dist2);
}

TEST(Dist, iid) {
    Dist  dist = Dist {"1/2*x + 1/2*x^2"} * Dist {"y^2"};
    dist = dist.updateIid("res", Fac::geometric("T","1/3"), "y");
    ASSERT_EQ(dist, Dist {"(1/2 * x + 1/2 * x^2) * y^2"} * Fac::geometric("res", "1/3") * Fac::geometric("res", "1/3"));

    dist = Dist {"x"};
    dist = dist.updateIid("x", Fac::geometric("T", "p"), "x");
    ASSERT_EQ(dist, Fac::geometric("x", "p"));
}

TEST(Dist, filter) {
    Dist dist = Fac::geometric("x", "1/2");
    ASSERT_EQ(dist.filterLess("x", "3").mass(), GiNaC::numeric(7,8));
    ASSERT_EQ(dist.filterLeq("x", "2").mass(), GiNaC::numeric(7,8));
}

TEST(Dist, filterUnknownVar) {
    Dist dist {"1"};
    ASSERT_TRUE(dist.isTrivial() == troolean::TRUE);
    ASSERT_EQ(dist.filterLess("x", "1").mass(), GiNaC::numeric(1));
}

TEST(Dist, gfSize) {
    ASSERT_EQ(Dist("(x+1)^2").gfSize(), Dist("(x+1)^100").gfSize());
    ASSERT_EQ(Dist("x").gfSize(), 1);
    ASSERT_EQ(Dist("x+y").gfSize(), 3);
}

TEST(Dist, coefficientIterator) {
    Dist dist = Fac::geometric("x", "q");
    Dist::CoefficientIterator it = dist.coefficientIterator("x");
    Dist dist2 = it.next() + it.next() + it.next();
    Dist dist3("q + (1-q)*q + (1-q)^2 * q", "q");
    EXPECT_ANY_THROW(dist.coefficientIterator("y"));
    ASSERT_EQ(dist2, dist3);
}
// r - Dist(str(r)) is zero only if the gf of r uses the registered symbols (issue #3)
static void expectLinked(const Dist& r) {
    std::stringstream s;
    s << r;
    std::vector<std::string> params;
    for (const auto& sym : r.getSymbols()) {
        const auto& name = GiNaC::ex_to<GiNaC::symbol>(sym).get_name();
        if (r.isKnownAsParam(name)) params.push_back(name);
    }
    EXPECT_TRUE((r - Dist(s.str(), params)).isZero() == troolean::TRUE) << s.str();
}

TEST(Dist, registeredSymbols) {
    for (const auto& part : Dist("1/(2-x)").arithmeticProgression("x", "2")) {
        expectLinked(part);
    }
    expectLinked(Dist("1/2*x^3+1/2*x^4").updateModulo("x", "x", "2", {"x"}, true, std::nullopt));
    expectLinked(Dist("1/2*x+1/2", "p").updateIid("y", Dist("p*n+(1-p)", "p"), "x"));
    expectLinked(Fac::dirac("x", "3"));
    expectLinked(Dist("x^2").setAssumptions({}));
    Dist::CoefficientIterator it = Fac::geometric("x", "1/2").coefficientIterator("x");
    it.next();
    expectLinked(it.rest());
}

// variables are possymbols, parameters are plain symbols (issue #5), also in results
static void expectPossymbols(const Dist& r) {
    for (auto it = r.gf().preorder_begin(); it != r.gf().preorder_end(); ++it) {
        if (GiNaC::is_a<GiNaC::symbol>(*it)) {
            const auto& name = GiNaC::ex_to<GiNaC::symbol>(*it).get_name();
            EXPECT_EQ(GiNaC::is_a<GiNaC::possymbol>(*it), r.isKnownAsVar(name)) << *it << " in " << r;
        }
    }
    for (const auto& sym : r.getSymbols()) {
        const auto& name = GiNaC::ex_to<GiNaC::symbol>(sym).get_name();
        EXPECT_EQ(GiNaC::is_a<GiNaC::possymbol>(sym), r.isKnownAsVar(name)) << sym;
    }
}

TEST(Dist, possymbols) {
    expectPossymbols(Dist("x*y^2"));
    expectPossymbols(Dist("p*x + (1-p)", "p"));
    expectPossymbols(Fac::geometric("x", "p"));
    expectPossymbols(Fac::uniform("x", "a", "b"));
    expectPossymbols(Fac::dirac("x", "n"));
    expectPossymbols(Fac::catalan("x", "p"));
    Dist dist = Fac::geometric("x", "p");
    expectPossymbols(dist.updateSum("y", "x", "z"));
    expectPossymbols(dist.updateIid("y", Dist("q*n + (1-q)", "q"), "x"));
    expectPossymbols(dist.coefficientIterator("x").next());
    expectPossymbols(dist * "r");
    expectPossymbols(dist.substituteParam("p", "s"));
    Dist fresh = dist;
    expectPossymbols(dist.updateVar(fresh.getFreshVariable(), "x"));
    expectPossymbols(dist.updateRoot("r", "4", "2", std::nullopt));
}

// simplifications that hold for positive variables are applied eagerly, before any substitution (issue #5)
TEST(Dist, positiveVariables) {
    Dist dist("sqrt(x^2)");
    GiNaC::ex x = dist.getSymbols().at(0);
    ASSERT_TRUE(dist.gf().is_equal(x));
    EXPECT_TRUE(dist.gf().subs(x == -x).is_equal(-x));
    EXPECT_TRUE(dist.gf().subs(x == GiNaC::I * x).is_equal(GiNaC::I * x));
    // parameters may be negative
    EXPECT_FALSE(Dist("sqrt(a^2)", "a") == Dist("a", "a"));
}

// GiNaC's series does not truncate products of polynomials (issue #6)
TEST(Dist, filterProductOfPolynomials) {
    Dist dist("(1+c*a)*(1+c^5)");
    EXPECT_EQ(dist.filterLess("c", "6"), Dist("1+a*c+c^5"));
    EXPECT_EQ(dist.filterGreater("c", "5"), Dist("a*c^6"));
    EXPECT_EQ(dist.filterEq("c", "6"), Dist("a*c^6"));
    EXPECT_EQ(Dist("(1+c)*(1+c^5)").filterGreater("c", "5"), Dist("c^6"));
}
