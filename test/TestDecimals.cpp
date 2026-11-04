// Decimal literals must be parsed as exact rationals: GiNaC's parser reads them as floating-point numbers with about
// 17 digits, so different numbers compared equal (issue #9).

#include <DistributionFactory.h>
#include "gtest/gtest.h"

using namespace prodigy;
using Fac = DistributionFactory;

TEST(Decimals, issue9Examples) {
    ASSERT_FALSE(Dist("1.00000000000000000001*x") == Dist("x"));
    ASSERT_FALSE(Dist("(1e20+1)*x") == Dist("1e20*x"));
    ASSERT_TRUE(Dist("(1e20+1)*x") == Dist("100000000000000000001*x"));
    ASSERT_EQ((Dist("(1e20+1)*x") - Dist("1e20*x")).isZero(), troolean::FALSE);
}

TEST(Decimals, exactValues) {
    ASSERT_EQ(Dist("0.1*x"), Dist("1/10*x"));
    ASSERT_EQ(Dist(".5*x"), Dist("1/2*x"));
    ASSERT_EQ(Dist("5.*x"), Dist("5*x"));
    ASSERT_EQ(Dist("2.5e-3*x"), Dist("1/400*x"));
    ASSERT_EQ(Dist("2.5E-3*x"), Dist("1/400*x"));
    ASSERT_EQ(Dist("25e-4*x"), Dist("1/400*x"));
    ASSERT_EQ(Dist("1e20*x"), Dist("100000000000000000000*x"));
    ASSERT_EQ(Dist("1e+2*x"), Dist("100*x"));
    ASSERT_EQ(Dist("1.5e1*x"), Dist("15*x"));
    ASSERT_EQ(Dist("0.1*x + 0.2*x"), Dist("0.3*x"));
    ASSERT_EQ(Dist("x^0.5"), Dist("x^(1/2)"));
    ASSERT_EQ(Dist("-0.5*x+1.5"), Dist("-1/2*x+3/2"));
    ASSERT_FALSE(Dist("0.50000000000000000001*x") == Dist("0.5*x"));
}

TEST(Decimals, rewriting) {
    ASSERT_EQ(Dist::exactDecimals("0.1"), "(1/10)");
    ASSERT_EQ(Dist::exactDecimals("2.5e-3"), "(1/400)");
    ASSERT_EQ(Dist::exactDecimals(".5"), "(1/2)");
    ASSERT_EQ(Dist::exactDecimals("1e3"), "(1000)");
    ASSERT_EQ(Dist::exactDecimals("007.50"), "(15/2)");
    // integers, identifiers with digits or 'e', function names and comments are unchanged
    ASSERT_EQ(Dist::exactDecimals("12*x1 + e + exp(x2e5) + E1 + a_1e2 # 0.1"), "12*x1 + e + exp(x2e5) + E1 + a_1e2 # 0.1");
    ASSERT_EQ(Dist::exactDecimals("x1*0.25e1"), "x1*(5/2)");
}

TEST(Decimals, identifiersWithDigits) {
    Dist d("0.5*x1 + 0.5*e2");
    ASSERT_TRUE(d.isKnownAsVar("x1"));
    ASSERT_TRUE(d.isKnownAsVar("e2"));
    ASSERT_EQ(d, Dist("1/2*x1 + 1/2*e2"));
    ASSERT_EQ(d.mass(), 1);
}

TEST(Decimals, malformed) {
    ASSERT_THROW(Dist::exactDecimals("1.2.3"), std::invalid_argument);
    ASSERT_THROW(Dist::exactDecimals("1e"), std::invalid_argument);
    ASSERT_THROW(Dist::exactDecimals("1e-"), std::invalid_argument);
    ASSERT_THROW(Dist::exactDecimals(". * x"), std::invalid_argument);
    ASSERT_THROW(Dist::exactDecimals("1e1000000"), std::invalid_argument);
}

TEST(Decimals, parameters) {
    ASSERT_EQ(Fac::geometric("x", "0.5"), Fac::geometric("x", "1/2"));
    ASSERT_FALSE(Fac::geometric("x", "0.50000000000000000001") == Fac::geometric("x", "1/2"));
    ASSERT_EQ(Fac::uniform("x", "0", "3e0"), Fac::uniform("x", "0", "3"));
    ASSERT_EQ(Dist("p*x + (1-p)", "p") * "0.25", Dist("1/4*p*x + 1/4*(1-p)", "p"));
}

TEST(Decimals, inexactNumbers) {
    GiNaC::symbol x("x");
    ASSERT_TRUE(Dist::hasInexactNumbers(GiNaC::numeric(0.5) * x));
    ASSERT_FALSE(Dist::hasInexactNumbers(GiNaC::numeric(1, 2) * x + GiNaC::I));
}

TEST(Decimals, otherEntryPoints) {
    ASSERT_EQ(Dist::evaluate("2.5e1 * n", {{"n", 2}}), 50);
    ASSERT_EQ(Dist::find_symbols("0.5*x1 + 1e2*y"), (std::vector<std::string>{"x1", "y"}));
}
