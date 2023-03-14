//
// Created by Tobias Winkler on 22.10.21.
//

#include <DistributionFactory.h>
#include "gtest/gtest.h"

using namespace std;
using namespace prodigy;
using Fac = DistributionFactory;

TEST(Dist, mass) {
    Dist::resetSymbolCache();
    Dist dist = Fac::geometric("x", "1/2");
    ASSERT_EQ(dist.mass(), 1);
    dist = Fac::geometric("x", "p");
    ASSERT_EQ(dist.mass(), 1);
}

TEST(Dist, marginal) {
    Dist::resetSymbolCache();
    Dist dist = Fac::geometric("x", "1/2") * Fac::geometric("y", "1/2");
    dist = dist.marginal("x");
    ASSERT_EQ(dist, Fac::geometric("x", "1/2"));
}

TEST(Dist, E) {
    Dist::resetSymbolCache();
    Dist dist = Fac::geometric("x", "1/2");
    ASSERT_EQ(dist.E("x"), 1);
}

TEST(Dist, update) {
    Dist::resetSymbolCache();
    Dist dist = Fac::geometric("x", "1/2");
    dist = dist.updateSum("x", "x", "5");
    ASSERT_EQ(dist.E("x"), 6);
}

TEST(Dist, updateMonus) {
    Dist::resetSymbolCache();
    Dist dist {"1/2 + 1/2 * x^2"};
    ASSERT_EQ(dist.updateSubtraction("x", "x", "1"), Dist {"1/2 + 1/2 * x"});
}

TEST(Dist, iid) {
    Dist::resetSymbolCache();
    Dist  dist = Dist {"1/2*x + 1/2*x^2"} * Dist {"y^2"};
    cout << dist.updateIid("res", Fac::geometric("T", "1/3"), "y").marginal("res") << endl;
}

TEST(Dist, filter) {
    Dist::resetSymbolCache();
    Dist dist = Fac::geometric("x", "1/2");
    ASSERT_EQ(dist.filterLess("x", "3").mass(), GiNaC::numeric(7,8));
    ASSERT_EQ(dist.filterLeq("x", "2").mass(), GiNaC::numeric(7,8));
}

TEST(Dist, filterUnknownVar) {
    Dist::resetSymbolCache();
    Dist dist {"1"};
    ASSERT_TRUE(dist.isTrivial() == troolean::TRUE);
    ASSERT_EQ(dist.filterLess("x", "1").mass(), GiNaC::numeric(1));
}

TEST(Dist, gfSize) {
    Dist::resetSymbolCache();
    ASSERT_EQ(Dist("(x+1)^2").gfSize(), Dist("(x+1)^100").gfSize());
    ASSERT_EQ(Dist("x").gfSize(), 1);
    ASSERT_EQ(Dist("x+y").gfSize(), 3);
}

TEST(Dist, coefficientIterator) {
    Dist::resetSymbolCache();
    Dist dist = Fac::geometric("x", "q");
    Dist::CoefficientIterator it = dist.coefficientIterator("x");
    ASSERT_EQ(it.next() + it.next() + it.next(), Dist("1-q - q^2 + q + q^2 - q^3", "q"));
}