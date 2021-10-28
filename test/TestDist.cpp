//
// Created by Tobias Winkler on 22.10.21.
//

#include <DistributionFactory.h>
#include "gtest/gtest.h"

using namespace std;
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
    dist = dist.update("x", "x + 5");
    ASSERT_EQ(dist.E("x"), 6);
}

TEST(Dist, filter) {
    Dist::resetSymbolCache();
    Dist dist = Fac::geometric("x", "1/2");
    ASSERT_EQ(dist.filterLess("x", "3").mass(), GiNaC::numeric(7,8));
    ASSERT_EQ(dist.filterLeq("x", "2").mass(), GiNaC::numeric(7,8));
}