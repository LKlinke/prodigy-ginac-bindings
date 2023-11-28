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

TEST(Dist, iid) {
    Dist  dist = Dist {"1/2*x + 1/2*x^2"} * Dist {"y^2"};
    dist = dist.updateIid("res", Fac::geometric("T","1/3"), "y");
    ASSERT_EQ(dist, Dist {"(1/2 * x + 1/2 * x^2) * y^2"} * Fac::geometric("res", "1/3") * Fac::geometric("res", "1/3"));
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