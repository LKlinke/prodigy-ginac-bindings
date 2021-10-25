//
// Created by Tobias Winkler on 22.10.21.
//

#include <DistributionFactory.h>
#include "gtest/gtest.h"

TEST(ParametricDist, mass) {
    ParametricDist dist = DistributionFactory::geometric(symbol("x"), numeric(1,2));
    ASSERT_EQ(dist.mass(), 1);
}

TEST(ParametricDist, E) {
    symbol p("p");
    ParametricDist dist = DistributionFactory::geometric(symbol("x"), p);
    ASSERT_EQ(dist.E(dist.getVariableByName("x")), - p / (-1+p));
}

TEST(ParametricDist, getVariables) {
    symbol x("x");
    symbol p("p");
    ParametricDist dist = DistributionFactory::geometric(x, p);
    ASSERT_EQ(dist.getVariables(), std::vector<symbol>({x}));

    symbol xByName = dist.getVariableByName("x");
    ASSERT_EQ(xByName, x);
}

TEST(ParametricDist, multiplication) {
    ParametricDist dist = DistributionFactory::geometric(symbol("x"), numeric(1, 2));
    ASSERT_EQ((dist * numeric(2, 3)).mass(), numeric(2, 3));

    ParametricDist distp = dist * std::string("p");
    ASSERT_TRUE(distp.hasParameter("p"));
    symbol p = distp.getParameterByName("p");
    EXPECT_ANY_THROW(distp.E(p)); // expected value in parameter should throw

    ParametricDist distHalf = dist * std::string("1/2");
    ASSERT_FALSE(distHalf.hasParameter("p"));

    EXPECT_ANY_THROW(dist * std::string("1/2 + p"));
}