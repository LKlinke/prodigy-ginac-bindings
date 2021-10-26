//
// Created by Tobias Winkler on 22.10.21.
//

#include <DistributionFactory.h>
#include "gtest/gtest.h"

TEST(Dist, mass) {
    Dist dist = DistributionFactory::geometric(symbol("x"), numeric(1, 2));
    ASSERT_EQ(dist.mass(), 1);
}

TEST(Dist, E) {
    symbol p("p");
    Dist dist = DistributionFactory::geometric(symbol("x"), p);
    ASSERT_EQ(dist.E(dist.getVariableByName("x")), - p / (-1+p));
}

TEST(Dist, getVariables) {
    symbol x("x");
    symbol p("p");
    Dist dist = DistributionFactory::geometric(x, p);
    ASSERT_EQ(dist.getVariables(), std::vector<symbol>({x}));

    symbol xByName = dist.getVariableByName("x");
    ASSERT_EQ(xByName, x);
}

TEST(Dist, multiplication) {
    Dist dist = DistributionFactory::geometric(symbol("x"), numeric(1, 2));
    ASSERT_EQ((dist * numeric(2, 3)).mass(), numeric(2, 3));

    Dist distp = dist * std::string("p");
    ASSERT_TRUE(distp.hasParameter("p"));
    symbol p = distp.getParameterByName("p");
    EXPECT_ANY_THROW(distp.E(p)); // expected value in parameter should throw

    Dist distHalf = dist * std::string("1/2");
    ASSERT_FALSE(distHalf.hasParameter("p"));

    EXPECT_ANY_THROW(dist * std::string("1/2 + p"));
}

TEST(Dist, update) {
    symbol x("x");
    symbol y("y");
    Dist dist = DistributionFactory::geometric(x, numeric(1, 2));

    ASSERT_EQ(dist.update(x, x), dist);
    ASSERT_EQ(dist.update("x", "x"), dist);
    ASSERT_EQ(dist.update(x,1).gf(), x );
    ASSERT_EQ(dist.update("x","1").gf(), x );
    ASSERT_EQ(dist.update(x,x + 1).gf(), x * dist.gf() );
    ASSERT_EQ(dist.update("x","x + 1").gf(), x * dist.gf() );

    dist.update("y", "x");
}