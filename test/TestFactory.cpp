//
// Created by Tobias Winkler on 22.10.21.
//

#include "gtest/gtest.h"
#include <DistributionFactory.h>

TEST(DistributionFactory, geometric) {
    EXPECT_NO_THROW(DistributionFactory::geometric(symbol("x"), symbol("p")));
    EXPECT_NO_THROW(DistributionFactory::geometric(symbol("x"), numeric(1,2)));

    // EXPECT_ANY_THROW(DistributionFactory::geometric(symbol("x"), symbol("x"))); // todo
    EXPECT_ANY_THROW(DistributionFactory::geometric(symbol("x"), 0.5)); // cannot use float as param
    EXPECT_ANY_THROW(DistributionFactory::geometric(symbol("x"), numeric(3,2))); // param must be in correct interval
}
