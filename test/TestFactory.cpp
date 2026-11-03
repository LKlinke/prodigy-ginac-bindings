//
// Created by Tobias Winkler on 22.10.21.
//

#include "gtest/gtest.h"
#include <DistributionFactory.h>

using namespace std;
using namespace prodigy;

TEST(DistributionFactory, geometric) {
    EXPECT_NO_THROW(DistributionFactory::geometric("x", "1/2"));
    EXPECT_NO_THROW(DistributionFactory::geometric("x", "p"));

    EXPECT_ANY_THROW(DistributionFactory::geometric("y", "1-p")); // param must be numeric or single symbol
    EXPECT_ANY_THROW(DistributionFactory::geometric("p", "p")); // var and param cannot be the same
    EXPECT_ANY_THROW(DistributionFactory::geometric("123", "r")); // var must be a symbol
    EXPECT_ANY_THROW(DistributionFactory::geometric("2*c", "r")); // var must be a symbol
    EXPECT_THROW(DistributionFactory::geometric("#01-3==", "r"), GiNaC::parse_error);
}

TEST(DistributionFactory, uniform) {

    EXPECT_NO_THROW(DistributionFactory::uniform("x", "0", "1"));
    EXPECT_NO_THROW(DistributionFactory::uniform("x", "0", "0"));
    EXPECT_NO_THROW(DistributionFactory::uniform("x", "a", "b"));
    EXPECT_NO_THROW(DistributionFactory::uniform("x", "0", "b"));
    EXPECT_NO_THROW(DistributionFactory::uniform("x", "a", "10"));

    EXPECT_ANY_THROW(DistributionFactory::uniform("x", "x", "2")); // x was a variable before
    EXPECT_ANY_THROW(DistributionFactory::uniform("2*x", "1", "2")); // var must be a symbol
    EXPECT_ANY_THROW(DistributionFactory::uniform("x", "3", "2")); // invalid interval
    EXPECT_ANY_THROW(DistributionFactory::uniform("x", "a + b", "2")); // expression not allowed a param

    Dist unif = DistributionFactory::uniform("x", "1", "3");
    EXPECT_EQ(unif.mass(), 1);
    EXPECT_EQ(unif.filterEq("x", "1").mass(), GiNaC::numeric(1,3));
    EXPECT_EQ(unif.filterEq("x", "2").mass(), GiNaC::numeric(1,3));
    EXPECT_EQ(unif.filterEq("x", "3").mass(), GiNaC::numeric(1,3));
    EXPECT_EQ(unif.filterEq("x", "4").mass(), 0);
}
