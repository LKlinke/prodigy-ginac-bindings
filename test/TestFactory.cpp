//
// Created by Tobias Winkler on 22.10.21.
//

#include "gtest/gtest.h"
#include <DistributionFactory.h>

using namespace std;
using namespace prodigy;

TEST(DistributionFactory, geometric) {
    Dist::resetSymbolCache();
    EXPECT_NO_THROW(DistributionFactory::geometric("x", "1/2"));
    EXPECT_NO_THROW(DistributionFactory::geometric("x", "p"));

    EXPECT_ANY_THROW(DistributionFactory::geometric("z", "x")); // x was used as variable before
    EXPECT_ANY_THROW(DistributionFactory::geometric("y", "1-p")); // param must be numeric or single symbol
    EXPECT_ANY_THROW(DistributionFactory::geometric("p", "p")); // var and param cannot be the same
    EXPECT_ANY_THROW(DistributionFactory::geometric("123", "r")); // var must be a symbol
    EXPECT_ANY_THROW(DistributionFactory::geometric("2*c", "r")); // var must be a symbol
    EXPECT_THROW(DistributionFactory::geometric("#01-3==", "r"), GiNaC::parse_error);
}
