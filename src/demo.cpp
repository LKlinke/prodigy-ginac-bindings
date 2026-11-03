#include <iostream>
#include "DistributionFactory.h"
#include <ginac/ginac.h>
using namespace prodigy;

int main() {

    /* example 1 */
    std::cout << "Example 1" << std::endl;
    Dist dist = DistributionFactory::geometric("x", "p") * Dist {"y^3"};
    std::cout << dist << std::endl;
    std::cout << "Total mass " << dist.mass() << std::endl;
    std::cout << "Expected value " << dist.E("x") << std::endl;
    std::cout << std::endl;
    
    std::cout << "Update " << dist.updateSubtraction("x", "x", "y") << std::endl;

    /* example 2 */
    std::cout << "Example 2" << std::endl;
    Dist dist2 = DistributionFactory::geometric("y", "1/3");

    std::cout << dist2 << std::endl;
    std::cout << "Total mass " << dist2.mass() << std::endl;
    std::cout << "Expected value " << dist2.E("y") << std::endl;

    dist2 = Dist("p", "p") * dist2;
    std::cout << dist2 << std::endl;
    std::cout << "has parameter p " << dist2.isKnownAsParam("p") << std::endl;
    std::cout << "Total mass " << dist2.mass() << std::endl;
    std::cout << "Expected value " << dist2.E("y") << std::endl;


    /* Example 3*/
    std::cout << "Example 3" << std::endl;
    Dist dist3 = DistributionFactory::uniform("x", "a", "b");
    Dist dist4 = DistributionFactory::uniform("x", "a", "b");

    std::cout << (dist3 == dist4) << std::endl;


    return 0;
}
