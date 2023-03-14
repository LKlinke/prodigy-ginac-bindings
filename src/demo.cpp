#include <iostream>
#include "DistributionFactory.h"

using namespace prodigy;

int main() {

    /* example 1 */


    Dist dist = DistributionFactory::geometric("x", "p") * Dist {"y^3"};
    std::cout << dist << std::endl;
    std::cout << "Total mass " << dist.mass() << std::endl;
    std::cout << "Expected value " << dist.E("x") << std::endl;
    std::cout << std::endl;
    
    std::cout << "Update " << dist.updateSubtraction("x", "x", "y") << std::endl;

    /* example 2 */

    Dist dist2 = DistributionFactory::geometric("y", "1/3");

    std::cout << dist2 << std::endl;
    std::cout << "Total mass " << dist2.mass() << std::endl;
    std::cout << "Expected value " << dist2.E("y") << std::endl;

    dist2 = Dist("p", "p") * dist2;
    std::cout << dist2 << std::endl;
    std::cout << "has parameter p " << Dist::isKnownAsParam("p") << std::endl;
    std::cout << "Total mass " << dist2.mass() << std::endl;
    std::cout << "Expected value " << dist2.E("y") << std::endl;


    return 0;
}
