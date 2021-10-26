#include <iostream>
#include "Formula.h"
#include "DistributionFactory.h"

int main() {

    /* example 1 */

    symbol x("x");
    symbol p("p");

    Dist dist = DistributionFactory::geometric(x, p);

    std::cout << dist << std::endl;
    std::cout << "Total mass " << dist.mass() << std::endl;
    std::cout << "Expected value " << dist.E(x) << std::endl;
    std::cout << std::endl;

    /* example 2 */

    symbol y("y");
    Dist dist2 = DistributionFactory::geometric(y, numeric("1/3"));

    std::cout << dist2 << std::endl;
    std::cout << "Total mass " << dist2.mass() << std::endl;
    std::cout << "Expected value " << dist2.E(y) << std::endl;

    dist2 = std::string("p") * dist2;
    std::cout << dist2 << std::endl;
    std::cout << "has parameter p " << dist.hasParameter("p") << std::endl;
    std::cout << "Total mass " << dist2.mass() << std::endl;
    std::cout << "Expected value " << dist2.E(y) << std::endl;


    return 0;
}
