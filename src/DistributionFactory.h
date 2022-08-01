//
// Created by Tobias Winkler on 22.10.21.
//

#ifndef PRODIGY_DISTRIBUTIONFACTORY_H
#define PRODIGY_DISTRIBUTIONFACTORY_H

#include "Dist.h"

namespace prodigy {

    class DistributionFactory {

    public:
        /*
         * constructs a distribution where x is p-geometrically distributed
         */
        static Dist geometric(const std::string& x, const std::string& p);

        /*
         * constructs a distribution where x is uniformly distributed over [a,b]
         */
        static Dist uniform(const std::string& x, const std::string& a, const std::string& b);

        /*
         * constructs a distribution where x is n with probability 1
         */
        static Dist dirac(const std::string& x, const std::string& n);

        /*
         * constructs a Catalan distribution with parameter p
         */
        static Dist catalan(const std::string& x, const std::string& p);
    };

} // namespace prodigy

#endif //PRODIGY_DISTRIBUTIONFACTORY_H
