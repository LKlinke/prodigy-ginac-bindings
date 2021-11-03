//
// Created by Tobias Winkler on 22.10.21.
//

#ifndef PRODIGY_DISTRIBUTIONFACTORY_H
#define PRODIGY_DISTRIBUTIONFACTORY_H

#include "Dist.h"

namespace prodigy {

    class DistributionFactory {

    public:
        static Dist geometric(const std::string& x, const std::string& p);
    };

} // namespace prodigy




#endif //PRODIGY_DISTRIBUTIONFACTORY_H
