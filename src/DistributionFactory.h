//
// Created by Tobias Winkler on 22.10.21.
//

#ifndef PRODIGY_DISTRIBUTIONFACTORY_H
#define PRODIGY_DISTRIBUTIONFACTORY_H

#include "Dist.h"

class DistributionFactory {

public:
    static Dist geometric(const symbol& x, const ex& p);
};




#endif //PRODIGY_DISTRIBUTIONFACTORY_H
