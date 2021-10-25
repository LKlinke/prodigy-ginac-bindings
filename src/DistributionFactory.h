//
// Created by Tobias Winkler on 22.10.21.
//

#ifndef PRODIGY_DISTRIBUTIONFACTORY_H
#define PRODIGY_DISTRIBUTIONFACTORY_H

#include "ParametricDist.h"

class DistributionFactory {

public:
    static ParametricDist geometric(const symbol& x, const ex& p);
};




#endif //PRODIGY_DISTRIBUTIONFACTORY_H
