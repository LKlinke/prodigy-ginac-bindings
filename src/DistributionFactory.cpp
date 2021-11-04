//
// Created by Tobias Winkler on 22.10.21.
//

#include <sstream>
#include "DistributionFactory.h"

using namespace GiNaC;

namespace prodigy {

    Dist DistributionFactory::geometric(const std::string &x, const std::string &p) {
        EXPECTS(!Dist::isKnownAsParam(x), "given variable is a parameter");
        parser parse;
        ex xExp = parse(x);
        EXPECTS(is_a<symbol>(xExp), "variable of geom distr must be a symbol");
        std::ostringstream ss;
        ss << "( 1 - (" << p << ") ) / (1 - (" << x << ") * (" << p << "))";
        ex pExp = parse(p);
        if (is_a<symbol>(pExp)) {
            return Dist{ss.str(), p};
        } else if (is_a<numeric>(pExp)) {
            return Dist{ss.str()};
        } else {
            EXPECTS(false, "parameter of geom distr must be a number or a symbol");
        }
    }

} // namespace prodigy