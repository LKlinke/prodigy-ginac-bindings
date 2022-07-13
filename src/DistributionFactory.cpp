//
// Created by Tobias Winkler on 22.10.21.
//

#include <sstream>
#include "DistributionFactory.h"

using namespace GiNaC;

namespace prodigy {

    Dist DistributionFactory::geometric(const std::string &x, const std::string &p) {
        EXPECTS(!Dist::isKnownAsParam(x), "given variable " + x + " is a parameter");
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

    Dist DistributionFactory::uniform(const std::string& x, const std::string& a, const std::string& b) {
        EXPECTS(!Dist::isKnownAsParam(x), "given variable " + x + " was previously declared a parameter");
        EXPECTS(!Dist::isKnownAsVar(a), "given parameter " + a + " was previously declared a variable");
        EXPECTS(!Dist::isKnownAsVar(b), "given parameter " + b + " was previously declared a variable");
        ex xExp = Dist::parseWithKnownSymbs(x);
        ex aExp = Dist::parseWithKnownSymbs(a);
        ex bExp = Dist::parseWithKnownSymbs(b);
        EXPECTS(is_a<symbol>(xExp), "variable of unif distr must be a symbol, got " + x);
        EXPECTS(is_a<symbol>(aExp) || aExp.info(info_flags::nonnegint), "parameters of unif distr must be a non-neg int or a symbol, got " + a);
        EXPECTS(is_a<symbol>(bExp) || bExp.info(info_flags::nonnegint), "parameters of unif distr must be a non-neg int or a symbol, got " + b);
        EXPECTS(!(aExp.info(info_flags::integer) && aExp.info(info_flags::integer)) || ex_to<numeric>(aExp).to_int() <= ex_to<numeric>(bExp).to_int(), "interval [" + a + "," + b + "] is invalid");
        // register variable/parameters in case we don't know them yet
        if(is_a<symbol>(aExp) && !Dist::isKnownAsParam(a)) {
            Dist::registerAsParam(ex_to<symbol>(aExp));
        }
        if(is_a<symbol>(bExp) && !Dist::isKnownAsParam(b)) {
            Dist::registerAsParam(ex_to<symbol>(bExp));
        }
        if(!Dist::isKnownAsVar(x)) Dist::registerAsVar(ex_to<symbol>(xExp));
        ex gf = pow(xExp, aExp) / (bExp - aExp + 1) * (1 - pow(xExp, bExp - aExp + 1)) / (1 - xExp);
        return Dist(gf);
    }

} // namespace prodigy