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
        EXPECTS(is_a<possymbol>(xExp), "variable of unif distr must be a symbol, got " + x);
        EXPECTS(is_a<possymbol>(aExp) || aExp.info(info_flags::nonnegint), "parameters of unif distr must be a non-neg int or a symbol, got " + a);
        EXPECTS(is_a<possymbol>(bExp) || bExp.info(info_flags::nonnegint), "parameters of unif distr must be a non-neg int or a symbol, got " + b);
        EXPECTS(!(aExp.info(info_flags::integer) && aExp.info(info_flags::integer)) || ex_to<numeric>(aExp).to_int() <= ex_to<numeric>(bExp).to_int(), "interval [" + a + "," + b + "] is invalid");
        // register variable/parameters in case we don't know them yet
        if(is_a<symbol>(aExp) && !Dist::isKnownAsParam(a)) {
            EXPECTS(is_a<possymbol>(aExp), "not a possymbol");
            Dist::registerAsParam(ex_to<possymbol>(aExp));
        }
        if(is_a<symbol>(bExp) && !Dist::isKnownAsParam(b)) {
            EXPECTS(is_a<possymbol>(aExp), "not a possymbol");
            Dist::registerAsParam(ex_to<possymbol>(bExp));
        }
        if(!Dist::isKnownAsVar(x)) Dist::registerAsVar(ex_to<possymbol>(xExp));
        ex gf = pow(xExp, aExp) / (bExp - aExp + 1) * (1 - pow(xExp, bExp - aExp + 1)) / (1 - xExp);
        return Dist(gf);
    }

    Dist DistributionFactory::dirac(const std::string &x, const std::string &n) {
        EXPECTS(!Dist::isKnownAsParam(x), "given variable " + x + " was previously declared a parameter");
        EXPECTS(!Dist::isKnownAsVar(n), "given parameter " + n + " was previously declared a variable");
        ex xExp = Dist::parseWithKnownSymbs(x);
        ex nExp = Dist::parseWithKnownSymbs(n);
        EXPECTS(is_a<possymbol>(xExp), "variable of Dirac distr must be a symbol, got " + x);
        EXPECTS(is_a<symbol>(nExp) || nExp.info(info_flags::nonnegint), "parameter of Dirac distribution must be a non-neg int or a symbol, got " + n);
        // register variable/parameters in case we don't know them yet
        if(is_a<symbol>(nExp) && !Dist::isKnownAsParam(n)) {
            EXPECTS(is_a<possymbol>(nExp), "not a possymbol");
            Dist::registerAsParam(ex_to<possymbol>(nExp));
        }
        if(!Dist::isKnownAsVar(x)) Dist::registerAsVar(ex_to<possymbol>(xExp));
        ex gf = pow(xExp, nExp);
        return Dist(gf);
    }

    Dist DistributionFactory::catalan(const std::string &x, const std::string &p) {
        EXPECTS(!Dist::isKnownAsParam(x), "given variable " + x + " was previously declared a parameter");
        EXPECTS(!Dist::isKnownAsVar(p), "given parameter " + p + " was previously declared a variable");
        ex xExp = Dist::parseWithKnownSymbs(x);
        ex pExp = Dist::parseWithKnownSymbs(p);
        EXPECTS(is_a<possymbol>(xExp), "variable of Catalan distr must be a symbol, got " + x);
        EXPECTS(is_a<symbol>(pExp) || pExp.info(info_flags::positive), "parameter of Catalan distribution must be a positive real or a symbol, got " + p);
        // register variable/parameters in case we don't know them yet
        if(is_a<symbol>(pExp) && !Dist::isKnownAsParam(p)) {
            EXPECTS(is_a<possymbol>(pExp), "not a possymbol");
            Dist::registerAsParam(ex_to<possymbol>(pExp));
        }
        if(!Dist::isKnownAsVar(x)) Dist::registerAsVar(ex_to<possymbol>(xExp));
        ex gf = (1 - sqrt(1 - 4 * pExp * (1-pExp) * xExp)) / (2 * pExp);
        return Dist(gf);
    }

} // namespace prodigy