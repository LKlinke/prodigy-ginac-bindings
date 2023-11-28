//
// Created by Tobias Winkler on 22.10.21.
//

#include <sstream>
#include "DistributionFactory.h"

using namespace GiNaC;

namespace prodigy {

    Dist DistributionFactory::geometric(const std::string &x, const std::string &p) {
        Dist res;

        parser parse;
        ex xExp = parse(x);
        EXPECTS(is_a<symbol>(xExp), "variable of geom distr must be a symbol");
        res.registerAsVar(ex_to<symbol>(xExp));
        ex pExp = parse(p);
        if (is_a<symbol>(pExp)) {
            res.registerAsParam(ex_to<possymbol>(pExp));
        }
        else EXPECTS(pExp.info(info_flags::nonnegative), "parameter " + p + " is either a non-negative number or not a symbol.");
        ex gf = (pExp) / (1 - (1 - (pExp)) * xExp);
        res._gf = gf;
        return res;
    }

    Dist DistributionFactory::uniform(const std::string& x, const std::string& a, const std::string& b) {
        parser pars;
        ex xExp = pars(x);
        ex aExp = pars(a);
        ex bExp = pars(b);
        EXPECTS(is_a<symbol>(xExp), "variable of unif distr must be a symbol, got " + x);
        EXPECTS(is_a<symbol>(aExp) || aExp.info(info_flags::nonnegint), "parameters of unif distr must be a non-neg int or a symbol, got " + a);
        EXPECTS(is_a<symbol>(bExp) || bExp.info(info_flags::nonnegint), "parameters of unif distr must be a non-neg int or a symbol, got " + b);
        EXPECTS(!(aExp.info(info_flags::integer) && aExp.info(info_flags::integer)) || ex_to<numeric>(aExp).to_int() <= ex_to<numeric>(bExp).to_int(), "interval [" + a + "," + b + "] is invalid");

        Dist res;
        xExp = possymbol(x);
        if(is_a<symbol>(aExp) && is_a<symbol>(bExp)){
            aExp = possymbol(a);
            bExp = possymbol(b);
            res._vars[x] = xExp;
            res._params[a] = aExp;
            res._params[b] = bExp;
        } else if (is_a<symbol>(aExp)) {
            aExp = possymbol(a);
            res._vars[x] = xExp;
            res._params[a] = aExp;
        } else if (is_a<symbol>(bExp)) {
            bExp = possymbol(b);
            res._vars[x] = xExp;
            res._params[b] = bExp;
        }
        ex gf = pow(xExp, aExp) / (bExp - aExp + 1) * (1 - pow(xExp, bExp - aExp + 1)) / (1 - xExp);
        res._gf = gf;
        return res;
    }

    Dist DistributionFactory::dirac(const std::string &x, const std::string &n) {
        parser pars;
        ex xExp = pars(x);
        ex nExp = pars(n);
        EXPECTS(is_a<symbol>(xExp), "variable of Dirac distr must be a symbol, got " + x);
        EXPECTS(is_a<symbol>(nExp) || nExp.info(info_flags::nonnegint), "parameter of Dirac distribution must be a non-neg int or a symbol, got " + n);
        // register variable/parameters in case we don't know them yet
        Dist result;
        if(is_a<symbol>(nExp)) {
            result.registerAsVar(ex_to<symbol>(xExp));
            result.registerAsParam(ex_to<symbol>(nExp));
        }
        ex gf = pow(xExp, nExp);
        result._gf = gf;
        return result;
    }

    Dist DistributionFactory::catalan(const std::string &x, const std::string &p) {
        parser pars;
        ex xExp = pars(x);
        ex pExp = pars(p);
        EXPECTS(is_a<symbol>(xExp), "variable of Catalan distr must be a symbol, got " + x);
        EXPECTS(is_a<symbol>(pExp) || pExp.info(info_flags::positive), "parameter of Catalan distribution must be a positive real or a symbol, got " + p);
        
        Dist res;
        if(is_a<symbol>(pExp)) {
            res.registerAsVar(ex_to<symbol>(xExp));
            res.registerAsParam(ex_to<possymbol>(pExp));
        }
        ex gf = (1 - sqrt(1 - 4 * pExp * (1-pExp) * xExp)) / (2 * pExp);
        res._gf = gf;
        return res;
    }

} // namespace prodigy