//
// Created by Tobias Winkler on 22.10.21.
//

#include <algorithm>
#include "Dist.h"
#include "macros.h"

using namespace GiNaC;

/* static (global) portion of Dist */

// initialize static variables of class Dist
Dist::str2sym Dist::_params;
Dist::str2sym Dist::_vars;

Dist::str2sym Dist::paramsAndVars() {
    str2sym res = _params;
    res.insert(_vars.begin(), _vars.end());
    return res;
}

bool Dist::isSymbKnown(const symbol &s) {
    return _params.count(s.get_name()) > 0 || _vars.count(s.get_name()) > 0;
}

bool Dist::isSymbKnown(const std::string &s) {
    return  _params.count(s) > 0 || _vars.count(s) > 0;
}

bool Dist::isKnownAsVar(const symbol &s) {
    return _vars.count(s.get_name()) > 0;
}

bool Dist::isKnownAsVar(const std::string &s) {
    return _vars.count(s) > 0;
}

bool Dist::isKnownAsParam(const symbol &s) {
    return _params.count(s.get_name()) > 0;
}

bool Dist::isKnownAsParam(const std::string  &s) {
    return _params.count(s) > 0;
}

// new symbols will be created for unknown names, but not registered
ex Dist::parseWithKnownSymbs(const std::string &s) {
    parser pars(paramsAndVars());
    ex parsed = pars(s); // might throw parse error
    return parsed;
}

void Dist::registerAllUnknownSymbsAsVars(const ex &e) {
    for(auto it = e.preorder_begin(); it != e.preorder_end(); ++it) {
        if(is_a<symbol>(*it)) {
            const auto& s = ex_to<symbol>(*it);
            if(!isSymbKnown(s)) {
                _vars[s.get_name()] = s;
            }
        }
    }
}

void Dist::registerAllUnknownSymbsAsParams(const ex &e) {
    for(auto it = e.preorder_begin(); it != e.preorder_end(); ++it) {
        if(is_a<symbol>(*it)) {
            const auto& s = ex_to<symbol>(*it);
            if(!isSymbKnown(s)) {
                _params[s.get_name()] = s;
            }
        }
    }
}

void Dist::registerAsVar(const symbol &p) {
    assert(!isKnownAsParam(p));
    _vars[p.get_name()] = p;
}

void Dist::registerAsParam(const symbol &p) {
    assert(!isKnownAsVar(p));
    _params[p.get_name()] = p;
}


/* public interface of dist */

Dist::Dist(const std::string &gfStr, const std::string &pStr) {
    ex p = parseWithKnownSymbs(pStr); // might throw parse error
    PRECOND(is_a<symbol>(p), "given parameter is not a symbol");
    const auto& pSym = ex_to<symbol>(p);
    PRECOND(!isKnownAsVar(pSym), "given parameter was registered as variable before");
    registerAsParam(pSym);
    ex gf = parseWithKnownSymbs(gfStr); // might throw parse error
    // todo if parseWithKnownSymbs(gfStr) throws then global state was modified without creating object successfully
    // todo issue warning if gf does not have p as subexpression
    registerAllUnknownSymbsAsVars(gf);
    _gf = gf;
}

Dist::Dist(const std::string &gfStr) {
    ex gf = parseWithKnownSymbs(gfStr); // might throw parse error
    registerAllUnknownSymbsAsVars(gf);
    _gf = gf;
}

ex Dist::mass() const {
    ex res {_gf};
    for(const auto& v : _vars) {
        res = res.subs(v.second == 1);
    }
    return res;
}

Dist Dist::marginal(const std::string& x) const {
    PRECOND(!isKnownAsParam(x), "given symbol is considered a parameter");
    ex res {_gf};
    for(const auto& v : _vars) {
        if(v.first != x) {
            res = res.subs(v.second == 1);
        }
    }
    return Dist {res};
}

ex Dist::E(const std::string &x) const {
    PRECOND(!isKnownAsParam(x), "given symbol is considered a parameter");
    if(isKnownAsVar(x)) {
        Dist marg = marginal(x);
        ex res {marg._gf};
        const auto& xSym = ex_to<symbol>(_vars.at(x));
        return res.diff(xSym).subs(xSym == 1);
    }
    else {
        // x is some string, neither a known param nor a known var
        parser pars;
        ex parsed = pars(x); // might throw parser error
        if(is_a<symbol>(parsed)) {
            // x is a symbol, interpret as trivial var
            return 0;
        }
        // expected value for expressions other than symbols not implemented
        throw prodigy::not_implemented_error();
    }
}

Dist Dist::update(const std::string &x, const std::string &e) const {
    PRECOND(!isKnownAsParam(x), "lhs of update may not be a parameter");
    ex xExp = parseWithKnownSymbs(x); // might throw parse error
    ex eExp = parseWithKnownSymbs(e);
    PRECOND(is_a<symbol>(xExp), "lhs of update must be a symbol");
    const auto& xSym = ex_to<symbol>(xExp);
    PRECOND(eExp.info(info_flags::polynomial), "only polynomial arithmetic is supported in updates");
    return _update(xSym, eExp); // might throw errors
}

Dist Dist::_update(const symbol &x, const ex &e) const {
    /*
     * algorithm for evaluating given assignment x = e = term_1 + ... + term_n :
     *
     * 1) tmp = 0
     * 2) tmp += term_1
     *    ...
     *    tmp += term_n
     * 3) x = tmp
     *
     */
    ex pol = e.expand();
    Dist res {*this};
    // 1) tmp = 0
    symbol tmp;
    // 2) tmp += term_1 ... tmp += term_n
    if(is_a<numeric>(pol) || is_a<symbol>(pol) || is_a<mul>(pol)) {
        // eg pol = 1 || pol = x || pol = 2*y
        res = res._incrAssgn(tmp, pol);
    } else {
        assert(is_a<add>(pol));
        // eg pol = 1 + x
        for(const auto& term : pol) {
            res = res._incrAssgn(tmp, term);
        }
    }
    // 3) x = tmp
    res._gf = res._gf.subs(x == 1);
    res._gf = res._gf.subs(tmp == x);
    // consider x a variable from now on
    registerAsVar(x);
    return res;
}

Dist Dist::filterLess(const std::string &x, const std::string &e) const {
    PRECOND(!isKnownAsParam(x), "given variable is parameter");
    parser parse;
    ex eExp = parse(e);
    PRECOND(eExp.info(info_flags::nonnegint), "given expression must be non-neg int");
    int order = ex_to<numeric>(eExp).to_int() ;
    PRECOND(isKnownAsVar(x), "given variable is not known");
    // todo one might accept unknown vars as well
    const auto& xExp = _vars.at(x);
    ex series = _gf.series(xExp, order);
    return Dist {series_to_poly(series)};
}

Dist Dist::_incrAssgn(const symbol &x, const ex &term) const {
    ex res {_gf};
    if(term.info(info_flags::nonnegint)) {
        // v += k, with k a non-negative int
        res *= pow(x, term);
    } else if(term.info(info_flags::symbol)) {
        // v += y, with y a variable
        res = res.subs(term == x * term);
    } else {
        throw prodigy::not_implemented_error{};
    }
    return Dist {res};
}


/* operators */

std::ostream& operator<<(std::ostream &os, const Dist &dist) {
    os << dist._gf;
    return os;
}

bool operator==(const Dist &lhs, const Dist &rhs) {
    return lhs._gf == rhs._gf;
}

Dist operator*(const Dist &lhs, const std::string &rhs) {
    ex parsed = Dist::parseWithKnownSymbs(rhs); // might throw parser error
    // todo error if parsed contained known variables
    // parsed contains just unknown symbs and known params
    // no error from now on
    Dist::registerAllUnknownSymbsAsParams(parsed);
    return Dist {parsed * lhs._gf};
}


