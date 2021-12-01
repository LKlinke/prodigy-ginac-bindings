//
// Created by Tobias Winkler on 22.10.21.
//

#include <algorithm>
#include <numeric>
#include "Dist.h"
#include "macros.h"

using namespace GiNaC;

namespace prodigy {

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
        return _params.count(s) > 0 || _vars.count(s) > 0;
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

    bool Dist::isKnownAsParam(const std::string &s) {
        return _params.count(s) > 0;
    }

// new symbols will be created for unknown names, but not registered
    ex Dist::parseWithKnownSymbs(const std::string &s) {
        parser pars(paramsAndVars());
        ex parsed = pars(s); // might throw parse error
        return parsed;
    }

    void Dist::registerAllUnknownSymbsAsVars(const ex &e) {
        for (auto it = e.preorder_begin(); it != e.preorder_end(); ++it) {
            if (is_a<symbol>(*it)) {
                const auto &s = ex_to<symbol>(*it);
                if (!isSymbKnown(s)) {
                    _vars[s.get_name()] = s;
                }
            }
        }
    }

    void Dist::registerAllUnknownSymbsAsParams(const ex &e) {
        for (auto it = e.preorder_begin(); it != e.preorder_end(); ++it) {
            if (is_a<symbol>(*it)) {
                const auto &s = ex_to<symbol>(*it);
                if (!isSymbKnown(s)) {
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

    void Dist::deregister(const std::string &s) {
        if (isKnownAsParam(s)) {
            _params.erase(s);
        }
        if (isKnownAsVar(s)) {
            _vars.erase(s);
        }
    }


/* public interface of dist */

    Dist::Dist(const std::string &gfStr, const std::string &pStr) {
        ex p = parseWithKnownSymbs(pStr); // might throw parse error
        EXPECTS(is_a<symbol>(p), "given parameter is not a symbol");
        const auto &pSym = ex_to<symbol>(p);
        EXPECTS(!isKnownAsVar(pSym), "given parameter was registered as variable before");
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

    const ex &Dist::gf() const {
        return _gf;
    }

    ex Dist::mass() const {
        ex res{_gf};
        /* TODO
         * normal() cancels common factors, this is needed since we are going to substitute stuff to avoid div by 0
         * encapsulate this somehow because substitution occurs elsewhere too
         */
        res = res.normal();
        for (const auto &v : _vars) {
            res = res.subs(v.second == 1);
        }
        return res;
    }

    Dist Dist::normalize() const {
        return Dist{_gf / mass()};
    }

    Dist Dist::marginal(const std::string &x) const {
        EXPECTS(!isKnownAsParam(x), "given symbol is considered a parameter");
        ex res{_gf};
        for (const auto &v : _vars) {
            if (v.first != x) {
                res = res.subs(v.second == 1);
            }
        }
        return Dist{res};
    }

    Dist Dist::marginal(std::initializer_list<std::string> vars) const {
        for (const auto &v : vars) {
            EXPECTS(!isKnownAsParam(v), "given symbol is considered a parameter");
        }
        ex res{_gf};
        for (const auto &v : _vars) {
            if (std::find(vars.begin(), vars.end(), v.first) == vars.end()) {
                res = res.subs(v.second == 1);
            }
        }
        return Dist{res};
    }

    troolean Dist::isZero() const {
        // todo look deeper into what is the best way in ginac for testing whether expression is zero
        bool ginacThinksIsZero = _gf.normal().expand().is_zero();
        return makeTroolean(ginacThinksIsZero);
    }

    troolean Dist::isTrivial() const {
        ex gfNormal = _gf.normal();
        if(gfNormal.info(info_flags::numeric) && 0 <= ex_to<numeric>(gfNormal) && ex_to<numeric>(gfNormal) <= 1) {
            return troolean::YES;
        }
        return troolean::UNKNOWN;
    }

    troolean Dist::isUnivariate(std::string& res) const {
        ex gfNormal = _gf.normal();
        str2sym found;
        for(auto it = gfNormal.preorder_begin(); it != gfNormal.preorder_end(); ++it) {
            if(is_a<symbol>(*it)) {
                const auto& s = ex_to<symbol>(*it);
                if(isKnownAsParam(s)) continue;
                assert(isKnownAsVar(s));
                found[s.get_name()] = s;
                res = s.get_name();
            }
        }
        if(found.size() > 1) {
            return troolean::UNKNOWN;
        } else {
            return troolean::YES;
        }
    }

    troolean Dist::isUnivariate() const {
        std::string trash;
        return isUnivariate(trash);
    }

    troolean Dist::areIndependent(const std::string &x, const std::string &y) const {
        EXPECTS(!isKnownAsParam(x) && !isKnownAsParam(y), "given symbol is considered a parameter");
        // todo think about normalization
        Dist mx = marginal(x).normalize();
        Dist my = marginal(y).normalize();
        Dist mxy = marginal({x, y}).normalize();
        if (mx * my == mxy) {
            return troolean::YES;
        } else {
            return troolean::UNKNOWN;
        }
    }


    ex Dist::E(const std::string &x) const {
        EXPECTS(!isKnownAsParam(x), "given symbol is considered a parameter");
        if (isKnownAsVar(x)) {
            Dist marg = marginal(x);
            ex res{marg._gf};
            const auto &xSym = ex_to<symbol>(_vars.at(x));
            return res.diff(xSym).subs(xSym == 1);
        } else {
            // x is some string, neither a known param nor a known var
            parser pars;
            ex parsed = pars(x); // might throw parser error
            if (is_a<symbol>(parsed)) {
                // x is a symbol, interpret as trivial var
                return 0;
            }
            // expected value for expressions other than symbols not implemented
            throw prodigy::not_implemented_error();
        }
    }

    Dist Dist::update(const std::string &x, const std::string &e) const {
        EXPECTS(!isKnownAsParam(x), "lhs of update may not be a parameter");
        ex xExp = parseWithKnownSymbs(x); // might throw parse error
        ex eExp = parseWithKnownSymbs(e);
        EXPECTS(is_a<symbol>(xExp), "lhs of update must be a symbol");
        const auto &xSym = ex_to<symbol>(xExp);
        EXPECTS(eExp.info(info_flags::polynomial), "only polynomial arithmetic is supported in updates");
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
        Dist res{*this};
        // 1) tmp = 0
        symbol tmp;
        // 2) tmp += term_1 ... tmp += term_n
        if (is_a<numeric>(pol) || is_a<symbol>(pol) || is_a<mul>(pol)) {
            // eg pol = 1 || pol = x || pol = 2*y
            res = res._incrAssgn(tmp, pol);
        } else {
            assert(is_a<add>(pol));
            // eg pol = 1 + x
            for (const auto &term : pol) {
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

    Dist Dist::_incrAssgn(const symbol &x, const ex &term) const {
        ex res{_gf};
        if (term.info(info_flags::nonnegint)) {
            // x += k, with k a non-negative int
            res *= pow(x, term);
        } else if (term.info(info_flags::negint)) {
            // x += k, with k a neg int
            // we implement ``monus''
            int k = ex_to<numeric>(term).to_int();
            ex smallerK = series_to_poly(res.series(x, -k));
            ex rest = res - smallerK;
            res = rest*pow(x, term) + smallerK.subs(x == 1);
        } else if (term.info(info_flags::symbol)) {
            // x += y, with y a variable
            res = res.subs(term == x * term);
        } else if (is_a<mul>(term)) {
            throw prodigy::not_implemented_error{};
            // todo implement multiplication for finite-support case
        } else {
            EXPECTS(false, "requested update operation is not supported");
        }
        return Dist{res};
    }

    Dist Dist::updateIid(const std::string &x, const Dist &d, const std::string &y) const {
        EXPECTS(!isKnownAsParam(x) && !isKnownAsParam(y), "given symbol is considered a parameter, not a variable");
        EXPECTS(d.mass() == 1, "the distribution wrapped into iid must not be a sub distribution");
        EXPECTS(d.isUnivariate() == troolean::YES, "the distribution wrapped into iid must be univariate");
        if (!isKnownAsVar(y)) {
            // then this is like x = iid(y, 0) so x simply becomes 0
            registerAsVar(symbol(y)); // consider z variable from now on
            return update(x, "0");
        }
        if (!isKnownAsVar(x)) {
            registerAsVar(symbol(x));
        }
        // first marginalize x away in *this
        ex res = _gf.normal().subs(_vars[x] == 1);
        // ySubs will be gf of y with (single) var replaced by x
        // do this by simply substituting all known vars for x
        // currently, gathering the symbols that actually occur in y is not implemented
        ex ySubs = d._gf;
        for (const auto &v : _vars) {
            ySubs = ySubs.subs(v.second == _vars[x]);
        }
        // do the actual substitution
        res = res.subs(_vars[y] == _vars[y] * ySubs);
        return Dist{res};
    }

    Dist Dist::filterLess(const std::string &x, const std::string &e) const {
        EXPECTS(!isKnownAsParam(x), "given symbol is a parameter");
        parser parse;
        ex eExp = parse(e);
        EXPECTS(eExp.info(info_flags::nonnegint), "given expression must be non-neg int");
        int order = ex_to<numeric>(eExp).to_int();
        ex xExp;
        if(!isKnownAsVar(x)) {
            xExp = symbol {x};
            Dist::registerAsVar(ex_to<symbol>(xExp));
        } else {
            xExp = _vars.at(x);
        }
        ex series = _gf.series(xExp, order);
        return Dist{series_to_poly(series)};
    }

    Dist Dist::substituteParam(const std::string &p, const std::string &val) const {
        EXPECTS(!isKnownAsVar(p), "given symbol is considered a variable");
        if (!isKnownAsParam(p)) {
            return *this;
        }
        ex pSym = _params[p];
        ex parsed = Dist::parseWithKnownSymbs(val);
        Dist::registerAllUnknownSymbsAsParams(parsed);
        // todo error if parsed contains known variables
        ex res = _gf.subs(pSym == parsed);
        return Dist{res};
    }


/* operators */

    std::ostream &operator<<(std::ostream &os, const Dist &dist) {
        os << dist._gf.normal();
        return os;
    }

    bool operator==(const Dist &lhs, const Dist &rhs) {
        return (lhs - rhs).isZero() == troolean::YES;

    }

    Dist operator*(const Dist &lhs, const std::string &rhs) {
        ex parsed = Dist::parseWithKnownSymbs(rhs); // might throw parser error
        // todo error if parsed contained known variables
        // parsed contains just unknown symbs and known params
        // no error from now on
        Dist::registerAllUnknownSymbsAsParams(parsed);
        return Dist{parsed * lhs._gf};
    }

} // namespace prodigy