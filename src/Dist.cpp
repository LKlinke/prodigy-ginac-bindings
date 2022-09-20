//
// Created by Tobias Winkler on 22.10.21.
//

#include <algorithm>
#include <numeric>
#include "Dist.h"
#include "macros.h"
#include <cassert>

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

    std::string Dist::getFreshVariable() {
        uint i = 0;
        do {
            std::stringstream s;
            s << "sym_" << i;
            std::string str = s.str();
            if(!isKnownAsParam(str) && !isKnownAsVar(str)) {
                registerAsVar(symbol(str));
                return str;
            }
            i++;
        } while (i > 0);

        throw std::overflow_error("cannot generate more fresh symbols");
    }


/* public interface of dist */

    Dist::Dist(const std::string &gfStr) {
        ex gf = parseWithKnownSymbs(gfStr); // might throw parse error
        registerAllUnknownSymbsAsVars(gf);
        _gf = gf;
    }

    Dist::Dist(const std::string &gfStr, const std::vector<std::string> &paramsStr) {
        for(const auto& pStr : paramsStr) {
            ex p = parseWithKnownSymbs(pStr); // might throw parse error
            EXPECTS(is_a<symbol>(p), "given parameter is not a symbol");
            const auto &pSym = ex_to<symbol>(p);
            EXPECTS(!isKnownAsVar(pSym), "given parameter was registered as variable before");
            // TODO we might want to not register the param if gf does not contain it at all
            registerAsParam(pSym);
        }
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

    int Dist::gfSize() const {
        int result = 0;
        for(auto it = _gf.preorder_begin(); it != _gf.preorder_end(); ++it) {
            result += 1;
        }
        return result;
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

    Dist Dist::marginal(const std::vector<std::string>& vars) const {
        for (const auto &v : vars) {
            EXPECTS(!isKnownAsParam(v), "given symbol is considered a parameter");
        }
        ex res{_gf};
        // go over all variables v we currently know
        for (const auto &v : _vars) {
            // if v is not in the given list vars then project it away
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
            return troolean::TRUE;
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
            return troolean::TRUE;
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
            return troolean::TRUE;
        } else {
            return troolean::UNKNOWN;
        }
    }

    troolean Dist::isPolynomial(const std::set<std::string>& vars) const {
        for (const auto &var : vars){
            EXPECTS(!isKnownAsParam(var), "given symbol is considered a parameter");
            // an expression is a polynomial in a list of variables iff it is a polynomial in all individual variables
            ex expr;
            if (isKnownAsVar(var)) {
                expr = _vars.at(var);
            } else {
                expr = symbol(var);
            }
            if(!_gf.is_polynomial(expr)) {
                return troolean::FALSE;
            }
        }
        return troolean::TRUE;
    }

    std::vector<Dist::State> Dist::getTerms(const std::set<std::string>& vars) const {
        EXPECTS(this->isPolynomial(vars) == troolean::TRUE, "can only iterate over terms of a finite polynomial");
        std::vector<ex> varEx;
        varEx.reserve(vars.size());
        for(const auto &var : vars) {
            EXPECTS(!isKnownAsParam(var), "given symbol is considered a parameter");
            if(isKnownAsVar(var)) {
                varEx.push_back(_vars.at(var));
            } else {
                varEx.push_back(symbol(var));
            }
        }

        auto handleProduct = [&](ex &term) {
            Dist::State state;
            for(ex var : varEx) {
                state.vals.push_back(std::pair<ex, int>(var, ex_to<numeric>(term.degree(var)).to_int()));
                term = term.subs(var == 1);
            }
            state.prob = term;
            return state;
        };

        ex poly = _gf.expand();
        auto nops = poly.nops();
        std::vector<Dist::State> res;
        res.reserve(nops);
        if(!is_a<add>(poly)) {
            res.push_back(handleProduct(poly));
            return res;
        }
        for(auto i = 0; i < nops; i++) {
            ex term = poly.op(i);
            res.push_back(handleProduct(term));
        }
        return res;
    }

    std::vector<Dist> Dist::arithmeticProgression(const std::string &var, const std::string &modulus) const {
        ex exMod = parseWithKnownSymbs(modulus);
        EXPECTS(is_a<numeric>(exMod), "given modulus is not a number");
        EXPECTS(ex_to<numeric>(exMod).is_pos_integer(), "given modulus is not a positive integer");
        ex exVar = parseWithKnownSymbs(var);
        EXPECTS(is_a<symbol>(exVar), "given variable is not a symbol");
        EXPECTS(!isKnownAsParam(var), "given variable is known as a parameter");

        ex root = exp(2 * Pi * I / exMod);
        std::vector<Dist> res;
        res.reserve(ex_to<numeric>(exMod).to_int());
        for(int rem = 0; rem < exMod; rem++) {
            ex psum = 0;
            for(int m = 0; m < exMod; m++) {
                psum += pow(root, -m * rem) * 
                    _gf.subs(exVar == pow(root, m) * exVar);
            }
            res.push_back(Dist(psum/exMod));
        }
        return res;
    }

    std::vector<symbol> Dist::getSymbols() const {
        std::vector<symbol> res;
        for(auto i = _gf.preorder_begin(); i != _gf.preorder_end(); i++) {
            if(is_a<symbol>(*i)) {
                res.push_back(ex_to<symbol>(*i));
            }
        }
        return res;
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
            size_t subexpr_count = term.nops();
            if (subexpr_count == 2){
                // todo this case probably has some errors and needs debugging!!
                // something like x*y || 2*x
                ex first=term.op(0);
                ex second=term.op(1);
                if (first.info(info_flags::symbol) && second.info(info_flags::nonnegint)){
                    // e.g. y*2
                    res = res.subs(term == x*pow(first, second));
                } else if (second.info(info_flags::symbol) && first.info(info_flags::nonnegint)){
                    // e.g. 2*y
                    res = res.subs(term == x * pow(second, first));
                } else if (first.info(info_flags::symbol) && second.info(info_flags::negint)){
                    // e.g. -y || -1 * y
                    //std::cerr << "Warning! This might lead to undesired behavior!" << std::endl;
                    res = res.subs(first == first*pow(x, second));
                } else if (second.info(info_flags::symbol) && first.info(info_flags::negint)){
                    // e.g. y * -1
                    //std::cerr << "Warning! This might lead to undesired behavior!" << std::endl;
                    res = res.subs(term == x*pow(second, first));
                } else{
                    EXPECTS(false, "multiplication of two variables not supported");
                }
            }
            else throw prodigy::not_implemented_error{};
            // todo implement multiplication for finite-support case
        } else {
            EXPECTS(false, "requested update operation is not supported");
        }
        return Dist{res};
    }

    Dist Dist::updateVar(const std::string &updatedVar, const std::string &assignVar) const {
        EXPECTS(!isKnownAsParam(updatedVar), "lhs is known as a parameter");
        EXPECTS(is_a<symbol>(parseWithKnownSymbs(updatedVar)), "lhs is not a symbol");
        registerAllUnknownSymbsAsVars(symbol(updatedVar));
        if(updatedVar == assignVar) {
            return *this;
        }

        EXPECTS(!isKnownAsParam(assignVar), "given rhs is a parameter");
        ex rhs = parseWithKnownSymbs(assignVar);
        ex lhs = _vars.at(updatedVar);
        if(is_a<symbol>(rhs)) {
            return Dist(_gf.subs(lhs == 1).subs(rhs == lhs * rhs));
        } else {
            return Dist(_gf.subs(lhs == 1) * pow(lhs, rhs));
        }
    }

    Dist Dist::updateSum(const std::string &updatedVar, const std::string &firstSummand, const std::string secondSummand) const {
        EXPECTS(!isKnownAsParam(updatedVar), "lhs is known as a parameter");
        EXPECTS(is_a<symbol>(parseWithKnownSymbs(updatedVar)), "lhs is not a symbol");
        registerAllUnknownSymbsAsVars(symbol(updatedVar));
        ex lhs = parseWithKnownSymbs(updatedVar);
        ex sum1 = parseWithKnownSymbs(firstSummand);
        ex sum2 = parseWithKnownSymbs(secondSummand);

        // we add two variables
        if(isKnownAsVar(firstSummand) && isKnownAsVar(secondSummand)) {
            if(sum2 == lhs) {
                ex tmp = sum1;
                sum1 = sum2;
                sum2 = tmp;
            }
            if(sum1 == lhs) {
                if(sum2 == lhs) {
                    return Dist(_gf.subs(lhs == pow(lhs, 2)));
                } else {
                    return Dist(_gf.subs(sum2 == sum2*lhs));
                }
            } else {
                return Dist(_gf.subs(lhs == 1).subs(sum1 == sum1 * lhs).subs(sum2 == sum2 * lhs));
            }
        }
        // we add a variable and a literal / parameter
        else if(isKnownAsVar(firstSummand) || isKnownAsVar(secondSummand)) {
            ex var = sum1;
            ex lit = sum2;
            if(isKnownAsVar(secondSummand)) {
                var = sum2;
                lit = sum1;
            }
            ex res = _gf;
            if(var != lhs) {
                res = res.subs(lhs == 1).subs(var == lhs * var);
            }
            return Dist(res * pow(lhs, lit));
        }
        // we add two literals / parameters
        else {
            return Dist(_gf.subs(lhs == 1) * pow(lhs, sum1 + sum2));
        }
    }

    Dist Dist::updateIid(const std::string &x, const Dist &d, const std::string &y) const {
        EXPECTS(!isKnownAsParam(x) && !isKnownAsParam(y), "given symbol is considered a parameter, not a variable");
        EXPECTS(d.mass() == 1, "the distribution wrapped into iid must not be a sub distribution");
        EXPECTS(d.isUnivariate() == troolean::TRUE, "the distribution wrapped into iid must be univariate");
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
        // todo why not use parseWithKnownSymbols?
        if(!isKnownAsVar(x)) {
            xExp = symbol {x};
            Dist::registerAsVar(ex_to<symbol>(xExp));
        } else {
            xExp = _vars.at(x);
        }
        ex series = _gf.series(xExp, order);
        return Dist{series_to_poly(series)};
    }

    Dist Dist::CoefficientIterator::next() {
        ex nextTerm = series_to_poly(_rest.series(_var == 0, _n));
        _n++;
        _rest = (_rest - nextTerm).normal();
        ex nextCoeff = nextTerm.subs(_var == 1);
        return Dist(nextCoeff);
    }

    Dist::CoefficientIterator Dist::coefficientIterator(const std::string &var) const {
        EXPECTS(isKnownAsVar(var), "can only expand series in known variables");
        return CoefficientIterator(*this, ex_to<symbol>(_vars.at(var)));
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
        return (lhs - rhs).isZero() == troolean::TRUE;

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
