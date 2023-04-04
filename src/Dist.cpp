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

    bool Dist::isSymbKnown(const possymbol &s) {
        return _params.count(s.get_name()) > 0 || _vars.count(s.get_name()) > 0;
    }

    bool Dist::isSymbKnown(const std::string &s) {
        return _params.count(s) > 0 || _vars.count(s) > 0;
    }

    bool Dist::isKnownAsVar(const possymbol &s) {
        return _vars.count(s.get_name()) > 0;
    }

    bool Dist::isKnownAsVar(const std::string &s) {
        return _vars.count(s) > 0;
    }

    bool Dist::isKnownAsParam(const possymbol &s) {
        return _params.count(s.get_name()) > 0;
    }

    bool Dist::isKnownAsParam(const std::string &s) {
        return _params.count(s) > 0;
    }

    ex Dist::parseWithKnownSymbs(const std::string &s) {
        auto syms = paramsAndVars();
        for(auto found : find_symbols(s)) {
            if(!isKnownAsVar(found) && !isKnownAsParam(found)) {
                syms[found] = possymbol(found);
            }
        }
        parser pars(syms);
        ex parsed = pars(s); // might throw parse error
        return parsed;
    }

    void Dist::registerAllUnknownSymbsAsVars(const ex &e) {
        for (auto it = e.preorder_begin(); it != e.preorder_end(); ++it) {
            if (is_a<symbol>(*it)) {
                EXPECTS(is_a<possymbol>(*it), "variable is not a possymbol");
                const auto &s = ex_to<possymbol>(*it);
                if (!isSymbKnown(s)) {
                    _vars[s.get_name()] = s;
                }
            }
        }
    }

    void Dist::registerAllUnknownSymbsAsParams(const ex &e) {
        for (auto it = e.preorder_begin(); it != e.preorder_end(); ++it) {
            if (is_a<symbol>(*it)) {
                EXPECTS(is_a<possymbol>(*it), "parameter is not a possymbol");
                const auto &s = ex_to<possymbol>(*it);
                if (!isSymbKnown(s)) {
                    _params[s.get_name()] = s;
                }
            }
        }
    }

    void Dist::registerAsVar(const possymbol &p) {
        assert(!isKnownAsParam(p));
        _vars[p.get_name()] = p;
    }

    void Dist::registerAsParam(const possymbol &p) {
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
                registerAsVar(possymbol(str));
                return str;
            }
            i++;
        } while (i > 0);

        throw std::overflow_error("cannot generate more fresh symbols");
    }

    int Dist::evaluate(const std::string &expr, const std::map<std::string, int> &valuations) {
        parser p;
        ex res = p(expr);
        for(auto sym : p.get_syms()) {
            res = res.subs(sym.second == valuations.at(sym.first));
        }
        EXPECTS(res.info(info_flags::integer), "result of valuation is not an integer");
        return ex_to<numeric>(res).to_int();
    }

    std::vector<std::string> Dist::find_symbols(const std::string &s) {
        // this function does not use the symbol registry, as opposed to Dist::getSymbols()
        // to use the modulo operator %, we need to replace it with something ginac can parse
        std::string str = s;
        std::replace(str.begin(), str.end(), '%', '/');

        GiNaC::parser p;
        ex parsed = p(str);
        std::vector<std::string> res;
        for(auto i : p.get_syms()) {
            res.push_back(i.first);
        }
        return res;
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
            EXPECTS(is_a<possymbol>(p), "given parameter is not a symbol");
            const auto &pSym = ex_to<possymbol>(p);
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

    ex Dist::mass(const std::optional<std::set<std::string>> vars) const {
        ex res{_gf};
        /* TODO
         * normal() cancels common factors, this is needed since we are going to substitute stuff to avoid div by 0
         * encapsulate this somehow because substitution occurs elsewhere too
         */
        res = res.normal();
        if(!vars.has_value()) {
            for (const auto &v : _vars) {
                res = res.subs(v.second == 1);
            }
        } else {
            for (auto v : vars.value()) {
                res = res.subs(_vars.at(v) == 1);
            }
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
        ex res{_gf.expand().normal()};
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
        ex res{_gf.normal()};
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
                EXPECTS(is_a<possymbol>(*it), "symbol is not a possymbol");
                const auto& s = ex_to<possymbol>(*it);
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
        lst varEx;
        for (const auto &var : vars){
            EXPECTS(!isKnownAsParam(var), "given symbol is considered a parameter");
            if (isKnownAsVar(var)) {
                varEx.append(_vars.at(var));
            } else {
                varEx.append(possymbol(var));
            }
        }
        if(_gf.normal().is_polynomial(varEx)) {
            return troolean::TRUE;
        }
        return troolean::FALSE;
    }

    troolean Dist::isFiniteSupport() const {
        lst vars;
        for(auto var : this->getSymbols()) {
            if(isKnownAsVar(var)) {
                vars.append(var);
            }
        }
        if(_gf.is_polynomial(vars)) {
            return troolean::TRUE;
        }
        return troolean::FALSE;
    }

    std::vector<Dist::State> Dist::getTerms(const std::set<std::string>& vars) const {
        lst varEx;
        for(const auto &var : vars) {
            EXPECTS(!isKnownAsParam(var), "given symbol is considered a parameter");
            if(isKnownAsVar(var)) {
                varEx.append(_vars.at(var));
            } else {
                varEx.append(possymbol(var));
            }
        }
        return getTerms(varEx);
    }

    std::vector<Dist::State> Dist::getTerms() const {
        lst vars;
        for(auto var : this->getSymbols()) {
            if(isKnownAsVar(var)) {
                vars.append(var);
            }
        }
        return getTerms(vars);
    }

    std::vector<Dist::State> Dist::getTerms(const lst &vars) const {
        auto normal = _gf.normal();
        EXPECTS(normal.is_polynomial(vars), "can only iterate over terms of a finite polynomial");
        auto handleProduct = [&](ex &term) {
            Dist::State state;
            for(ex var : vars) {
                state.vals[ex_to<symbol>(var).get_name()] = ex_to<numeric>(term.degree(var)).to_int();
                term = term.subs(var == 1);
            }
            state.prob = term;
            return state;
        };

        ex poly = normal.expand().collect(vars, true);
        auto nops = poly.nops();
        std::vector<Dist::State> res;
        res.reserve(nops);
        if(!is_a<add>(poly)) {
            res.push_back(handleProduct(poly));
        } else {
            for(auto i = 0; i < nops; i++) {
                ex term = poly.op(i);
                res.push_back(handleProduct(term));
            }
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
                psum += GiNaC::pow(root, -m * rem) * 
                    _gf.subs(exVar == GiNaC::pow(root, m) * exVar);
            }
            res.push_back(Dist(psum/exMod));
        }
        return res;
    }

    std::vector<possymbol> Dist::getSymbols() const {
        std::vector<possymbol> res;
        for(auto i = _gf.preorder_begin(); i != _gf.preorder_end(); i++) {
            if(is_a<symbol>(*i)) {
                EXPECTS(is_a<possymbol>(*i), "symbol is not a possymbol");
                if(res.empty() || std::find(res.begin(), res.end(), *i) == res.end()) {
                    res.push_back(ex_to<possymbol>(*i));
                }
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
            return res.diff(xSym).normal().subs(xSym == 1);
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

    Dist Dist::updateVar(const std::string &updatedVar, const std::string &assignVar) const {
        EXPECTS(!isKnownAsParam(updatedVar), "lhs is known as a parameter");
        EXPECTS(is_a<possymbol>(parseWithKnownSymbs(updatedVar)), "lhs is not a possymbol");
        registerAllUnknownSymbsAsVars(possymbol(updatedVar));
        if(updatedVar == assignVar) {
            return *this;
        }

        EXPECTS(!isKnownAsParam(assignVar), "given rhs is a parameter");
        ex rhs = parseWithKnownSymbs(assignVar);
        ex lhs = _vars.at(updatedVar);
        if(is_a<symbol>(rhs)) {
            EXPECTS(is_a<possymbol>(rhs), "rhs is not a possymbol");
            return Dist(_gf.subs(lhs == 1).subs(rhs == lhs * rhs));
        } else {
            return Dist((rhs != 0 ? _gf : _gf.normal()).subs(lhs == 1) * GiNaC::pow(lhs, rhs));
        }
    }

    Dist Dist::updateSum(const std::string &updatedVar, const std::string &firstSummand, const std::string secondSummand) const {
        EXPECTS(!isKnownAsParam(updatedVar), "lhs is known as a parameter");
        EXPECTS(is_a<symbol>(parseWithKnownSymbs(updatedVar)), "lhs is not a symbol");
        registerAllUnknownSymbsAsVars(possymbol(updatedVar));
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
                    return Dist(_gf.subs(lhs == GiNaC::pow(lhs, 2)));
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
            return Dist(res * GiNaC::pow(lhs, lit));
        }
        // we add two literals / parameters
        else {
            return Dist(_gf.subs(lhs == 1) * GiNaC::pow(lhs, sum1 + sum2));
        }
    }

    Dist Dist::updateProduct(const std::string &updatedVar, const std::string &firstFactor, const std::string &secondFactor, 
            const std::set<std::string> &variables, const bool &finite, const std::optional<std::string> &approx) const {
        EXPECTS(!isKnownAsParam(updatedVar), "lhs is known as a parameter");
        EXPECTS(is_a<symbol>(parseWithKnownSymbs(updatedVar)), "lhs is not a symbol");
        EXPECTS(!isKnownAsParam(firstFactor) && !isKnownAsParam(secondFactor), "cannot perform multiplication of parameters");
        registerAllUnknownSymbsAsVars(possymbol(updatedVar));
        ex lhs = parseWithKnownSymbs(updatedVar);
        ex fac1 = parseWithKnownSymbs(firstFactor);
        ex fac2 = parseWithKnownSymbs(secondFactor);

        // we multiply two variables
        if(isKnownAsVar(firstFactor) && isKnownAsVar(secondFactor)) {
            // infinite distribution
            if(!finite) {
                auto marginal_l = this->marginal(firstFactor);
                auto finite_left = marginal_l.isPolynomial({firstFactor}) == troolean::TRUE;
                auto marginal_r = this->marginal(secondFactor);
                auto finite_right = marginal_r.isPolynomial({secondFactor}) == troolean::TRUE;

                if(!finite_left && !finite_right) {
                    if(!approx.has_value()) {
                        throw std::invalid_argument("cannot perform multiplication of two variables with infinite range");
                    }
                    marginal_l = marginal_l.approximateUnilaterally(firstFactor, approx.value());
                    finite_left = true;
                }

                Dist finite;
                std::string finite_var, infinite_var;
                bool other_finite;
                if(finite_left) {
                    finite = marginal_l;
                    finite_var = firstFactor;
                    infinite_var = secondFactor;
                    other_finite = finite_right;
                } else {
                    finite = marginal_r;
                    finite_var = secondFactor;
                    infinite_var = firstFactor;
                    other_finite = finite_left; // == false
                }
                ex res = 0;
                for(auto state : finite.getTerms({finite_var})) {
                    // state.vals contains only one item
                    auto val = std::to_string(state.vals[finite_var]);
                    res += this->filterEq(finite_var, val)
                        .updateProduct(updatedVar, val, infinite_var, variables, other_finite, approx)._gf;
                }
                return Dist(res);
            }
            // finite distribution
            else {
                ex res = _gf;
                for(auto state : this->getTerms(variables)) {
                    ex term = parseWithKnownSymbs(state.toString());
                    res -= term;
                    term = term.subs(lhs == 1) * GiNaC::pow(lhs, state.vals.at(firstFactor) * state.vals.at(secondFactor));
                    res += term;
                }
                return Dist(res);
            }
        }
        // we multiply a variable with a literal
        else if(isKnownAsVar(firstFactor) || isKnownAsVar(secondFactor)) {
            ex var, lit;
            if(isKnownAsVar(firstFactor)) {
                var = fac1;
                lit = fac2;
            } else {
                var = fac2;
                lit = fac1;
            }
            if(var == lhs) {
                return Dist(_gf.subs(lhs == GiNaC::pow(lhs, lit)));
            } else {
                return Dist(_gf.subs(lhs == 1).subs(var == var * GiNaC::pow(lhs, lit)));
            }
        }
        // we multiply two literals
        else {
            return Dist(_gf.subs(lhs == 1) * GiNaC::pow(lhs, fac1 * fac2));
        }
    }

    Dist Dist::updateSubtraction(const std::string updatedVar, const std::string subFrom, const std::string sub) const {
        EXPECTS(!isKnownAsParam(updatedVar), "lhs is known as a parameter");
        EXPECTS(is_a<symbol>(parseWithKnownSymbs(updatedVar)), "lhs is not a symbol");
        EXPECTS(!isKnownAsParam(subFrom), "cannot subtract from a parameter");
        EXPECTS(!isKnownAsParam(sub), "cannot subtract a parameter");
        registerAllUnknownSymbsAsVars(possymbol(updatedVar));
        ex lhs = parseWithKnownSymbs(updatedVar);
        ex sub1 = parseWithKnownSymbs(subFrom);
        ex sub2 = parseWithKnownSymbs(sub);

        ex res = _gf;
        // we subtract a variable from another variable
        if(isKnownAsVar(sub) && isKnownAsVar(subFrom)) {
            if(this->isFiniteSupport() == troolean::TRUE) {
                for(auto state: this->getTerms()) {
                    EXPECTS(state.vals[subFrom] >= state.vals[sub], "result can be negative");
                }
            }
            else {
                // TODO how to check for this case? We probably need limits
                // ex testFun = dist.marginal(updatedVar)._gf.normal().subs(lhs == 0);
            }

            if(sub2 == lhs) {
                if(sub1 == lhs) {
                    res = _gf.subs(lhs == 1);
                } else {
                    res = _gf.subs(lhs == GiNaC::pow(lhs, -1)).subs(sub1 == sub1 * lhs);
                }
            } else {
                if(sub1 != lhs) {
                    res = res.subs(lhs == 1).subs(sub1 == sub1 * lhs);
                }
                res = res.subs(sub2 == sub2 * GiNaC::pow(lhs, -1));
            }
        }
        //we subtract a literal from a variable
        else if(isKnownAsVar(subFrom)) {
            EXPECTS(this->filterLess(subFrom, sub) == Dist("0"), "result can be negative");
            if(lhs != sub1) {
                res = res.subs(lhs == 1).subs(sub1 == sub1 * lhs);
            }
            res = res * GiNaC::pow(lhs, -sub2);
        }
        // we subtract a variable from a literal
        else if(isKnownAsVar(sub)) {
            EXPECTS(this->filterGreater(sub, subFrom) == Dist("0"), "result can be negative");
            if(sub2 == lhs) {
                res = _gf.subs(lhs == GiNaC::pow(lhs, -1)) * GiNaC::pow(lhs, sub1);
            } else {
                res = (_gf.subs(lhs == 1) * GiNaC::pow(lhs, sub1)).subs(sub2 == sub2 * GiNaC::pow(lhs, -1));
            }
        }
        // we subtract two literals from each other
        else {
            EXPECTS(sub1 >= sub2, "result can be negative");
            res = _gf.subs(lhs == 1) * GiNaC::pow(lhs, sub1 - sub2);
        }

        return Dist(res);
            
    }

    Dist Dist::updateModulo(const std::string &updatedVar, const std::string &left, const std::string &right, 
            const std::set<std::string> &variables, const std::optional<bool> &finite, const std::optional<std::string> &approx) const {
        EXPECTS(!isKnownAsParam(updatedVar), "lhs is known as a parameter");
        EXPECTS(is_a<symbol>(parseWithKnownSymbs(updatedVar)), "lhs is not a symbol");
        EXPECTS(!isKnownAsParam(left) && !isKnownAsParam(right), "cannot perform modulo of parameters");
        registerAllUnknownSymbsAsVars(possymbol(updatedVar));
        ex lhs = parseWithKnownSymbs(updatedVar);

        // this is not the same as finite.value_or, which will always evaluate isFiniteSupport()
        bool is_finite = finite.has_value() ? finite.value() : this->isFiniteSupport() == troolean::TRUE;

        // on finite distributions, iterate over all states
        if(is_finite) {
            ex res = 0;
            for(auto state : this->getTerms(variables)) {
                numeric leftVar, rightVar;
                if(isKnownAsVar(left)) {
                    leftVar = state.vals.at(left);
                } else {
                    leftVar = ex_to<numeric>(parseWithKnownSymbs(left));
                }
                if(isKnownAsVar(right)) {
                    rightVar = state.vals.at(right);
                } else {
                    rightVar = ex_to<numeric>(parseWithKnownSymbs(right));
                }
                res += parseWithKnownSymbs(state.toString()).subs(lhs == 1) * GiNaC::pow(lhs, mod(leftVar, rightVar));
            }
            return Dist(res);
        }
        // if the distribution is infinite and right is a variable, it needs to have finite range
        else if(isKnownAsVar(right)) {
            auto marginal_r = this->marginal(right);
            if(marginal_r.isPolynomial({right}) == troolean::FALSE) {
                if(!approx.has_value()) {
                    throw std::invalid_argument("rhs needs to have a finite marginal");
                }
                marginal_r = marginal_r.approximateUnilaterally(right, approx.value());
            }
            ex res = 0;
            for(auto state : marginal_r.getTerms({right})) {
                auto right_val = std::to_string(state.vals.at(right));
                res += this->filterEq(right, right_val)
                    .updateModulo(updatedVar, left, right_val, variables, is_finite ? true : std::optional<bool>(), approx)._gf;
            }
            return Dist(res);
        }
        // if left is a variable, it doesn't have to have finite range
        else if(isKnownAsVar(left)) {
            auto marginal_l = this->marginal(left);
            ex res = 0;
            if(marginal_l.isPolynomial({left}) == troolean::TRUE) {
                for(auto state : marginal_l.getTerms({left})) {
                    auto left_val = std::to_string(state.vals.at(left));
                    res += this->filterEq(left, left_val)
                        .updateModulo(updatedVar, left_val, right, variables, is_finite ? true : std::optional<bool>(), approx)._gf;
                }
            } else {
                auto prog = this->arithmeticProgression(left, right);
                for(int i = 0; i < prog.size(); i++) {
                    res += prog.at(i)._gf.subs(lhs == 1) * GiNaC::pow(lhs, i);
                }
            }
            return Dist(res);
        }
        // if both are not variables, simply compute the result
        else {
            numeric mod1 = ex_to<numeric>(parseWithKnownSymbs(left));
            numeric mod2 = ex_to<numeric>(parseWithKnownSymbs(right));
            return this->updateVar(updatedVar, std::to_string(GiNaC::to_int(mod(mod1, mod2))));
        }
    }

    Dist Dist::updateDivision(const std::string &updatedVar, const std::string &numerator, const std::string &denominator, const std::optional<std::string> &approx) const {
        EXPECTS(!isKnownAsParam(updatedVar), "lhs is known as a parameter");
        EXPECTS(is_a<symbol>(parseWithKnownSymbs(updatedVar)), "lhs is not a symbol");
        EXPECTS(!isKnownAsParam(numerator) && !isKnownAsParam(denominator), "cannot perform division of parameters");
        registerAllUnknownSymbsAsVars(possymbol(updatedVar));
        ex lhs = parseWithKnownSymbs(updatedVar);
        ex div1 = parseWithKnownSymbs(numerator);
        ex div2 = parseWithKnownSymbs(denominator);

        Dist marginal_l, marginal_r;
        if(isKnownAsVar(numerator)) {
            marginal_l = this->marginal(numerator);
            if(marginal_l.isPolynomial({numerator}) == troolean::FALSE) {
                if(!approx.has_value()) {
                    throw std::invalid_argument("numerator must have a finite marginal");
                }
                marginal_l = marginal_l.approximateUnilaterally(numerator, approx.value());
            }
        } else {
            marginal_l = Dist(div1);
        }
        if(isKnownAsVar(denominator)) {
            marginal_r = this->marginal(denominator);
            if(marginal_r.isPolynomial({denominator}) == troolean::FALSE) {
                if(!approx.has_value()) {
                    throw std::invalid_argument("denominator must have a finite marginal");
                }
                marginal_r = marginal_r.approximateUnilaterally(denominator, approx.value());
            }
        } else {
            marginal_r = Dist(div2);
        }

        ex res = 0;
        if(isKnownAsVar(numerator)) {
            for(auto state_l : marginal_l.getTerms({numerator})) {
                if(isKnownAsVar(denominator)) {
                    for(auto state_r : marginal_r.getTerms({denominator})) {
                        int val_l = state_l.vals.at(numerator);
                        int val_r = state_r.vals.at(denominator);
                        ex x = this->filterEq(numerator, std::to_string(val_l))
                            .filterEq(denominator, std::to_string(val_r))._gf;
                        EXPECTS(mod(val_l, val_r) == 0 || x == 0, "numerator is not always divisible by denominator");
                        res += x.subs(lhs == 1) * GiNaC::pow(lhs, val_l / val_r);
                    }
                } else {
                    int val_l = state_l.vals.at(numerator);
                    ex x = this->filterEq(numerator, std::to_string(val_l))._gf;
                    EXPECTS(mod(val_l, ex_to<numeric>(div2)) == 0 || x == 0, "numerator is not always divisible by denominator");
                    res += x.subs(lhs == 1) * GiNaC::pow(lhs, val_l / div2);
                }
            }
        } else {
            if(isKnownAsVar(denominator)) {
                for(auto state_r : marginal_r.getTerms({denominator})) {
                    int val_r = state_r.vals.at(denominator);
                    ex x = this->filterEq(denominator, std::to_string(val_r))._gf;
                    EXPECTS(mod(ex_to<numeric>(div1), val_r) == 0 || x == 0, "numerator is not always divisible by denominator");
                    res += x.subs(lhs == 1) * GiNaC::pow(lhs, div1 / val_r);
                } 
            } else {
                EXPECTS(mod(ex_to<numeric>(div1), ex_to<numeric>(div2)) == 0, "numerator is not always divisible by denominator");
                res = _gf.subs(lhs == 1) * GiNaC::pow(lhs, div1 / div2);
            }
        }

        return Dist(res);
    }

    Dist Dist::updatePower(const std::string &updatedVar, const std::string &base, const std::string &exp, const std::optional<std::string> &approx) const {
        EXPECTS(!isKnownAsParam(updatedVar), "lhs is known as a parameter");
        EXPECTS(is_a<symbol>(parseWithKnownSymbs(updatedVar)), "lhs is not a symbol");
        EXPECTS(!isKnownAsParam(base) && !isKnownAsParam(exp), "cannot perform exponentiation of parameters");
        registerAllUnknownSymbsAsVars(possymbol(updatedVar));
        ex lhs = parseWithKnownSymbs(updatedVar);
        ex pow1 = parseWithKnownSymbs(base);
        ex pow2 = parseWithKnownSymbs(exp);
        ex res = this->_gf;

        // variable to the power of a variable
        if(isKnownAsVar(base) && isKnownAsVar(exp)) {
            auto marginal_l = this->marginal(base);
            auto marginal_r = this->marginal(exp);
            if(marginal_l.isPolynomial({base}) == troolean::FALSE) {
                if(!approx.has_value()) {
                    throw std::invalid_argument("base does not have a finite marginal");
                }
                marginal_l = marginal_l.approximateUnilaterally(base, approx.value());
            }
            if(marginal_l.isPolynomial({exp}) == troolean::FALSE) {
                if(!approx.has_value()) {
                    throw std::invalid_argument("exponent does not have a finite marginal");
                }
                marginal_r = marginal_r.approximateUnilaterally(exp, approx.value());
            }

            for(auto state_l : marginal_l.getTerms({base})) {
                for(auto state_r : marginal_r.getTerms({exp})) {
                    int b = state_l.vals.at(base);
                    int e = state_r.vals.at(exp);
                    ex x = this->filterEq(base, std::to_string(b)).filterEq(exp, std::to_string(e))._gf;
                    res -= x;
                    res += x.subs(lhs == 1) * GiNaC::pow(lhs, GiNaC::pow(b, e));
                }
            }
        }
        // variable to the power of a literal
        else if(isKnownAsVar(base)) {
            auto marginal = this->marginal(base);
            EXPECTS(marginal.isPolynomial({base}) == troolean::TRUE, 
                    "Can only perform exponentiation if the base has a finite marginal");

            for(auto state : marginal.getTerms({base})) {
                int b = state.vals.at(base);
                ex x = this->filterEq(base, std::to_string(b))._gf;
                res -= x;
                res += x.subs(lhs == 1) * GiNaC::pow(lhs, GiNaC::pow(b, pow2));
            }
        }
        // literal to the power of a variable
        else if(isKnownAsVar(exp)) {
            auto marginal = this->marginal(exp);
            EXPECTS(marginal.isPolynomial({exp}) == troolean::TRUE, 
                    "Can only perform exponentiation if the base has a finite marginal");

            for(auto state : marginal.getTerms({exp})) {
                int e = state.vals.at(exp);
                ex x = this->filterEq(exp, std::to_string(e))._gf;
                res -= x;
                res += x.subs(lhs == 1) * GiNaC::pow(lhs, GiNaC::pow(pow1, e));
            }
        }
        // literal to the power of a literal
        else {
            res = res.subs(lhs == 1) * GiNaC::pow(lhs, GiNaC::pow(pow1, pow2));
        }

        return Dist(res);
    }

    Dist Dist::approximateUnilaterally(const std::string &variable, const std::string &prob_mass) const {
        ex mass = parseWithKnownSymbs(prob_mass).evalf();
        EXPECTS(is_a<numeric>(mass), "given probability mass is not a number");
        if(mass == 0) {
            return Dist("0");
        }
        EXPECTS(mass < this->mass(), "given probability mass is too large");
        EXPECTS(mass > 0, "given probability mass must be non-negative");
        EXPECTS(isKnownAsVar(variable), "given variable is not known as variable");
        ex res = 0;
        ex mass_res = 0;

        int i = 0;
        while(mass_res.evalf() < mass) {
            auto next = this->filterEq(variable, std::to_string(i));
            i++;
            res += next._gf;
            mass_res += next.mass();
        }
        return Dist(res);
    }

    Dist Dist::updateIid(const std::string &x, const Dist &d, const std::string &y) const {
        EXPECTS(!isKnownAsParam(x) && !isKnownAsParam(y), "given symbol is considered a parameter, not a variable");
        EXPECTS(d.mass() == 1, "the distribution wrapped into iid must not be a sub distribution");
        EXPECTS(d.isUnivariate() == troolean::TRUE, "the distribution wrapped into iid must be univariate");
        if (!isKnownAsVar(y)) {
            // then this is like x = iid(y, 0) so x simply becomes 0
            registerAsVar(possymbol(y)); // consider z variable from now on
            return updateVar(x, "0");
        }
        if (!isKnownAsVar(x)) {
            registerAsVar(possymbol(x));
        }
        // first marginalize x away in *this
        ex res = _gf.normal();
        if (x != y) {
            res = res.subs(_vars[x] == 1);
        }
        // ySubs will be gf of y with (single) var replaced by x
        // do this by simply substituting all known vars for x
        // currently, gathering the symbols that actually occur in y is not implemented
        ex ySubs = d._gf;
        for (const auto &v : _vars) {
            ySubs = ySubs.subs(v.second == _vars[x]);
        }
        // do the actual substitution
        if (x != y) {
            res = res.subs(_vars[y] == _vars[y] * ySubs);
        } else {
            res = res.subs(_vars[y] == ySubs);
        }
        return Dist{res};
    }

    Dist Dist::filterLess(const std::string &x, const std::string &e) const {
        EXPECTS(!isKnownAsParam(x), "given symbol is a parameter");
        parser parse;
        ex eExp = parse(e);
        EXPECTS(eExp.info(info_flags::nonnegint), "given expression must be non-neg int");
        int order = ex_to<numeric>(eExp).to_int();
        if(order == 0) {
            // nothing can be less than 0
            return Dist("0");
        }
        ex xExp;
        // todo why not use parseWithKnownSymbols?
        if(!isKnownAsVar(x)) {
            xExp = possymbol {x};
            Dist::registerAsVar(ex_to<possymbol>(xExp));
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
        return CoefficientIterator(*this, ex_to<possymbol>(_vars.at(var)));
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
