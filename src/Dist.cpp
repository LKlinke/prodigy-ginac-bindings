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

    std::vector<std::string> Dist::find_symbols(const std::string &s) {
        // this function does not use the symbol registry, as opposed to Dist::getSymbols()
        // to use the modulo operator %, we need to replace it with something GiNaC can parse
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

    Dist::str2sym Dist::paramsAndVars() const {
        str2sym res = _params;
        res.insert(_vars.begin(), _vars.end());
        return res;
    }

    const bool Dist::isSymbKnown(const symbol &s) const {
        return _params.find(s.get_name()) != _params.end() || _vars.find(s.get_name()) != _vars.end();
    }

    const bool Dist::isSymbKnown(const std::string &s) const {
        return _params.find(s) != _params.end() || _vars.find(s) != _vars.end();
    }

    const bool Dist::isKnownAsVar(const symbol &s) const {
        return _vars.find(s.get_name()) != _vars.end();
    }

    const bool Dist::isKnownAsVar(const std::string &s) const {
        return _vars.find(s) != _vars.end();
    }

    const bool Dist::isKnownAsParam(const symbol &s) const {
        return _params.find(s.get_name()) != _params.end();
    }

    const bool Dist::isKnownAsParam(const std::string &s) const {
        return _params.find(s) != _params.end();
    }

    ex Dist::parseWithKnownSymbs(const std::string &s) const {
        str2sym syms = this->paramsAndVars();
        parser pars(syms);
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
        EXPECTS(!isKnownAsParam(p), "'variable' " + p.get_name() + " is already known as a parameter.");
        _vars[p.get_name()] = p;
    }

    void Dist::registerAsParam(const symbol &p) {
        EXPECTS(!isKnownAsVar(p), "'parameter' " + p.get_name() + " is already known as a variable.");
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
        
        parser pars;

        ex res = pars(expr);
        auto syms = pars.get_syms();

        lst names;
        lst values;
        for(const auto& var_val : syms) {
            EXPECTS(valuations.find(var_val.first) != valuations.end(), "Missing valuation for variable " + var_val.first);
            names.append(pars.get_syms().at(var_val.first));
            values.append(valuations.at(var_val.first));
        }
        res = res.subs(names, values);
        EXPECTS(res.info(info_flags::integer), "result of valuation is not an integer");
        return ex_to<numeric>(res).to_int();
    }

    uint64_t Dist::nthRoot(const uint64_t number, const uint64_t root) {
        ex power = pow(number, 1/ex(root));
        if (is_a<numeric>(power)) {
            numeric res = ex_to<numeric>(power);
            if (res.is_integer())
                return res.to_int();
        }
        throw std::invalid_argument("result is not an integer");
    }


/* public interface of dist */

    Dist::Dist(const std::string &gfStr) {
        parser pars_vars;
        ex e = pars_vars(gfStr);
        for (auto it = e.preorder_begin(); it != e.preorder_end(); ++it) {
            if (is_a<symbol>(*it)) {
                const auto &s = ex_to<symbol>(*it);
                if (!isSymbKnown(s)) {
                    _vars[s.get_name()] = s;
                }
            }
        }
        _gf=e;
        
        assumptions = {};
    }

    Dist::Dist(const std::string &gfStr, const std::vector<std::string> &paramsStr) {
        parser pars_vars_and_params;
        ex gf = pars_vars_and_params(gfStr); // might throw parse error
        registerAllUnknownSymbsAsVars(gf);
        for(const auto& pStr : paramsStr) {
            ex param = pars_vars_and_params(pStr); // might throw parse error
            const auto &pSym = ex_to<symbol>(param);
            deregister(pStr);
            registerAsParam(pSym);
        }
        _gf = gf;
        assumptions = {};
    }

    const ex &Dist::gf() const {
        return _gf;
    }

    const std::vector<std::shared_ptr<Assumption>> Dist::getAssumptions() const {
        return assumptions;
    }

    const Dist Dist::setAssumptions(const std::vector<std::shared_ptr<Assumption>>& assumptions) const {
        Dist res(_gf);
        res.assumptions = assumptions;
        return res;
    }

    const bool Dist::assumptionsHold(const std::vector<std::shared_ptr<Assumption>>& assumptions) const {
        for (auto assumption : assumptions) {
            bool isImplied = false;
            for (auto holds : this->assumptions) {
                if(holds->implies(assumption)) {
                    isImplied = true;
                    break;
                }
            }
            if(!isImplied)
                return false;
        }
        return true;
    }

    ex Dist::mass(const std::optional<std::set<std::string>> vars) const {
        ex res{_gf};
        /* TODO
         * normal() cancels common factors, this is needed since we are going to substitute stuff to avoid div by 0
         * encapsulate this somehow because substitution occurs elsewhere too
         */
        res = res.normal();
        lst names;
        lst values;
        if(!vars.has_value()) {
            for (const auto &v : _vars) {
                names.append(v.second);
                values.append(1);
            } 
        } else {
            for (auto v : vars.value()) {
                names.append(_vars.at(v));
                values.append(1);
            }
        }
        res = res.subs(names, values);
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
        Dist result(*this);
        result._gf = result._gf / result.mass();
        return result;
    }

    Dist Dist::marginal(const std::string &x) const {
        EXPECTS(!this->isKnownAsParam(x), "given symbol is considered a parameter");
        Dist result(*this);
        ex res{_gf.expand().normal()};
        lst names, values;
        for (const auto &v : _vars) {
            if (v.first != x) {
                names.append(v.second);
                values.append(1);
                result._vars.erase(v.first);
            }
        }
        res = res.subs(names, values);
        result._gf = res;
        return result;
    }

    Dist Dist::marginal(const std::vector<std::string>& vars) const {
        for (const auto &v : vars) {
            EXPECTS(!isKnownAsParam(v), "given symbol is considered a parameter");
        }
        Dist result(*this);
        ex res{_gf.normal()};
        // go over all variables v we currently know
        lst names, values;
        for (const auto &v : _vars) {
            // if v is not in the given list vars then project it away
            if (std::find(vars.begin(), vars.end(), v.first) == vars.end()) {
                names.append(v.second);
                values.append(1);
            }
        }
        res = res.subs(names, values);
        result._gf = res;
        return result;
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
        lst varEx;
        for (const auto &var : vars){
            EXPECTS(!isKnownAsParam(var), "given symbol is considered a parameter");
            if (isKnownAsVar(var)) {
                varEx.append(_vars.at(var));
            } else {
                varEx.append(symbol(var));
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
                varEx.append(symbol(var));
            }
        }
        return getTerms(varEx);
    }

    std::vector<Dist::State> Dist::getTerms() const {
        std::set<std::string> vars;
        for(auto var : _vars)
            vars.insert(var.first);
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

    std::vector<symbol> Dist::getSymbols() const {
        
        std::vector<symbol> res;

        for(const auto& item : this->paramsAndVars() ) {
            res.push_back(ex_to<symbol>(item.second));
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
        // basically encoding a skip statement.
        if(updatedVar == assignVar) {
            return *this;
        }
        
        // validate input
        Dist result = Dist(*this);
        ex exUpdatedVar = parseWithKnownSymbs(updatedVar);
        ex exAssignVar = parseWithKnownSymbs(assignVar);
        EXPECTS(!isKnownAsParam(updatedVar), "'variable' " + updatedVar + " is known as a parameter");
        EXPECTS(!isKnownAsParam(assignVar), "'variable' " + assignVar + " is a parameter");
        EXPECTS(is_a<symbol>(exAssignVar) || exAssignVar.info(info_flags::nonnegint), "The rhs mus be a variable or a non-negative integer, but was " + assignVar);
        EXPECTS(is_a<symbol>(exUpdatedVar), "'variable' " + updatedVar + " is not a symbol");
        
        if (!result.isKnownAsVar(ex_to<symbol>(exUpdatedVar))){
            result.registerAsVar(ex_to<symbol>(exUpdatedVar));
        }

        if(is_a<symbol>(exAssignVar)) {
            //EXPECTS(isKnownAsVar(assignVar), "cannot substitute " + updatedVar + " := " + assignVar + " -- " + assignVar + " is unknown.");
            if(!result.isKnownAsVar(ex_to<symbol>(exAssignVar))){
                result.registerAsVar(ex_to<symbol>(exAssignVar));
                result._gf = result._gf.subs(exUpdatedVar == exAssignVar);
            } else {
                result._gf = result._gf.subs(exUpdatedVar == 1).subs(exAssignVar == exAssignVar * exUpdatedVar);
            }
        } else {
            result._gf = (exAssignVar != 0 ? result._gf : result._gf.normal()).subs(exUpdatedVar == 1) * GiNaC::pow(exUpdatedVar, exAssignVar);
        }
        return result;
    }

    Dist Dist::updateSum(const std::string &updatedVar, const std::string &firstSummand, const std::string secondSummand) const {
        EXPECTS(!isKnownAsParam(updatedVar), "'variable' " + updatedVar + " is known as a parameter");

        auto exUpdatedVar = parseWithKnownSymbs(updatedVar);
        auto exFirstSummand = parseWithKnownSymbs(firstSummand);
        auto exSecondSummand = parseWithKnownSymbs(secondSummand);

        
        EXPECTS(is_a<symbol>(exUpdatedVar), "'variable' " + updatedVar + " is not a symbol");
        
        
        Dist res(*this);
        res.registerAllUnknownSymbsAsVars(exUpdatedVar);

        // we add two variables
        if(res.isKnownAsVar(firstSummand) && res.isKnownAsVar(secondSummand)) {
            if(exSecondSummand == exUpdatedVar) {
                ex tmp = exFirstSummand;
                exFirstSummand = exSecondSummand;
                exSecondSummand = tmp;
            }
            if(exFirstSummand == exUpdatedVar) {
                if(exSecondSummand == exUpdatedVar) {
                    res._gf = res._gf.subs(exUpdatedVar == GiNaC::pow(exUpdatedVar, 2));
                } else {
                    res._gf = res._gf.subs(exSecondSummand == exSecondSummand * exUpdatedVar);
                }
            } else {
                res._gf = res._gf.subs(exUpdatedVar == 1).subs(exFirstSummand == exFirstSummand * exUpdatedVar).subs(exSecondSummand == exSecondSummand * exUpdatedVar);
            }
            return res;
        }
        // we add a variable and a literal / parameter
        else if(res.isKnownAsVar(firstSummand) || res.isKnownAsVar(secondSummand)) {
            ex var = exFirstSummand;
            ex lit = exSecondSummand;
            if(res.isKnownAsVar(secondSummand)) {
                var = exSecondSummand;
                lit = exFirstSummand;
            }
            if(var != exUpdatedVar) {
                res._gf = res._gf.subs(exUpdatedVar == 1).subs(var == exUpdatedVar * var);
            }
            res._gf *= GiNaC::pow(exUpdatedVar, lit);
            return res;
        }
        // we add two literals / parameters
        else {
            res._gf = res._gf.subs(exUpdatedVar ==1) * GiNaC::pow(exUpdatedVar, exFirstSummand + exSecondSummand);
            return res;
        }
    }

    Dist Dist::updateProduct(const std::string &updatedVar, const std::string &firstFactor, const std::string &secondFactor, 
            const std::set<std::string> &variables, const bool &finite, const std::optional<std::string> &approx) const {
        EXPECTS(!isKnownAsParam(updatedVar), "'variable' " + updatedVar + " is known as a parameter");
        EXPECTS(!isKnownAsParam(firstFactor) && !isKnownAsParam(secondFactor), "cannot perform multiplication of parameters");
       
        auto exUpdatedVar = parseWithKnownSymbs(updatedVar);
        EXPECTS(is_a<symbol>(exUpdatedVar), "'variable' " + updatedVar + " is not a symbol");
        
        Dist res(*this);
        res.registerAllUnknownSymbsAsVars(exUpdatedVar);
        ex lhs = res.parseWithKnownSymbs(updatedVar);
        ex fac1 = res.parseWithKnownSymbs(firstFactor);
        ex fac2 = res.parseWithKnownSymbs(secondFactor);

        // we multiply two variables
        if(res.isKnownAsVar(firstFactor) && res.isKnownAsVar(secondFactor)) {
            // infinite distribution
            if(!finite) {
                auto marginal_l = res.marginal(firstFactor);
                auto finite_left = marginal_l.isPolynomial({firstFactor}) == troolean::TRUE;
                auto marginal_r = res.marginal(secondFactor);
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
                ex result = 0;
                for(auto state : finite.getTerms({finite_var})) {
                    // state.vals contains only one item
                    auto val = std::to_string(state.vals[finite_var]);
                    result += res.filterEq(finite_var, val)
                        .updateProduct(updatedVar, val, infinite_var, variables, other_finite, approx)._gf;
                }
                res._gf = result;
                return res;
            }
            // finite distribution
            else {
                ex result = res._gf;
                for(auto state : res.getTerms(variables)) {
                    ex term = res.parseWithKnownSymbs(state.toString());
                    result -= term;
                    term = term.subs(lhs == 1) * GiNaC::pow(lhs, state.vals.at(firstFactor) * state.vals.at(secondFactor));
                    result += term;
                }
                res._gf = result;
                return res;
            }
        }
        // we multiply a variable with a literal
        else if(res.isKnownAsVar(firstFactor) || res.isKnownAsVar(secondFactor)) {
            ex var, lit;
            if(res.isKnownAsVar(firstFactor)) {
                var = fac1;
                lit = fac2;
            } else {
                var = fac2;
                lit = fac1;
            }
            if(var == lhs) {
                res._gf = res._gf.subs(lhs == GiNaC::pow(lhs, lit));
                return res;
            } else {
                res._gf = res._gf.subs(lhs == 1).subs(var == var * GiNaC::pow(lhs, lit));
                return res;
            }
        }
        // we multiply two literals
        else {
            res._gf = res._gf.subs(lhs == 1) * GiNaC::pow(lhs, fac1 * fac2);
            return res;
        }
    }

    Dist Dist::updateSubtraction(const std::string updatedVar, const std::string subFrom, const std::string sub) const {
        EXPECTS(!isKnownAsParam(subFrom), "cannot subtract from parameter " + subFrom);
        EXPECTS(!isKnownAsParam(sub), "cannot subtract a parameter " + sub);

        EXPECTS(!isKnownAsParam(updatedVar), "'variable' " + updatedVar + " is known as a parameter");
        auto exUpdatedVar = parseWithKnownSymbs(updatedVar);
        EXPECTS(is_a<symbol>(exUpdatedVar), "'variable' " + updatedVar + " is not a symbol");
        
        Dist res(*this);
        res.registerAllUnknownSymbsAsVars(exUpdatedVar);
        ex lhs = res.parseWithKnownSymbs(updatedVar);
        ex sub1 = res.parseWithKnownSymbs(subFrom);
        ex sub2 = res.parseWithKnownSymbs(sub);

        ex result = res._gf;
        // we subtract a variable from another variable
        if(res.isKnownAsVar(sub) && res.isKnownAsVar(subFrom)) {
            if(res.isFiniteSupport() == troolean::TRUE) {
                for(auto state: res.getTerms()) {
                    std::cout << state.toString() << std::endl;
                    EXPECTS(state.vals[subFrom] >= state.vals[sub], "result can be negative");
                }
            }
            else {
                // TODO how to check for this case? We probably need limits
                // ex testFun = dist.marginal(updatedVar)._gf.normal().subs(lhs == 0);
            }

            if(sub2 == lhs) {
                if(sub1 == lhs) {
                    result = res._gf.subs(lhs == 1);
                } else {
                    result = res._gf.subs(lhs == GiNaC::pow(lhs, -1)).subs(sub1 == sub1 * lhs);
                }
            } else {
                if(sub1 != lhs) {
                    result = result.subs(lhs == 1).subs(sub1 == sub1 * lhs);
                }
                result = result.subs(sub2 == sub2 * GiNaC::pow(lhs, -1));
            }
        }
        //we subtract a literal from a variable
        else if(res.isKnownAsVar(subFrom)) {
            EXPECTS(res.filterLess(subFrom, sub) == Dist("0"), "result can be negative");
            if(lhs != sub1) {
                result = result.subs(lhs == 1).subs(sub1 == sub1 * lhs);
            }
            result = result * GiNaC::pow(lhs, -sub2);
        }
        // we subtract a variable from a literal
        else if(res.isKnownAsVar(sub)) {
            EXPECTS(res.filterGreater(sub, subFrom) == Dist("0"), "result can be negative");
            if(sub2 == lhs) {
                result = res._gf.subs(lhs == GiNaC::pow(lhs, -1)) * GiNaC::pow(lhs, sub1);
            } else {
                result = (res._gf.subs(lhs == 1) * GiNaC::pow(lhs, sub1)).subs(sub2 == sub2 * GiNaC::pow(lhs, -1));
            }
        }
        // we subtract two literals from each other
        else {
            EXPECTS(sub1 >= sub2, "result can be negative");
            result = res._gf.subs(lhs == 1) * GiNaC::pow(lhs, sub1 - sub2);
        }
        res._gf = result;
        return res;
            
    }

    Dist Dist::updateModulo(const std::string &updatedVar, const std::string &left, const std::string &right, 
            const std::set<std::string> &variables, const std::optional<bool> &finite, const std::optional<std::string> &approx) const {
        EXPECTS(!isKnownAsParam(updatedVar), "'variable' " + updatedVar + " is known as a parameter");
        EXPECTS(!isKnownAsParam(left) && !isKnownAsParam(right), "cannot perform modulo of parameters " + left + " " + right);
        
        auto exUpdatedVar = parseWithKnownSymbs(updatedVar);
        EXPECTS(is_a<symbol>(exUpdatedVar), "'variable' " + updatedVar + " is not a symbol");
        
        Dist res(*this);
        res.registerAllUnknownSymbsAsVars(exUpdatedVar);
        ex lhs = res.parseWithKnownSymbs(updatedVar);

        // this is not the same as finite.value_or, which will always evaluate isFiniteSupport()
        bool is_finite = finite.has_value() ? finite.value() : res.isFiniteSupport() == troolean::TRUE;

        // on finite distributions, iterate over all states
        if(is_finite) {
            ex result = 0;
            for(auto state : res.getTerms(variables)) {
                numeric leftVar, rightVar;
                if(res.isKnownAsVar(left)) {
                    leftVar = state.vals.at(left);
                } else {
                    leftVar = ex_to<numeric>(res.parseWithKnownSymbs(left));
                }
                if(res.isKnownAsVar(right)) {
                    rightVar = state.vals.at(right);
                } else {
                    rightVar = ex_to<numeric>(res.parseWithKnownSymbs(right));
                }
                result += res.parseWithKnownSymbs(state.toString()).subs(lhs == 1) * GiNaC::pow(lhs, mod(leftVar, rightVar));
            }
            return Dist(res);
        }
        // if the distribution is infinite and right is a variable, it needs to have finite range
        else if(res.isKnownAsVar(right)) {
            auto marginal_r = res.marginal(right);
            if(marginal_r.isPolynomial({right}) == troolean::FALSE) {
                if(!approx.has_value()) {
                    throw std::invalid_argument("rhs needs to have a finite marginal");
                }
                marginal_r = marginal_r.approximateUnilaterally(right, approx.value());
            }
            ex result = 0;
            for(auto state : marginal_r.getTerms({right})) {
                auto right_val = std::to_string(state.vals.at(right));
                result += res.filterEq(right, right_val)
                    .updateModulo(updatedVar, left, right_val, variables, is_finite ? true : std::optional<bool>(), approx)._gf;
            }
            res._gf = result;
            return res;
        }
        // if left is a variable, it doesn't have to have finite range
        else if(res.isKnownAsVar(left)) {
            auto marginal_l = res.marginal(left);
            ex result = 0;
            if(marginal_l.isPolynomial({left}) == troolean::TRUE) {
                for(auto state : marginal_l.getTerms({left})) {
                    auto left_val = std::to_string(state.vals.at(left));
                    result += res.filterEq(left, left_val)
                        .updateModulo(updatedVar, left_val, right, variables, is_finite ? true : std::optional<bool>(), approx)._gf;
                }
            } else {
                auto prog = res.arithmeticProgression(left, right);
                for(int i = 0; i < prog.size(); i++) {
                    result += prog.at(i)._gf.subs(lhs == 1) * GiNaC::pow(lhs, i);
                }
            }
            res._gf = result;
            return res;
        }
        // if both are not variables, simply compute the result
        else {
            numeric mod1 = ex_to<numeric>(parseWithKnownSymbs(left));
            numeric mod2 = ex_to<numeric>(parseWithKnownSymbs(right));
            return res.updateVar(updatedVar, std::to_string(GiNaC::to_int(mod(mod1, mod2))));
        }
    }

    Dist Dist::updateDivision(const std::string &updatedVar, const std::string &numerator, const std::string &denominator, const std::optional<std::string> &approx) const {
        EXPECTS(!isKnownAsParam(updatedVar), "'variable' " + updatedVar + " is known as a parameter");
        EXPECTS(!isKnownAsParam(numerator) && !isKnownAsParam(denominator), "cannot perform division of parameters " + numerator + " " + denominator);
        auto exUpdatedVar = parseWithKnownSymbs(updatedVar);
        EXPECTS(is_a<symbol>(exUpdatedVar), "'variable'" + updatedVar + " is not a symbol");
        
        Dist result(*this);
        result.registerAllUnknownSymbsAsVars(exUpdatedVar);
        ex lhs = result.parseWithKnownSymbs(updatedVar);
        ex div1 = result.parseWithKnownSymbs(numerator);
        ex div2 = result.parseWithKnownSymbs(denominator);

        ex res;

        if(result.isKnownAsVar(numerator) && is_a<numeric>(div2)) {
            if(result.assumptionsHold({std::make_shared<Divisible>(numerator, ex_to<numeric>(div2).to_int())})){
                if(lhs == div1) {
                    res = result._gf.subs(div1 == pow(div1, 1/div2));
                } else {
                    res = result._gf.subs(lhs == 1).subs(div1 == div1*pow(lhs, 1/div2));
                }
                result._gf = res;
                return result;
            }
        }

        Dist marginal_l, marginal_r;
        if(result.isKnownAsVar(numerator)) {
            marginal_l = result.marginal(numerator);
            if(marginal_l.isPolynomial({numerator}) == troolean::FALSE) {
                if(!approx.has_value()) {
                    throw std::invalid_argument("numerator must have a finite marginal");
                }
                marginal_l = marginal_l.approximateUnilaterally(numerator, approx.value());
            }
        } else {
            marginal_l = Dist(div1);
        }
        if(result.isKnownAsVar(denominator)) {
            marginal_r = result.marginal(denominator);
            if(marginal_r.isPolynomial({denominator}) == troolean::FALSE) {
                if(!approx.has_value()) {
                    throw std::invalid_argument("denominator must have a finite marginal");
                }
                marginal_r = marginal_r.approximateUnilaterally(denominator, approx.value());
            }
        } else {
            marginal_r = Dist(div2);
        }

        if(result.isKnownAsVar(numerator)) {
            for(auto state_l : marginal_l.getTerms({numerator})) {
                if(result.isKnownAsVar(denominator)) {
                    for(auto state_r : marginal_r.getTerms({denominator})) {
                        int val_l = state_l.vals.at(numerator);
                        int val_r = state_r.vals.at(denominator);
                        ex x = result.filterEq(numerator, std::to_string(val_l))
                            .filterEq(denominator, std::to_string(val_r))._gf;
                        EXPECTS(mod(val_l, val_r) == 0 || x == 0, "numerator is not always divisible by denominator");
                        res += x.subs(lhs == 1) * GiNaC::pow(lhs, val_l / val_r);
                    }
                } else {
                    int val_l = state_l.vals.at(numerator);
                    ex x = result.filterEq(numerator, std::to_string(val_l))._gf;
                    EXPECTS(mod(val_l, ex_to<numeric>(div2)) == 0 || x == 0, "numerator is not always divisible by denominator");
                    res += x.subs(lhs == 1) * GiNaC::pow(lhs, val_l / div2);
                }
            }
        } else {
            if(result.isKnownAsVar(denominator)) {
                for(auto state_r : marginal_r.getTerms({denominator})) {
                    int val_r = state_r.vals.at(denominator);
                    ex x = result.filterEq(denominator, std::to_string(val_r))._gf;
                    EXPECTS(mod(ex_to<numeric>(div1), val_r) == 0 || x == 0, "numerator is not always divisible by denominator");
                    res += x.subs(lhs == 1) * GiNaC::pow(lhs, div1 / val_r);
                } 
            } else {
                EXPECTS(mod(ex_to<numeric>(div1), ex_to<numeric>(div2)) == 0, "numerator is not always divisible by denominator");
                res = result._gf.subs(lhs == 1) * GiNaC::pow(lhs, div1 / div2);
            }
        }
        result._gf = res;
        return result;
    }

    Dist Dist::updatePower(const std::string &updatedVar, const std::string &base, const std::string &exp, const std::optional<std::string> &approx) const {
        EXPECTS(!isKnownAsParam(updatedVar), "'variable' " + updatedVar + " is known as a parameter");
        EXPECTS(!isKnownAsParam(base) && !isKnownAsParam(exp), "cannot perform exponentiation of parameters " + base + " " + exp);
        auto exUpdatedVar = parseWithKnownSymbs(updatedVar);
        EXPECTS(is_a<symbol>(exUpdatedVar), "'variable' " + updatedVar + " is not a symbol");
        
        Dist result(*this);
        result.registerAllUnknownSymbsAsVars(exUpdatedVar);
        ex lhs = result.parseWithKnownSymbs(updatedVar);
        ex pow1 = result.parseWithKnownSymbs(base);
        ex pow2 = result.parseWithKnownSymbs(exp);
        ex res = result._gf;

        // variable to the power of a variable
        if(result.isKnownAsVar(base) && result.isKnownAsVar(exp)) {
            auto marginal_l = result.marginal(base);
            auto marginal_r = result.marginal(exp);
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
                    ex x = result.filterEq(base, std::to_string(b)).filterEq(exp, std::to_string(e))._gf;
                    res -= x;
                    res += x.subs(lhs == 1) * GiNaC::pow(lhs, GiNaC::pow(b, e));
                }
            }
        }
        // variable to the power of a literal
        else if(result.isKnownAsVar(base)) {
            auto marginal = result.marginal(base);
            EXPECTS(marginal.isPolynomial({base}) == troolean::TRUE, 
                    "Can only perform exponentiation if the base has a finite marginal");

            for(auto state : marginal.getTerms({base})) {
                int b = state.vals.at(base);
                ex x = result.filterEq(base, std::to_string(b))._gf;
                res -= x;
                res += x.subs(lhs == 1) * GiNaC::pow(lhs, GiNaC::pow(b, pow2));
            }
        }
        // literal to the power of a variable
        else if(result.isKnownAsVar(exp)) {
            auto marginal = result.marginal(exp);
            EXPECTS(marginal.isPolynomial({exp}) == troolean::TRUE, 
                    "Can only perform exponentiation if the base has a finite marginal");

            for(auto state : marginal.getTerms({exp})) {
                int e = state.vals.at(exp);
                ex x = result.filterEq(exp, std::to_string(e))._gf;
                res -= x;
                res += x.subs(lhs == 1) * GiNaC::pow(lhs, GiNaC::pow(pow1, e));
            }
        }
        // literal to the power of a literal
        else {
            res = res.subs(lhs == 1) * GiNaC::pow(lhs, GiNaC::pow(pow1, pow2));
        }
        result._gf = res;
        return result;
    }

    Dist Dist::updateRoot(const std::string &updatedVar, const std::string &radicand, const std::string &index, const std::optional<std::string> &approx) const {
        EXPECTS(!isKnownAsParam(updatedVar), "lhs is known as a parameter");
        EXPECTS(is_a<symbol>(parseWithKnownSymbs(updatedVar)), "lhs is not a symbol");
        EXPECTS(!isKnownAsParam(radicand) && !isKnownAsParam(index), "cannot perform exponentiation of parameters");
        Dist result(*this);
        result.registerAllUnknownSymbsAsVars(possymbol(updatedVar));
        ex lhs = result.parseWithKnownSymbs(updatedVar);
        ex number = result.parseWithKnownSymbs(radicand);
        ex root = result.parseWithKnownSymbs(index);
        ex res;

        if(result.isKnownAsVar(radicand) || result.isKnownAsVar(index)) {
            std::vector<std::string> variables(2);
            if(result.isKnownAsVar(radicand))
                variables.push_back(radicand);
            else
                EXPECTS(is_a<numeric>(number), "radicand is neither a known symbol nor a numeric");
            if(result.isKnownAsVar(index))
                variables.push_back(index);
            else
                EXPECTS(is_a<numeric>(root), "index is neither a known symbol nor a numeric");
            auto marginal = result.marginal(variables);
            if(result.isKnownAsVar(radicand) && marginal.isPolynomial({radicand}) != troolean::TRUE) {
                EXPECTS(approx.has_value(), "radicand has infinite range");
                marginal = marginal.approximateUnilaterally(radicand, approx.value());
            }
            if(result.isKnownAsVar(index) && marginal.isPolynomial({index}) != troolean::TRUE) {
                EXPECTS(approx.has_value(), "index has infinite range");
                marginal = marginal.approximateUnilaterally(index, approx.value());
            }
            res = 0;
            for(auto state: marginal.getTerms()) {
                auto filtered = *this;
                int a, b;
                if(result.isKnownAsVar(radicand)) {
                    std::stringstream s;
                    s << state.vals[radicand];
                    filtered = filtered.filterEq(radicand, s.str());
                    a = state.vals[radicand];
                } else
                    a = ex_to<numeric>(number).to_int();
                if(result.isKnownAsVar(index)) {
                    std::stringstream s;
                    s << state.vals[index];
                    filtered = filtered.filterEq(index, s.str());
                    b = state.vals[index];
                } else
                    b = ex_to<numeric>(root).to_int();
                res += filtered._gf.subs(lhs == 1) * pow(lhs, nthRoot(a, b));
            }
        }
        else {
            EXPECTS(is_a<numeric>(number) && is_a<numeric>(root), "radicand or index is neither a known symbol nor a numeric");
            res = result._gf.subs(lhs == 1) * pow(lhs, nthRoot(ex_to<numeric>(number).to_int(), ex_to<numeric>(root).to_int()));
        }
        result._gf = res;
        return result;
    }

    Dist Dist::approximateUnilaterally(const std::string &variable, const std::string &prob_mass) const {
        ex mass = parseWithKnownSymbs(prob_mass).evalf();
        EXPECTS(is_a<numeric>(mass), "given probability mass is not a number");
        if(mass == 0) {
            return Dist(*this) * "0";
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
        Dist result(*this);
        result._gf = res;
        return result;
    }

    Dist Dist::updateIid(const std::string &x, const Dist &d, const std::string &y) const {
        // models assignments like x := iid(d, y).
        EXPECTS(!isKnownAsParam(x) && !isKnownAsParam(y), "given symbol is considered a parameter, not a variable");
        EXPECTS(d.mass() == 1, "the distribution wrapped into iid must not be a sub distribution");
        EXPECTS(d.isUnivariate() == troolean::TRUE, "the distribution wrapped into iid must be univariate");
        auto exX = parseWithKnownSymbs(x);
        auto exY = parseWithKnownSymbs(y);
        EXPECTS(is_a<symbol>(exX), "Cannot assign a value to " + x);
        EXPECTS(is_a<symbol>(exY) || exY.info(info_flags::nonnegint), "in " + x + ":= iid(D, " + y + "), " + y + " must be a symbol or a non-negative integer");

        Dist result(*this);

        auto substitute = d._gf;

        // make temp variable to actual x variable:
        for(auto sym : d._vars)
                substitute = substitute.subs(sym.second == exX);

        for(auto sym : d._params)
            if(!result.isKnownAsParam(ex_to<symbol>(sym.second)) && !result.isKnownAsVar(ex_to<symbol>(sym.second))){
                result.registerAsParam(ex_to<symbol>(sym.second));
            }

        // x is not known
        if(!result.isKnownAsVar(ex_to<symbol>(exX))){
            // x is a fresh variable.
            result.registerAsVar(ex_to<symbol>(exX));
            if(!result.isKnownAsVar(ex_to<symbol>(exY))){
                if(exY.info(info_flags::nonnegint)){
                    // y is just a number
                    result._gf *= pow(substitute, exY);
                }
                // y is a fresh variable (nothing to do in this case).
            } else {
                // y is a known variable
                result._gf = result._gf.subs(exY == (exY * substitute));
            }
        } else {
            // x is a known variable thus we might need to marginalize.
            if(!result.isKnownAsVar(ex_to<symbol>(exY))){
                if(exY.info(info_flags::nonnegint)){
                    // y is just a number -> marginalize
                    result._gf = result._gf.subs(exX == 1) * pow(substitute, exY);
                }
                // y is a fresh variable (nothing to do in this case).
            } else {
                // y is a known variable
                if (x == y){
                    // in this case we do not have to marginalize!
                    result._gf = result._gf.subs(exY == substitute);
                } else {
                    // z and y are different so we need to marginalize.
                    result._gf = result._gf.subs(exX == 1).subs(exY == (exY * substitute));
                }
                
            }
        }
        return result;
    }

    Dist Dist::filterLess(const std::string &x, const std::string &e) const {
        EXPECTS(!isKnownAsParam(x), "given symbol " + x + " is a parameter");
        ex eExp = parseWithKnownSymbs(e);
        EXPECTS(eExp.info(info_flags::nonnegint), "given expression must be non-neg int");

        int order = ex_to<numeric>(eExp).to_int();
        if(order == 0) {
            // nothing can be less than 0
            return Dist("0");
        }
        Dist result(*this);
        ex xExp = parseWithKnownSymbs(x);
        // todo why not use parseWithKnownSymbols?
        if(!result.isKnownAsVar(x)) {
            result.registerAsVar(ex_to<symbol>(xExp));
        } else {
            xExp = result._vars.at(x);
        }
        ex series = result._gf.series(xExp, order);
        result._gf = series_to_poly(series);
        return result;
    }

    Dist Dist::CoefficientIterator::next() {
        ex nextTerm = series_to_poly(_rest.series(_var == 0, _n));
        _n++;
        _rest = (_rest - nextTerm).normal();
        ex nextCoeff = nextTerm.subs(_var == 1);
        Dist result;
        result._gf = nextCoeff;
        result._vars = _vars;
        result._params = _params;
        return result;
    }

    Dist::CoefficientIterator Dist::coefficientIterator(const std::string &var) const {
        EXPECTS(isKnownAsVar(var), "can only expand series in known variables -- not in " + var);
        return CoefficientIterator(*this, ex_to<symbol>(_vars.at(var)));
    }

    Dist Dist::substituteParam(const std::string &p, const std::string &val) const {
        EXPECTS(!isKnownAsVar(p), "symbol " + p + " is considered a variable");
        if (!isKnownAsParam(p)) {
            return *this;
        }
        Dist result(*this);
        ex pSym = result._params[p];
        ex parsed = result.parseWithKnownSymbs(val);
        result.registerAllUnknownSymbsAsParams(parsed);
        // todo error if parsed contains known variables
        result._gf = result._gf.subs(pSym == parsed);
        return result;
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
        Dist result(lhs);
        auto ex_rhs = result.parseWithKnownSymbs(rhs);

        result.registerAllUnknownSymbsAsParams(ex_rhs);
        // todo error if parsed contained known variables
        // parsed contains just unknown symbs and known params
        // no error from now on
        result._gf *= ex_rhs;
        return result;
    }

} // namespace prodigy
