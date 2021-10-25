//
// Created by Tobias Winkler on 22.10.21.
//

#include <algorithm>
#include "ParametricDist.h"

/*
 * Visitor for collecting all symbols in a GiNaC::ex
 */
class gather_symbols_visitor : public visitor, public symbol::visitor {
    std::vector<symbol> symbols; // will contain the result
    std::vector<symbol> exclude; // symbols to exclude in result
    void visit(const symbol& s) override {
        // do not insert copies
        if(std::find(symbols.begin(), symbols.end(), s) == symbols.end()) {
            // do not insert explicitly excluded symbols
            if(exclude.empty() || std::find(exclude.begin(), exclude.end(), s) == exclude.end()) {
                symbols.push_back(s);
            }
        }
    }
public:
    const std::vector<symbol>& getResult() const {
        return symbols;
    }
    void setExclude(const std::vector<symbol>& toExclude) {
        exclude = toExclude;
    }
};

/* Public interface of ParametricDist */

ex ParametricDist::mass() const {
    ex res = _gf;
    for(const auto& v : getVariables()) {
        res = res.subs(v == 1);
    }
    return res;
}

ex ParametricDist::E(const symbol &var) const {
    if (hasParameter(var)) {
        throw std::invalid_argument("the given symbol is considered a parameter, not a variable");
    }
    ex res = _gf;
    for(const auto& v : getVariables()) {
        if(v != var) {
            res = res.subs(v == 1);
        }
    }
    return res.diff(var).subs(var == 1).normal();
}

ex ParametricDist::E(const std::string &var) const {
    if (hasParameter(var)) {
        throw std::invalid_argument("the given symbol is considered a parameter, not a variable");
    }
    try {
        symbol x = getVariableByName(var);
        return E(x);
    } catch (std::invalid_argument&) {
        // given variable is trivial, so expected value is 0
        return ex(0);
    }
}

std::vector<symbol> ParametricDist::getVariables() const {
    gather_symbols_visitor v;
    v.setExclude(params);
    _gf.traverse(v);
    return v.getResult();
}

symbol ParametricDist::getVariableByName(const std::string& v) const {
    std::vector<symbol> vars(getVariables());
    for(const auto& s : vars) {
        if(s.get_name() == v)
            return s;
    }
    throw std::invalid_argument("gf does not contain this variable");
}

bool ParametricDist::hasNonTrivialVariable(const symbol & s) const {
    std::vector<symbol> vars(getVariables());
    if(std::find(vars.begin(), vars.end(), s) != vars.end()) {
        return true;
    }
    return false;
}

bool ParametricDist::hasNonTrivialVariable(const std::string & s) const {
    try {
        getVariableByName(s);
        return true;
    } catch (std::invalid_argument&) {
        return false;
    }
}

symbol ParametricDist::getParameterByName(const std::string &v) const {
    for(const auto& s : params) {
        if(s.get_name() == v)
            return s;
    }
    throw std::invalid_argument("given symbol is no parameter");
}

bool ParametricDist::hasParameter(const std::string &p) const {
    try {
        getParameterByName(p);
        return true;
    } catch (std::invalid_argument&) {
        return false;
    }
}

bool ParametricDist::hasParameter(const symbol& p) const {
    if(std::find(params.begin(), params.end(), p) != params.end()) {
        return true;
    }
    return false;
}

/* operators */

ParametricDist ParametricDist::operator*(const std::string& rhs) const {
    try {
        // try to convert the given string to a numeric
        numeric scalar(rhs.c_str());
        return *this * scalar;
    } catch (std::exception& e) {
        // could not convert to numeric, try to convert to expression
        parser pars;
        ex parsed = pars(rhs);
        if (parsed.info(info_flags::symbol)) {
            // we have to decide if we interpret the symbol as a variable or a parameter
            if(hasParameter(rhs)) {
                return ParametricDist(parsed * _gf, params);
            }
            if(hasNonTrivialVariable(rhs)) {
                return ParametricDist(parsed * _gf, params);
            }
            // neither a parameter nor a non-trivial variable - interpret as parameter by convention
            ParametricDist res(*this);
            symbol newParam = symbol(rhs);
            res.params.push_back(newParam); // add the new parameter
            res._gf = newParam * res._gf;
            return res;
        } else {
            // only expressions that are a symbol are currently allowed
            throw prodigy::not_implemented_error();
        }
    } catch (parse_error& e) {
        throw std::invalid_argument("cannot interpret given string as number or expression");
    }
}






