//
// Created by Tobias Winkler on 22.10.21.
//

#include <algorithm>
#include "Dist.h"
#include "macros.h"

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

/* Public interface of Dist */

ex Dist::mass() const {
    ex res = _gf;
    for(const auto& v : getVariables()) {
        res = res.subs(v == 1);
    }
    return res;
}

ex Dist::E(const symbol &var) const {
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

ex Dist::E(const std::string &var) const {
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

std::vector<symbol> Dist::getVariables() const {
    std::vector<symbol> res;
    for(auto it = _gf.preorder_begin(); it != _gf.preorder_end(); ++it) {
        if(is_a<symbol>(*it)) {
            const auto& s = ex_to<symbol>(*it);
            if (std::find(res.begin(), res.end(), s) == res.end()) {
                if (std::find(_params.begin(), _params.end(), s) == _params.end()) {
                    res.push_back(s);
                }
            }
        }
    }
    return res;
}

symbol Dist::getVariableByName(const std::string& v) const {
    std::vector<symbol> vars(getVariables());
    for(const auto& s : vars) {
        if(s.get_name() == v)
            return s;
    }
    throw std::invalid_argument("gf does not contain this variable");
}

bool Dist::hasNonTrivialVariable(const symbol & s) const {
    std::vector<symbol> vars(getVariables());
    if(std::find(vars.begin(), vars.end(), s) != vars.end()) {
        return true;
    }
    return false;
}

bool Dist::hasNonTrivialVariable(const std::string & s) const {
    try {
        getVariableByName(s);
        return true;
    } catch (std::invalid_argument&) {
        return false;
    }
}

symbol Dist::getParameterByName(const std::string &v) const {
    for(const auto& s : _params) {
        if(s.get_name() == v)
            return s;
    }
    throw std::invalid_argument("given symbol is no parameter");
}

bool Dist::hasParameter(const std::string &p) const {
    try {
        getParameterByName(p);
        return true;
    } catch (std::invalid_argument&) {
        return false;
    }
}

bool Dist::hasParameter(const symbol& p) const {
    if(std::find(_params.begin(), _params.end(), p) != _params.end()) {
        return true;
    }
    return false;
}

Dist Dist::update(const std::string &v, const std::string &e) const {
    symtab tab = makeSymtab();
    parser pars(tab); // pass symtab of this distribution to parser so that it links symbols correctly
    ex parsed = pars(e);
    symbol vsymb;
    if(tab.find(v) != tab.end()) {
        // v is an existing symbol
        vsymb = ex_to<symbol>(tab[v]);
    } else {
        vsymb = symbol(v);
    }
    return update(vsymb, parsed);
}


Dist Dist::update(const symbol &v, const ex &e) const {
    /*
     * algorithm for evaluating given assignment v = term_1 + ... + term_n :
     *
     * 1) tmp = 0
     * 2) tmp += term_1
     *    ...
     *    tmp += term_n
     * 3) v = tmp
     *
     */
    PRECOND(!hasParameter(v), "parameters cannot be updated");
    PRECOND(e.info(info_flags::polynomial), "only polynomial arithmetic is supported in updates");
    ex pol = e.expand();
    Dist res {*this};
    // 1) tmp = 0
    symbol tmp;
    // 2) tmp += term_1 ... tmp += term_n
    if(is_a<numeric>(pol) || is_a<symbol>(pol) || is_a<mul>(pol)) {
        // e.g. pol = 1 || pol = x || pol = 2*y
        res = res.incrAssgn(tmp, pol);
    } else {
        assert(is_a<add>(pol));
        // e.g. pol = 1 + x
        for(const auto& term : pol) {
            res = res.incrAssgn(tmp, term);
        }
    }
    // 3) v = tmp
    res._gf = res._gf.subs(v == 1);
    res._gf = res._gf.subs(tmp == v);
    return res;
}



/* operators */

bool Dist::operator==(const Dist &rhs) const {
    return _gf == rhs._gf && _params == rhs._params;
}


Dist Dist::operator*(const std::string& rhs) const {
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
                return Dist {parsed * _gf, _params};
            }
            if(hasNonTrivialVariable(rhs)) {
                return Dist {parsed * _gf, _params};
            }
            // neither a parameter nor a non-trivial variable - interpret as parameter by convention
            Dist res(*this);
            symbol newParam = symbol(rhs);
            res._params.push_back(newParam); // add the new parameter
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

/* private member functions */

Dist Dist::incrAssgn(const symbol &v, const ex &term) const {
    Dist res {*this};
    // v += k, with k a non-negative int
    if(term.info(info_flags::nonnegint)) {
        res._gf *= pow(v, term);
    }
    // v += y, with y a variable
    else if(term.info(info_flags::symbol)) {
        res._gf = res._gf.subs(term == v * term);
    }
    else {
        throw prodigy::not_implemented_error{};
    }
    return res;
}

symtab Dist::makeSymtab() const {
    symtab res;
    for(auto it = _gf.preorder_begin(); it != _gf.preorder_end(); ++it) {
        if(is_a<symbol>(*it)) {
            const auto& s = ex_to<symbol>(*it);
            res[s.get_name()] = s;
        }
    }
    for(const auto& p : _params) {
        res[p.get_name()] = p;
    }
    return res;
}




