//
// Created by Tobias Winkler on 22.10.21.
//

#include <algorithm>
#include <numeric>
#include "Dist.h"
#include "macros.h"
#include <cassert>
#include <cctype>

using namespace GiNaC;

namespace prodigy {

    /*
     * substitutes the given relation(s), e.g. x == 1, into e
     * a gf may have a removable singularity there, e.g. x*(x^3-1)/(3*(x-1)) at x = 1, which makes GiNaC divide by zero;
     * in that case the gf is cancelled first (normal()) and substituted again, a genuine pole still throws
     */
    static ex subsCancelled(const ex &e, const lst &symbols, const lst &values) {
        try {
            return e.subs(symbols, values);
        } catch (const std::domain_error &) {
            return e.normal().subs(symbols, values);
        }
    }

    static ex subsCancelled(const ex &e, const ex &relations) {
        try {
            return e.subs(relations);
        } catch (const std::domain_error &) {
            return e.normal().subs(relations);
        }
    }

    static ex atOne(const ex &e, const ex &var) {
        return subsCancelled(e, var == 1);
    }

    /*
     * returns the sum of all terms of e of degree < order in x, where e is a gf (a power series in x)
     * GiNaC's series does not truncate products of polynomials, e.g. the series of (1+x)*(1+x^5) to order 3
     * is the whole product, so we drop the terms of degree >= order ourselves
     */
    static ex truncatedSeries(const ex &e, const ex &x, int order) {
        ex poly = series_to_poly(e.series(x == 0, order));
        // for a non-expanded polynomial, degree() is an upper bound of the degree
        if (poly.degree(x) < order) {
            return poly;
        }
        poly = poly.expand();
        ex res = 0;
        for (int k = 0; k < order; k++) {
            res += poly.coeff(x, k) * GiNaC::pow(x, k);
        }
        return res;
    }

    /*
     * multisection of a gf that is a rational function in x: computes the sections f_r = sum_{k = r mod q} [x^k]f x^k,
     * r = 0..q-1, exactly, without roots of unity (issue #10)
     *
     * Let f = P/Q with polynomials P, Q in x, and let w be a primitive q-th root of unity. Then
     * N(x) = prod_{m=0}^{q-1} Q(w^m x) is, up to sign, the resultant Res_y(Q(y), y^q - x^q), so a polynomial in x^q
     * whose coefficients lie in the coefficient domain of Q. N/Q = prod_{m=1}^{q-1} Q(w^m x) is a polynomial, hence so
     * is M = P*N/Q, and f = M/N. Since N(w x) = N(x),
     *     f_r(x) = 1/q sum_m w^(-m r) f(w^m x) = (1/q sum_m w^(-m r) M(w^m x)) / N(x) = section_r(M)/N,
     * where section_r(M) collects the terms of M whose degree in x is r modulo q. This is an identity of rational
     * functions, so the sections are equal to the sections of f as power series.
     *
     * Subexpressions that are not rational (e.g. exp(-l) or sqrt(p)) are replaced by temporary symbols first; if one of
     * them contains x, f is not rational in x and false is returned. false is also returned if the exact division
     * fails, which cannot happen mathematically; the caller then falls back to the roots of unity.
     */
    static bool rationalMultisection(const ex &f, const symbol &x, int q, std::vector<ex> &sections) {
        exmap repl;
        const ex rational = f.to_rational(repl);
        for (const auto &r : repl) {
            if (r.first.has(x) || r.second.has(x)) {
                return false;
            }
        }
        const ex nd = rational.normal().numer_denom();
        const ex P = nd.op(0).expand();
        const ex Q = nd.op(1).expand();
        if (!P.is_polynomial(x) || !Q.is_polynomial(x)) {
            return false;
        }

        ex N, M;
        if (!Q.has(x)) {
            N = Q;
            M = P;
        } else {
            const symbol y;
            N = resultant(Q.subs(x == y), GiNaC::pow(y, q) - GiNaC::pow(x, q), y).expand();
            if (!divide((P * N).expand(), Q, M)) {
                return false;
            }
            M = M.expand();
        }

        const int deg = M.degree(x);
        std::vector<ex> numerators(q, 0);
        for (int k = M.ldegree(x); k <= deg; k++) {
            numerators[k % q] += M.coeff(x, k) * GiNaC::pow(x, k);
        }
        sections.clear();
        sections.reserve(q);
        for (const ex &numerator : numerators) {
            sections.push_back((numerator / N).normal().subs(repl));
        }
        return true;
    }

/* static (global) portion of Dist */

    std::string Dist::exactDecimals(const std::string &s) {
        // mirrors the tokenizer of GiNaC's parser (ginac/parser/lexer.cpp):
        // identifiers are a letter followed by letters, digits and '_', numbers are [0-9.]+ with an optional
        // exponent [eE][+-]?[0-9]+, and '#' starts a comment until the end of the line
        const auto isDigit = [](char c) { return std::isdigit(static_cast<unsigned char>(c)) != 0; };
        const auto isAlpha = [](char c) { return std::isalpha(static_cast<unsigned char>(c)) != 0; };
        const auto isAlnum = [](char c) { return std::isalnum(static_cast<unsigned char>(c)) != 0; };
        const long maxExponent = 100000;

        std::string res;
        res.reserve(s.size());
        std::size_t i = 0;
        const std::size_t n = s.size();
        while (i < n) {
            const char c = s[i];
            std::size_t j = i + 1;
            if (isAlpha(c)) {
                while (j < n && (isAlnum(s[j]) || s[j] == '_')) j++;
                res.append(s, i, j - i);
            } else if (c == '#') {
                j = s.find_first_of("\n\r", i);
                if (j == std::string::npos) j = n;
                res.append(s, i, j - i);
            } else if (isDigit(c) || c == '.') {
                std::string digits; // all digits of the mantissa, without the dot
                std::size_t fracDigits = 0;
                bool dot = false;
                for (j = i; j < n && (isDigit(s[j]) || s[j] == '.'); j++) {
                    if (s[j] == '.') {
                        if (dot) throw std::invalid_argument("malformed number in '" + s + "'");
                        dot = true;
                    } else {
                        digits += s[j];
                        if (dot) fracDigits++;
                    }
                }
                if (digits.empty()) throw std::invalid_argument("malformed number in '" + s + "'");
                bool hasExponent = false;
                bool negativeExponent = false;
                std::string expDigits;
                if (j < n && (s[j] == 'e' || s[j] == 'E')) {
                    hasExponent = true;
                    j++;
                    if (j < n && (s[j] == '+' || s[j] == '-')) {
                        negativeExponent = s[j] == '-';
                        j++;
                    }
                    while (j < n && isDigit(s[j])) expDigits += s[j++];
                    if (expDigits.empty()) throw std::invalid_argument("malformed exponent in '" + s + "'");
                }
                if (!dot && !hasExponent) {
                    // integer literals are exact in GiNaC
                    res.append(s, i, j - i);
                } else {
                    expDigits.erase(0, std::min(expDigits.find_first_not_of('0'), expDigits.size()));
                    if (expDigits.size() > 6 || (!expDigits.empty() && std::stol(expDigits) > maxExponent)) {
                        throw std::invalid_argument("decimal exponent too large in '" + s + "'");
                    }
                    long exponent = expDigits.empty() ? 0 : std::stol(expDigits);
                    if (negativeExponent) exponent = -exponent;
                    exponent -= static_cast<long>(fracDigits);
                    numeric value = numeric(digits.c_str()) * GiNaC::pow(numeric(10), numeric(exponent));
                    std::ostringstream str;
                    str << '(' << value << ')';
                    res += str.str();
                }
            } else {
                res += c;
            }
            i = j;
        }
        return res;
    }

    bool Dist::hasInexactNumbers(const ex &e) {
        for (auto it = e.preorder_begin(); it != e.preorder_end(); ++it) {
            if (is_a<numeric>(*it) && !ex_to<numeric>(*it).is_crational()) {
                return true;
            }
        }
        return false;
    }

// initialize static variables of class Dist

    std::vector<std::string> Dist::find_symbols(const std::string &s) {
        // this function does not use the symbol registry, as opposed to Dist::getSymbols()
        // to use the modulo operator %, we need to replace it with something GiNaC can parse
        std::string str = s;
        std::replace(str.begin(), str.end(), '%', '/');

        GiNaC::parser p;
        ex parsed = p(exactDecimals(str));
        std::vector<std::string> res;
        for(auto i : p.get_syms()) {
            res.push_back(i.first);
        }
        return res;
    }

    ex Dist::parse(const std::string &s, str2sym &table, bool unknownAreVariables) {
        parser pars(table);
        ex e = pars(exactDecimals(s)); // might throw parse error
        exmap renaming;
        for (const auto &name : pars.get_syms()) {
            if (table.find(name.first) == table.end()) {
                // GiNaC's parser creates plain symbols for unknown names, variables are possymbols
                ex sym = name.second;
                if (unknownAreVariables) {
                    sym = possymbol(name.first);
                    renaming[name.second] = sym;
                }
                table[name.first] = sym;
            }
        }
        return renaming.empty() ? e : e.subs(renaming, subs_options::no_pattern);
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

    ex Dist::parseWithKnownSymbs(const std::string &s, bool unknownAreVariables) const {
        str2sym syms = this->paramsAndVars();
        return parse(s, syms, unknownAreVariables); // might throw parse error
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

    ex Dist::linkSymbs(const ex &e) {
        exmap renaming;
        for (auto it = e.preorder_begin(); it != e.preorder_end(); ++it) {
            if (is_a<symbol>(*it)) {
                const auto &s = ex_to<symbol>(*it);
                if (!isSymbKnown(s)) {
                    // variables are possymbols
                    if (is_a<possymbol>(s)) {
                        _vars[s.get_name()] = s;
                        continue;
                    }
                    _vars[s.get_name()] = possymbol(s.get_name());
                }
                const ex &known = isKnownAsParam(s) ? _params.at(s.get_name()) : _vars.at(s.get_name());
                if (!known.is_equal(s)) {
                    // same name, but a different GiNaC symbol
                    renaming[s] = known;
                }
            }
        }
        return renaming.empty() ? e : e.subs(renaming, subs_options::no_pattern);
    }

    void Dist::setGf(const ex &gf) {
        _gf = linkSymbs(gf);
    }

    Dist Dist::withGf(const ex &gf) const {
        Dist res;
        res._params = _params;
        res._vars = _vars;
        res.setGf(gf);
        return res;
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

        ex res = pars(exactDecimals(expr));
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
        str2sym table;
        ex e = parse(gfStr, table);
        registerAllUnknownSymbsAsVars(e);
        setGf(e);
        
        assumptions = {};
    }

    Dist::Dist(const std::string &gfStr, const std::vector<std::string> &paramsStr) {
        str2sym table;
        for(const auto& pStr : paramsStr) {
            ex param = parse(pStr, table, false); // might throw parse error
            EXPECTS(is_a<symbol>(param), "parameter " + pStr + " is not a symbol");
            registerAsParam(ex_to<symbol>(param));
        }
        ex gf = parse(gfStr, table); // might throw parse error
        registerAllUnknownSymbsAsVars(gf);
        setGf(gf);
        assumptions = {};
    }

    const ex &Dist::gf() const {
        return _gf;
    }

    const std::vector<std::shared_ptr<Assumption>> Dist::getAssumptions() const {
        return assumptions;
    }

    const Dist Dist::setAssumptions(const std::vector<std::shared_ptr<Assumption>>& assumptions) const {
        Dist res(*this);
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
        res = subsCancelled(res, names, values);
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
        result.setGf(result._gf / result.mass());
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
        result.setGf(res);
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
        result.setGf(res);
        return result;
    }

    troolean Dist::isZero() const {
        // todo look deeper into what is the best way in ginac for testing whether expression is zero
        // floating-point numbers are rounded, so neither zero nor non-zero is certain
        if (hasInexactNumbers(_gf)) {
            return troolean::UNKNOWN;
        }
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
        for(const auto& var : _vars) {
            vars.append(var.second);
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
                term = atOne(term, var);
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

        const int q = ex_to<numeric>(exMod).to_int();
        std::vector<Dist> res;
        res.reserve(q);

        std::vector<ex> sections;
        if (rationalMultisection(_gf, ex_to<symbol>(exVar), q, sections)) {
            for (const ex &section : sections) {
                res.push_back(withGf(section));
            }
            return res;
        }

        // the gf is not rational in var: use the roots of unity
        ex root = exp(2 * Pi * I / exMod);
        for(int rem = 0; rem < exMod; rem++) {
            ex psum = 0;
            for(int m = 0; m < exMod; m++) {
                psum += GiNaC::pow(root, -m * rem) * 
                    _gf.subs(exVar == GiNaC::pow(root, m) * exVar);
            }
            res.push_back(withGf(psum/exMod));
        }
        return res;
    }

    std::vector<ex> Dist::getSymbols() const {
        
        std::vector<ex> res;

        for(const auto& item : this->paramsAndVars() ) {
            res.push_back(item.second);
        }
        return res;
    }


    ex Dist::E(const std::string &x) const {
        EXPECTS(!isKnownAsParam(x), "given symbol is considered a parameter");
        ex parsed = parseWithKnownSymbs(x); // might throw parser error
        exset symbols;
        for (auto it = parsed.preorder_begin(); it != parsed.preorder_end(); ++it) {
            if (is_a<symbol>(*it)) {
                symbols.insert(*it);
            }
        }
        lst vars;
        for (const auto &sym : symbols) {
            const std::string &name = ex_to<symbol>(sym).get_name();
            if (isKnownAsVar(name)) {
                vars.append(sym);
            } else {
                // a symbol that is neither a variable nor a parameter has no known value: a result would be wrong
                EXPECTS(isKnownAsParam(name), "given symbol " + name + " is neither a variable nor a parameter");
            }
        }
        if (!parsed.is_polynomial(vars) || hasInexactNumbers(parsed)) {
            throw not_implemented_error("expected value only computable for polynomial expressions in the variables "
                                        "with exact coefficients, got " + x);
        }
        ex poly = parsed.expand();
        if (vars.nops() == 0) {
            return poly;
        }

        // E[c * x1^k1 * ... * xn^kn] = c * ((x1 d/dx1)^k1 ... (xn d/dxn)^kn G)(1, ..., 1) for the marginal gf G
        std::vector<std::string> varNames;
        for (const auto &v : vars) {
            varNames.push_back(ex_to<symbol>(v).get_name());
        }
        const ex gf = marginal(varNames)._gf;
        lst ones;
        for (size_t i = 0; i < vars.nops(); i++) {
            ones.append(1);
        }
        std::vector<ex> terms;
        if (is_a<add>(poly)) {
            terms.assign(poly.begin(), poly.end());
        } else {
            terms.push_back(poly);
        }
        ex res = 0;
        for (const ex &term : terms) {
            ex coeff = term;
            ex moment = gf;
            for (const auto &v : vars) {
                const symbol &sym = ex_to<symbol>(v);
                const int k = coeff.degree(sym);
                coeff = coeff.coeff(sym, k);
                for (int i = 0; i < k; i++) {
                    moment = (sym * moment.diff(sym)).normal();
                }
            }
            res += coeff * moment.subs(vars, ones);
        }
        return res.normal();
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
                result.setGf(result._gf.subs(exUpdatedVar == exAssignVar));
            } else {
                result.setGf(atOne(result._gf, exUpdatedVar).subs(exAssignVar == exAssignVar * exUpdatedVar));
            }
        } else {
            result.setGf(atOne(exAssignVar != 0 ? result._gf : result._gf.normal(), exUpdatedVar) * GiNaC::pow(exUpdatedVar, exAssignVar));
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
                    res.setGf(res._gf.subs(exUpdatedVar == GiNaC::pow(exUpdatedVar, 2)));
                } else {
                    res.setGf(res._gf.subs(exSecondSummand == exSecondSummand * exUpdatedVar));
                }
            } else {
                res.setGf(atOne(res._gf, exUpdatedVar).subs(exFirstSummand == exFirstSummand * exUpdatedVar).subs(exSecondSummand == exSecondSummand * exUpdatedVar));
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
                res.setGf(atOne(res._gf, exUpdatedVar).subs(var == exUpdatedVar * var));
            }
            res.setGf(res._gf * GiNaC::pow(exUpdatedVar, lit));
            return res;
        }
        // we add two literals / parameters
        else {
            res.setGf(atOne(res._gf, exUpdatedVar) * GiNaC::pow(exUpdatedVar, exFirstSummand + exSecondSummand));
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
                res.setGf(result);
                return res;
            }
            // finite distribution
            else {
                ex result = res._gf;
                for(auto state : res.getTerms(variables)) {
                    ex term = res.parseWithKnownSymbs(state.toString());
                    result -= term;
                    term = atOne(term, lhs) * GiNaC::pow(lhs, state.vals.at(firstFactor) * state.vals.at(secondFactor));
                    result += term;
                }
                res.setGf(result);
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
                res.setGf(res._gf.subs(lhs == GiNaC::pow(lhs, lit)));
                return res;
            } else {
                res.setGf(atOne(res._gf, lhs).subs(var == var * GiNaC::pow(lhs, lit)));
                return res;
            }
        }
        // we multiply two literals
        else {
            res.setGf(atOne(res._gf, lhs) * GiNaC::pow(lhs, fac1 * fac2));
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
                    EXPECTS(state.vals[subFrom] >= state.vals[sub], "result can be negative");
                }
            }
            else {
                // TODO how to check for this case? We probably need limits
                // ex testFun = dist.marginal(updatedVar)._gf.normal().subs(lhs == 0);
            }

            if(sub2 == lhs) {
                if(sub1 == lhs) {
                    result = atOne(res._gf, lhs);
                } else {
                    result = res._gf.subs(lhs == GiNaC::pow(lhs, -1)).subs(sub1 == sub1 * lhs);
                }
            } else {
                if(sub1 != lhs) {
                    result = atOne(result, lhs).subs(sub1 == sub1 * lhs);
                }
                result = result.subs(sub2 == sub2 * GiNaC::pow(lhs, -1));
            }
        }
        //we subtract a literal from a variable
        else if(res.isKnownAsVar(subFrom)) {
            EXPECTS(res.filterLess(subFrom, sub) == Dist("0"), "result can be negative");
            if(lhs != sub1) {
                result = atOne(result, lhs).subs(sub1 == sub1 * lhs);
            }
            result = result * GiNaC::pow(lhs, -sub2);
        }
        // we subtract a variable from a literal
        else if(res.isKnownAsVar(sub)) {
            EXPECTS(res.filterGreater(sub, subFrom) == Dist("0"), "result can be negative");
            if(sub2 == lhs) {
                result = res._gf.subs(lhs == GiNaC::pow(lhs, -1)) * GiNaC::pow(lhs, sub1);
            } else {
                result = (atOne(res._gf, lhs) * GiNaC::pow(lhs, sub1)).subs(sub2 == sub2 * GiNaC::pow(lhs, -1));
            }
        }
        // we subtract two literals from each other
        else {
            EXPECTS(sub1 >= sub2, "result can be negative");
            result = atOne(res._gf, lhs) * GiNaC::pow(lhs, sub1 - sub2);
        }
        res.setGf(result);
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
                result += atOne(res.parseWithKnownSymbs(state.toString()), lhs) * GiNaC::pow(lhs, mod(leftVar, rightVar));
            }
            return res.withGf(result);
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
            res.setGf(result);
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
                    result += atOne(prog.at(i)._gf, lhs) * GiNaC::pow(lhs, i);
                }
            }
            res.setGf(result);
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
                    res = atOne(result._gf, lhs).subs(div1 == div1*pow(lhs, 1/div2));
                }
                result.setGf(res);
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
            marginal_l = result.withGf(div1);
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
            marginal_r = result.withGf(div2);
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
                        res += atOne(x, lhs) * GiNaC::pow(lhs, val_l / val_r);
                    }
                } else {
                    int val_l = state_l.vals.at(numerator);
                    ex x = result.filterEq(numerator, std::to_string(val_l))._gf;
                    EXPECTS(mod(val_l, ex_to<numeric>(div2)) == 0 || x == 0, "numerator is not always divisible by denominator");
                    res += atOne(x, lhs) * GiNaC::pow(lhs, val_l / div2);
                }
            }
        } else {
            if(result.isKnownAsVar(denominator)) {
                for(auto state_r : marginal_r.getTerms({denominator})) {
                    int val_r = state_r.vals.at(denominator);
                    ex x = result.filterEq(denominator, std::to_string(val_r))._gf;
                    EXPECTS(mod(ex_to<numeric>(div1), val_r) == 0 || x == 0, "numerator is not always divisible by denominator");
                    res += atOne(x, lhs) * GiNaC::pow(lhs, div1 / val_r);
                } 
            } else {
                EXPECTS(mod(ex_to<numeric>(div1), ex_to<numeric>(div2)) == 0, "numerator is not always divisible by denominator");
                res = atOne(result._gf, lhs) * GiNaC::pow(lhs, div1 / div2);
            }
        }
        result.setGf(res);
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
                    res += atOne(x, lhs) * GiNaC::pow(lhs, GiNaC::pow(b, e));
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
                res += atOne(x, lhs) * GiNaC::pow(lhs, GiNaC::pow(b, pow2));
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
                res += atOne(x, lhs) * GiNaC::pow(lhs, GiNaC::pow(pow1, e));
            }
        }
        // literal to the power of a literal
        else {
            res = atOne(res, lhs) * GiNaC::pow(lhs, GiNaC::pow(pow1, pow2));
        }
        result.setGf(res);
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
                res += atOne(filtered._gf, lhs) * pow(lhs, nthRoot(a, b));
            }
        }
        else {
            EXPECTS(is_a<numeric>(number) && is_a<numeric>(root), "radicand or index is neither a known symbol nor a numeric");
            res = atOne(result._gf, lhs) * pow(lhs, nthRoot(ex_to<numeric>(number).to_int(), ex_to<numeric>(root).to_int()));
        }
        result.setGf(res);
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
        result.setGf(res);
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

        // the parameters of d are parameters of the result; setGf replaces d's parameter symbols by the result's
        for(const auto& param : d._params) {
            EXPECTS(!result.isKnownAsVar(param.first), "parameter " + param.first + " is known as variable.");
            if(!result.isKnownAsParam(param.first)) {
                result.registerAsParam(ex_to<symbol>(param.second));
            }
        }

        // x is not known
        if(!result.isKnownAsVar(ex_to<symbol>(exX))){
            // x is a fresh variable.
            result.registerAsVar(ex_to<symbol>(exX));
            if(!result.isKnownAsVar(y)){
                if(exY.info(info_flags::nonnegint)){
                    // y is just a number
                    result.setGf(result._gf * pow(substitute, exY));
                }
                // y is a fresh variable (nothing to do in this case).
            } else {
                // y is a known variable
                result.setGf(result._gf.subs(exY == (exY * substitute)));
            }
        } else {
            // x is a known variable thus we might need to marginalize.
            if(!result.isKnownAsVar(y)){
                if(exY.info(info_flags::nonnegint)){
                    // y is just a number -> marginalize
                    result.setGf(atOne(result._gf, exX) * pow(substitute, exY));
                }
                // y is a fresh variable (nothing to do in this case).
            } else {
                // y is a known variable
                if (x == y){
                    // in this case we do not have to marginalize!
                    result.setGf(result._gf.subs(exY == substitute));
                } else {
                    // z and y are different so we need to marginalize.
                    result.setGf(atOne(result._gf, exX).subs(exY == (exY * substitute)));
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
        result.setGf(truncatedSeries(result._gf, xExp, order));
        return result;
    }

    Dist Dist::CoefficientIterator::next() {
        ex nextTerm = truncatedSeries(_rest, _var, _n);
        _n++;
        _rest = (_rest - nextTerm).normal();
        ex nextCoeff = atOne(nextTerm, _var);
        return withSymbols(nextCoeff);
    }

    Dist Dist::CoefficientIterator::withSymbols(const ex &gf) const {
        Dist result;
        result._vars = _vars;
        result._params = _params;
        result.setGf(gf);
        return result;
    }

    Dist::CoefficientIterator Dist::coefficientIterator(const std::string &var) const {
        EXPECTS(isKnownAsVar(var), "can only expand series in known variables -- not in " + var);
        return CoefficientIterator(*this, _vars.at(var));
    }

    Dist Dist::substituteParam(const std::string &p, const std::string &val) const {
        EXPECTS(!isKnownAsVar(p), "symbol " + p + " is considered a variable");
        if (!isKnownAsParam(p)) {
            return *this;
        }
        Dist result(*this);
        ex pSym = result._params[p];
        ex parsed = result.parseWithKnownSymbs(val, false);
        result.registerAllUnknownSymbsAsParams(parsed);
        // todo error if parsed contains known variables
        result.setGf(result._gf.subs(pSym == parsed));
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
        auto ex_rhs = result.parseWithKnownSymbs(rhs, false);

        result.registerAllUnknownSymbsAsParams(ex_rhs);
        // todo error if parsed contained known variables
        // parsed contains just unknown symbs and known params
        // no error from now on
        result.setGf(result._gf * ex_rhs);
        return result;
    }

} // namespace prodigy
