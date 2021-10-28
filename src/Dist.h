//
// Created by Tobias Winkler on 22.10.21.
//

#ifndef PRODIGY_DIST_H
#define PRODIGY_DIST_H


#include <ginac/ginac.h>

#include <utility>
#include "Exceptions.h"
#include "macros.h"


/*
 * represents a (possibly parametric) probability distribution over non-negative integer variable valuations using gfs
 *
 * Notes:
 * - Variables and parameters are managed globally.
 * - Variables that are non-zero with positive probability are called non-trivial. There may be at most a finite number
 *   of non-trivial variables.
 * - The distribution knows its exact set of parameters. Everything that is no parameter is considered a variable.
 *   Variables and parameters are strictly separated.
 */
class Dist {

    using ex = GiNaC::ex;
    using symbol = GiNaC::symbol;
    using str2sym = std::map<std::string, ex>;

private:

    /* private member variables */

    ex _gf;

    /* static (global) portion */

    // global invariant: the two maps are disjoint, they map a string "s" to a symbol symb with symb.get_name() = "s" jhd d
    static str2sym _params;
    static str2sym _vars;

public:

    static void printKnownVars() {
        std::for_each(_vars.begin(), _vars.end(), [](const auto& entry) {std::cout << entry.first << " ";} );
        std::cout << std::endl;
    }

    static void printKnownParams() {
        std::for_each(_params.begin(), _params.end(), [](const auto& entry) {std::cout << entry.first << " ";} );
        std::cout << std::endl;
    }

    static void resetSymbolCache() {
        _params.clear();
        _vars.clear();
    }
private:

    static str2sym paramsAndVars();

public:
    static bool isSymbKnown(const std::string& s);
    static bool isKnownAsVar(const std::string& s);
    static bool isKnownAsParam(const std::string& s);
private:
    static bool isSymbKnown(const symbol& s);
    static bool isKnownAsVar(const symbol& s);
    static bool isKnownAsParam(const symbol& s);

    static ex parseWithKnownSymbs(const std::string& s);
    static void registerAllUnknownSymbsAsVars(const ex& e);
    static void registerAllUnknownSymbsAsParams(const ex& e);
    static void registerAsVar(const symbol& p);
    static void registerAsParam(const symbol& p);


public:

    /* constructors */

    Dist() = default;
    Dist(const std::string& gfStr, const std::string& pStr);
    explicit Dist(const std::string& gfStr);
    private: explicit Dist(ex gf) : _gf(std::move(gf)) { /* intentionally left empty */ } public:

    /* public member functions */

    /*
     * returns total probability mass, possibly non-numeric due to params
     */
    ex mass() const;

    /*
     *
     */
    Dist marginal(const std::string& x) const;

    /*
     * returns expected value, possibly non-numeric due to params
     * does not modify global state
     */
    ex E(const std::string& var) const;

    /*
     * update distribution according to assignment x := e
     */
    Dist update(const std::string& x, const std::string& e) const;
private:
    Dist _update(const symbol& x, const ex& e) const;
    Dist _incrAssgn(const symbol& x, const ex& term) const;
public:

    /*
     * less: x < e
     * leq: x <= e
     * greater: x > e
     * geq: x >= e
     * eq: x = e
     * neq: x != e
     */
    Dist filterLess(const std::string& x, const std::string& e) const;
    // todo make more efficient implementations
    Dist filterLeq(const std::string& x, const std::string& e) const { return filterLess(x, e + "+1"); }
    Dist filterGreater(const std::string& x, const std::string& e) const { return *this - filterLeq(x, e); }
    Dist filterGeq(const std::string& x, const std::string& e) const { return *this - filterLess(x, e); }
    Dist filterEq(const std::string& x, const std::string& e) const { return *this - filterLess(x, e) - filterGreater(x, e); }
    Dist filterNeq(const std::string& x, const std::string& e) const { return *this - filterEq(x, e); }


    /* operators */

    friend std::ostream& operator<<(std::ostream& os, const Dist& dist);;
    friend bool operator==(const Dist& lhs, const Dist& rhs);

    friend Dist operator*(const Dist& lhs, const std::string& rhs);
    friend Dist operator*(const std::string& lhs, const Dist& rhs) { return rhs * lhs; }
    friend Dist operator*(const Dist& lhs, const Dist& rhs) { return Dist {lhs._gf * rhs._gf}; }

    friend Dist operator+(const Dist& lhs, const Dist& rhs) { return Dist {lhs._gf + rhs._gf}; }
    friend Dist operator-(const Dist& lhs, const Dist& rhs) { return lhs + ("-1" * rhs); }

};

#endif //PRODIGY_DIST_H