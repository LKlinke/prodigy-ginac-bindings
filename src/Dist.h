//
// Created by Tobias Winkler on 22.10.21.
//

#ifndef PRODIGY_DIST_H
#define PRODIGY_DIST_H


#include <ginac/ginac.h>

#include <utility>
#include "Exceptions.h"
#include "macros.h"
#include "troolean.h"

namespace prodigy {

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

        friend class DistributionFactory;

    private:

        /* private member variables */

        ex _gf;

        /* static (global) part of class Dist */

        // global invariant:
        // _params and _vars are disjoint, they map a string "s" to a symbol symb with symb.get_name() = "s"
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

        /*
         * attempts to convert the given string to an expression
         * links all syntactic symbol names (e.g. "x") to the symbol objects registered in the tables _vars and _params
         * new symbols are created for unknown names, but they are not registered
         */
        static ex parseWithKnownSymbs(const std::string& s);
        static void registerAllUnknownSymbsAsVars(const ex& e);
        static void registerAllUnknownSymbsAsParams(const ex& e);
        static void registerAsVar(const symbol& p);
        static void registerAsParam(const symbol& p);
        static void deregister(const std::string& s);


    public:

        /* constructors */

        Dist() = default;

        /*
         * constructs a distribution with gf 'gfStr'
         * all symbols in 'gfStr' are registered as variables
         */
        explicit Dist(const std::string& gfStr);

        /*
         * constructs a distribution with gf 'gfStr' containing the symbols in 'paramsStr'
         * the symbols in 'paramsStr' are registered as parameters
         * all other symbols in 'gfStr' are registered as variables
         */
        Dist(const std::string& gfStr, const std::vector<std::string>& paramsStr);

        /*
         * constructs a distribution with gf 'gfStr' containing a symbol 'pStr'
         * the symbol 'pStr' is registered as a parameter
         * all other symbols in 'gfStr' are registered as variables
         */
        Dist(const std::string& gfStr, const std::string& pStr) : Dist(gfStr, std::vector<std::string>{pStr}) { /* intentionally left empty */ }



    private:
        explicit Dist(ex gf) : _gf(std::move(gf)) { /* intentionally left empty */ }

    public:

        /* public member functions */

        const ex& gf() const;

        /*
         * returns total probability mass, possibly non-numeric due to params
         */
        ex mass() const;

        /*
         * approximate size of underlying gf representation
         */
        int gfSize() const;

        /*
         * turns subdistribution into proper distribution by multiplying with 1/mass()
         */
        Dist normalize() const;

        /*
         * returns marginal distribution in given variable(s)
         */
        Dist marginal(const std::string& x) const;
        Dist marginal(const std::vector<std::string>& vars) const;
        Dist marginal(std::initializer_list<std::string> vars) const { return this->marginal(std::vector<std::string>{vars}); }

        /*
         * attempts to prove that *this is the zero distribution (mass 0)
         */
        troolean isZero() const;

        /*
         * attempts to prove that *this is a trivial distribution
         * (i.e. all vars trivial meaning that no state except (0,...,0) has positive mass)
         */
        troolean isTrivial() const;

        /*
         * attempts to prove that *this is univariate (at most one non-trivial var)
         * if *this is non-trivial and univariate, then optional res is the variable, otherwise res is undefined
         */
        troolean isUnivariate(std::string& res) const;
        troolean isUnivariate() const;

        /*
         * attempts to prove that x and y are independent variables in *this
         */
        troolean areIndependent(const std::string& x, const std::string& y) const;

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
         * update distribution according to assignment x := iid(D,y) where x and y are variables
         * D is a univariate proper distribution (no sub distribution)
         * the function iid(D,y) returns the sum of y independent D-distributed samples
         */
        Dist updateIid(const std::string& x, const Dist& d, const std::string& y) const;

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

        /*
         * replaces given parameter with given value or parametric expression
         */
        Dist substituteParam(const std::string& p, const std::string& val) const;

        /* operators */

        friend std::ostream& operator<<(std::ostream& os, const Dist& dist);

        friend bool operator==(const Dist& lhs, const Dist& rhs);
        friend bool operator!=(const Dist& lhs, const Dist& rhs) { return !(lhs == rhs); }

        /*
         * scale a distribution by a constant (possibly containing params)
         */
        friend Dist operator*(const Dist& lhs, const std::string& rhs);
        friend Dist operator*(const std::string& lhs, const Dist& rhs) { return rhs * lhs; }

        /*
         * convolution of distributions
         */
        friend Dist operator*(const Dist& lhs, const Dist& rhs) { return Dist {lhs._gf * rhs._gf}; }

        friend Dist operator+(const Dist& lhs, const Dist& rhs) { return Dist {lhs._gf + rhs._gf}; }
        friend Dist operator-(const Dist& lhs, const Dist& rhs) { return lhs + ("-1" * rhs); }

    };

} // namespace prodigy

#endif //PRODIGY_DIST_H