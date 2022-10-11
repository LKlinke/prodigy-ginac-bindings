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
     * represents a (possibly parametric) distribution over non-negative integer variable valuations using gfs
     * 'distribution' does not necessarily refer to probability distribution
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

        /**
         * returns a symbol name that is registered as neither a parameter nor a variable, nor contained in the exclude parameter
         * also registers the symbol as a variable
         */
        static std::string getFreshVariable();

        /**
         * evaluates the expression in a given state
         * does not read from or write to the symbol table, the only requirement is that valuations to all symbols in the expression are given
        */
        static int evaluate(const std::string &expr, const std::map<std::string, int> &valuations);

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
        Dist(const std::string& gfStr, std::initializer_list<std::string> paramsStr) : Dist(gfStr, std::vector<std::string>{paramsStr}) { }

        /*
         * constructs a distribution with gf 'gfStr' containing a symbol 'pStr'
         * the symbol 'pStr' is registered as a parameter
         * all other symbols in 'gfStr' are registered as variables
         */
        Dist(const std::string& gfStr, const std::string& pStr) : Dist(gfStr, std::vector<std::string>{pStr}) { }



    private:
        explicit Dist(ex gf) : _gf(std::move(gf)) { }

    public:

        /* public member functions */

        const ex& gf() const;

        /*
         * returns the specified raw moment (possibly non-numeric due to presence of parameters)
         * if *this is a distribution, then the raw moment is the expected value of the given monomial
         * this function does not modify the global state (//todo why?)
         */
        ex moment(const std::map<std::string, int>& monomial) const;

        /*
         * returns total probability mass, possibly non-numeric due to params
         * // todo implement in terms of moment
         */
        ex mass() const;

        /*
         * returns expected value, possibly non-numeric due to params
         * // todo implement in terms of moment
         */
        ex E(const std::string& var) const;

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
        Dist marginal(std::initializer_list<std::string> vars) const {
            return this->marginal(std::vector<std::string>{vars});
        }

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
         * attempts to prove that *this is finite-support, i.e., just finitely many coefficients are non-zero
         */
        troolean isFiniteSupport() const;

        /**
         * determines whether *this is a polynomial in the given variables
         * can be used as a heuristic for finite support
         */
        troolean isPolynomial(const std::set<std::string>& vars) const;

        typedef struct State {
            ex prob;
            std::map<std::string, int> vals;

            std::string toString() const {
                std::stringstream res;
                res << prob;
                for(auto item : vals) {
                    res << '*' << item.first << '^' << item.second;
                }
                return res.str();
            }
        } State;
        /**
         * returns a list of all terms if *this is a polynomial
         */
        std::vector<State> getTerms(const std::set<std::string>& vars) const;

        /**
         * creates a vector where at index i, the given variable is congruent
         * i modulo the given modulus
         */
        std::vector<Dist> arithmeticProgression(const std::string &var, const std::string &modulus) const;

        /**
         * returns all symbols that occur in this distribution (variables and parameters)
         */
        std::vector<symbol> getSymbols() const;

        /*
         * update distribution according to assignment x := e
         */
        Dist update(const std::string& x, const std::string& e) const;

        /**
         * performs the update updatedVar := assignVar
         * assignVar may be a variable or an integer (not a parameter)
         */
        Dist updateVar(const std::string &updatedVar, const std::string &assignVar) const;

        /**
         * applies the update updatedVar := firstSummand + secondSummand to *this
         */
        Dist updateSum(const std::string &updatedVar, const std::string &firstSummand, const std::string secondSummand) const;

        /**
         * applies the update updatedVar := firstFactor * secondFactor to *this
         */
        Dist updateProduct(const std::string &updatedVar, const std::string &firstFactor, const std::string secondFactor, 
            const std::set<std::string> variables, const bool finite) const;

        /**
         * applies the update updatedVar := subFrom - sub to *this
         * throws an exception if this might result in negative values
         */
        Dist updateSubtraction(const std::string updatedVar, const std::string subFrom, const std::string sub) const;

        /**
         * applies the update updatedVar := left % right to *this
         */
        Dist updateModulo(const std::string &updatedVar, const std::string &left, const std::string right, 
            const std::set<std::string> variables, const bool finite) const;

        /**
         * applies the update updatedVar := numerator / denominator to *this
         * throws an exception if this might result in a non-integer value
         */
        Dist updateDivision(const std::string updatedVar, const std::string numerator, const std::string denominator) const;
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
         * this iterator regards a k-dimensional distribution as 1-dimensional in the given variable with coefficients
         * that are k-1 dimensional, and it iterates of these coefficients
         *
         * in other words, the iterator ranges over the coefficients in the series expansion wrt. the given of the gf
         * underlying the given distribution
         */
        class CoefficientIterator {
        private:
            ex _rest;
            symbol _var;
            int _n;
        public:
            CoefficientIterator(Dist dist, symbol var) : _rest(std::move(dist._gf)), _var(std::move(var)), _n(1) {}
            Dist next();
            Dist rest() const { return Dist(_rest); }
        };

        /*
         * returns a coefficient iterator for the dimension specified by var
         * throws an exception if var is not registered as a variable
         */
        CoefficientIterator coefficientIterator(const std::string& var) const;

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
         * TODO we don't want this in the python anymore as the semantics regarding what is var/param is not consistent
         */
        friend Dist operator*(const Dist& lhs, const std::string& rhs);
        friend Dist operator*(const std::string& lhs, const Dist& rhs) { return rhs * lhs; }

        /*
         * convolution
         */
        friend Dist operator*(const Dist& lhs, const Dist& rhs) { return Dist {lhs._gf * rhs._gf}; }

        /*
         * point-wise sum/difference
         */
        friend Dist operator+(const Dist& lhs, const Dist& rhs) { return Dist {lhs._gf + rhs._gf}; }
        friend Dist operator-(const Dist& lhs, const Dist& rhs) { return lhs + ("-1" * rhs); }

    };

} // namespace prodigy

#endif //PRODIGY_DIST_H