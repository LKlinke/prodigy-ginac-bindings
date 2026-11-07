//
// Created by Tobias Winkler on 22.10.21.
//

#ifndef PRODIGY_DIST_H
#define PRODIGY_DIST_H


#include <ginac/ginac.h>

#include <utility>
#include <sstream>
#include "Exceptions.h"
#include "macros.h"
#include "troolean.h"
#include <optional>
#include "Assumptions.h"


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
        using possymbol = GiNaC::possymbol;
        using str2sym = std::map<std::string, ex>;

        friend class DistributionFactory;

    private:

        /* private member variables */

        ex _gf;

        std::vector<std::shared_ptr<Assumption>> assumptions;
        
        str2sym _params;
        str2sym _vars;

    public:

        /* static parts */

        static std::vector<std::string> find_symbols(const std::string &exp);

        /*
         * parses 's', linking symbol names to the symbols in 'table'
         * unknown names are added to 'table', as possymbols if 'unknownAreVariables', otherwise as plain symbols
         *
         * Variables are possymbols, i.e. GiNaC regards them as positive reals; parameters are plain symbols, since
         * they may be negative. This is sound for variables: a gf is a power series in its variables, and GiNaC only
         * applies simplifications that hold for all positive values of the variables, eagerly, e.g.
         * (x^2)^(9/2) -> x^9 or sqrt(x^2) -> x. If both sides of such a simplification are analytic at 0, they agree
         * on an interval (0, r) and therefore, by the identity theorem, have the same power series. Later
         * substitutions such as x -> -x or x -> exp(2*Pi*I/k)*x then apply to the simplified gf.
         */
        static ex parse(const std::string& s, str2sym& table, bool unknownAreVariables = true);

        /*
         * rewrites every decimal literal in 's' into an exact rational, e.g. "0.1" -> "(1/10)", "2.5e-3" -> "(1/400)",
         * "1e20" -> "(100000000000000000000)"; integer literals, identifiers (e.g. "x1", "e", "exp") and comments
         * are left unchanged
         * GiNaC's parser reads decimal literals as floating-point numbers with about 17 digits, so equal-looking
         * but different numbers would compare equal; every string that is handed to GiNaC's parser goes through
         * this function first
         * throws std::invalid_argument for malformed literals (e.g. "1.2.3", "1e", ".") and decimal exponents of
         * absolute value above 100000
         */
        static std::string exactDecimals(const std::string& s);

        /*
         * true iff 'e' contains an inexact (floating-point) number
         */
        static bool hasInexactNumbers(const ex& e);

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
        Dist(const std::string& gfStr, std::initializer_list<std::string> paramsStr)
            : Dist(gfStr, std::vector<std::string>{paramsStr}) {}

        /*
         * constructs a distribution with gf 'gfStr' containing a symbol 'pStr'
         * the symbol 'pStr' is registered as a parameter
         * all other symbols in 'gfStr' are registered as variables
         */
        Dist(const std::string& gfStr, const std::string& pStr)
            : Dist(gfStr, std::vector<std::string>{pStr}) {}
        

        void printKnownVars() const {
            std::for_each(_vars.begin(), _vars.end(), [](const auto& entry) {std::cout << entry.first << " ";} );
            std::cout << std::endl;
        }

        void printKnownParams() const {
            std::for_each(_params.begin(), _params.end(),[](const auto& entry) {std::cout <<  entry.first << " ";} );
            std::cout << std::endl;
        }

        void resetSymbolCache() {
            _params.clear();
            _vars.clear();
        }
        
        const bool isSymbKnown(const std::string& s) const;
        const bool isKnownAsVar(const std::string& s) const;
        const bool isKnownAsParam(const std::string& s) const;
        
        /**
         * returns a symbol name that is registered as neither a parameter nor a variable
         * also registers the symbol as a variable
         */
        std::string getFreshVariable();
        
        /**
         * evaluates the expression in a given state
         * does not read from or write to the symbol table, the only requirement is that valuations to all symbols in the expression are given
        */
        static int evaluate(const std::string &expr, const std::map<std::string, int> &valuations);



        static uint64_t nthRoot(const uint64_t number, const uint64_t root);
        
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

    private:

        str2sym paramsAndVars() const;
        std::vector<State> getTerms(const GiNaC::lst& vars) const;

        const bool isSymbKnown(const symbol& s) const;
        const bool isKnownAsVar(const symbol& s) const;
        const bool isKnownAsParam(const symbol& s) const;

        /*
         * attempts to convert the given string to an expression
         * links all syntactic symbol names (e.g. "x") to the symbol objects registered in the tables _vars and _params
         * new symbols are created for unknown names, but they are not registered
         * the new symbols are variables (possymbols) if 'unknownAreVariables', otherwise parameters (plain symbols)
         */
        ex parseWithKnownSymbs(const std::string& s, bool unknownAreVariables = true) const;

        void registerAllUnknownSymbsAsVars(const ex& e);
        void registerAllUnknownSymbsAsParams(const ex& e);
        void registerAsVar(const symbol& p);
        void registerAsParam(const symbol& p);
        void deregister(const std::string& s);

        /*
         * returns e with every symbol replaced by the registered symbol of the same name
         * unknown symbols are registered as variables
         */
        ex linkSymbs(const ex& e);

        /*
         * sets the gf to gf, linked to the symbol table (see linkSymbs)
         * every assignment to _gf goes through this function, so that the gf contains only registered symbols
         * and variables or parameters of the same name are the same GiNaC symbol
         */
        void setGf(const ex& gf);

        /*
         * returns a distribution with gf 'gf', the symbol table of *this and no assumptions
         */
        Dist withGf(const ex& gf) const;

    public:

        /* public member functions */

        const ex& gf() const;

        const std::vector<std::shared_ptr<Assumption>> getAssumptions() const;

        const Dist setAssumptions(const std::vector<std::shared_ptr<Assumption>>& assumptions) const;

        const bool assumptionsHold(const std::vector<std::shared_ptr<Assumption>>& assumptions) const;

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
        ex mass(const std::optional<std::set<std::string>> vars = std::optional<std::set<std::string>>()) const;

        /*
         * returns the expected value of a polynomial expression in the variables (e.g. "x + y", "2*x - y/3 + 1",
         * "x*y^2"), possibly non-numeric due to params; the coefficients may contain params
         * throws not_implemented_error for expressions that are not polynomial in the variables
         * and std::invalid_argument if the expression is a single param or contains a symbol that is neither a
         * variable nor a param
         */
        ex E(const std::string& expr) const;

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

        /**
         * returns a list of all terms if *this is a polynomial in the provided variables
         */
        std::vector<State> getTerms(const std::set<std::string>& vars) const;

        /**
         * returns a list of all terms if *this is a polynomial
         */
        std::vector<State> getTerms() const;

        /**
         * creates a vector where at index i, the given variable is congruent
         * i modulo the given modulus
         * if the gf is rational in var, the parts are computed exactly as rational functions without roots of unity
         * (via a resultant, see rationalMultisection in Dist.cpp); otherwise var is substituted by exp(2*Pi*I*m/modulus)*var
         */
        std::vector<Dist> arithmeticProgression(const std::string &var, const std::string &modulus) const;

        /**
         * returns all symbols that occur in this distribution (variables and parameters)
         * as ex: copying a possymbol into a symbol would slice it into a different symbol
         */
        std::vector<ex> getSymbols() const;

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
        Dist updateProduct(const std::string &updatedVar, const std::string &firstFactor, const std::string &secondFactor, 
            const std::set<std::string> &variables, const bool &finite, const std::optional<std::string> &approx) const;

        /**
         * applies the update updatedVar := subFrom - sub to *this
         * throws an exception if this might result in negative values
         */
        Dist updateSubtraction(const std::string updatedVar, const std::string subFrom, const std::string sub) const;

        /**
         * applies the update updatedVar := left % right to *this
         */
        Dist updateModulo(const std::string &updatedVar, const std::string &left, const std::string &right, 
            const std::set<std::string> &variables, const std::optional<bool> &finite, const std::optional<std::string> &approx) const;

        /**
         * applies the update updatedVar := numerator / denominator to *this
         * throws an exception if this might result in a non-integer value
         */
        Dist updateDivision(const std::string &updatedVar, const std::string &numerator, const std::string &denominator, const std::optional<std::string> &approx) const;

        /**
         * applies the update updatedVar := base^exp to *this
        */
        Dist updatePower(const std::string &updatedVar, const std::string &base, const std::string &exp, const std::optional<std::string> &approx) const;

        /**
         * applies the update updatedVar := radicand^(1/index) to *this
        */
        Dist updateRoot(const std::string &updatedVar, const std::string &radicand, const std::string &index, const std::optional<std::string> &approx) const;

        /**
         * approximates the distribution in one variable via its series expansion
        */
        Dist approximateUnilaterally(const std::string &variable, const std::string &prob_mass) const;

    public:
        
        
        /**
         * this iterator regards a k-dimensional distribution as 1-dimensional in the given variable with coefficients
         * that are k-1 dimensional, and it iterates of these coefficients
         *
         * in other words, the iterator ranges over the coefficients in the series expansion wrt. the given of the gf
         * underlying the given distribution
         */
        class CoefficientIterator {
        private:
            str2sym _params;
            str2sym _vars;
            ex _rest;
            ex _var;
            int _n;
        public:
            CoefficientIterator(Dist dist, ex var) : _params(std::move(dist._params)), _vars(std::move(dist._vars)), _rest(std::move(dist._gf)), _var(std::move(var)), _n(1) {}
            Dist next();
            Dist rest() const { return withSymbols(_rest); }
        private:
            Dist withSymbols(const ex& gf) const;
        };
        
        
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
        friend Dist operator*(const Dist& lhs, const Dist& rhs) {
            Dist result(lhs);
            for (const auto& var : rhs._vars){
                EXPECTS(!result.isKnownAsParam(var.first), "variable " + var.first + " is known as parameter.");
                if(!result.isKnownAsVar(var.first)){
                    result.registerAsVar(GiNaC::ex_to<symbol>(var.second));
                }
            }
            for (const auto& param: rhs._params){
                EXPECTS(!result.isKnownAsVar(param.first), "parameter " + param.first + " is known as variable.");
                if(!result.isKnownAsParam(param.first)){
                    result.registerAsParam(GiNaC::ex_to<symbol>(param.second));
                }
            }

            // links the symbols of rhs to the symbols of the same name in result
            result.setGf(result._gf * rhs._gf);
            return result; 
        }

        /*
         * point-wise sum/difference
         */
        friend Dist operator+(const Dist& lhs, const Dist& rhs) {
            Dist result(lhs);
            for (const auto& var : rhs._vars){
                EXPECTS(!result.isKnownAsParam(var.first), "variable is known as parameter.");
                if(!result.isKnownAsVar(var.first)){
                    result.registerAsVar(GiNaC::ex_to<symbol>(var.second));
                }
            }
            for (const auto& param: rhs._params){
                EXPECTS(!result.isKnownAsVar(param.first), "parameter is known as variable.");
                if(!result.isKnownAsParam(param.first)){
                    result.registerAsParam(GiNaC::ex_to<symbol>(param.second));
                }
            }

            // links the symbols of rhs to the symbols of the same name in result
            result.setGf(result._gf + rhs._gf);
            return result; 
        }

        friend Dist operator-(const Dist& lhs, const Dist& rhs) { return lhs + ("(-1)" * rhs); }

    };

} // namespace prodigy

#endif //PRODIGY_DIST_H
