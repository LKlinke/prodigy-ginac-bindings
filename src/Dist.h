//
// Created by Tobias Winkler on 22.10.21.
//

#ifndef PRODIGY_DIST_H
#define PRODIGY_DIST_H


#include <ginac/ginac.h>
#include "Exceptions.h"

using namespace GiNaC;

/*
 * represents a (possibly parametric) probability distribution over non-negative integer variable valuations using gfs
 *
 * Notes:
 * - Variables that are non-zero with positive probability are called non-trivial. There may be at most a finite number
 *   of non-trivial variables.
 * - The distribution knows its exact set of parameters. Everything that is no parameter is considered a variable.
 *   Variables and parameter are strictly separated.
 */
class Dist {

private:
    ex _gf;
    std::vector<symbol> _params;

public:

    /* Constructors */

    Dist() = default;
//
//    /*
//     * copy constructor
//     */
//    Dist(const Dist& other) {
//        _gf = other._gf;
//        params = other.params;
//    };

    template<typename Container>
    Dist(const ex &gf, const Container &params) {
        this->_gf = gf;
        this->_params = std::vector<symbol>(std::begin(params), std::end(params));
    }

    /* public member functions */

    /*
     * returns the underlying raw gf
     */
    const ex& gf() const {
        return _gf;
    };

    /*
     * returns total probability mass
     */
    ex mass() const;

    /*
     * returns expected value
     * note that this is possibly non-numeric due to params
     */
    ex E(const symbol& var) const;
    ex E(const std::string& var) const;


    /*
     * all non-trivial program variables
     */
    std::vector<symbol> getVariables() const;

    /*
     * returns variable with the given name, throws if not present
     */
    symbol getVariableByName(const std::string& v) const;

    /*
     * true iff this distribution is non-trivial in the given variable
     */
    bool hasNonTrivialVariable(const symbol&) const;
    bool hasNonTrivialVariable(const std::string&) const;

    /*
     * returns parameter with the given name, throws if not present
     */
    symbol getParameterByName(const std::string& v) const;

    /*
     * true iff this distribution considers given symbol as a parameter
     */
    bool hasParameter(const std::string& p) const;
    bool hasParameter(const symbol& p) const;

    /*
     * update distribution according to assignment v := e
     */
    Dist update(const symbol& v, const ex& e) const;
    Dist update(const std::string& v, const std::string& e) const;


    /* operators */

    friend std::ostream& operator<<(std::ostream& os, const Dist& dist) {
        os << dist._gf;
        return os;
    };

    bool operator==(const Dist& rhs) const;

    Dist operator*(const numeric& rhs) const { return Dist(_gf * rhs, _params); }
    friend Dist operator*(const numeric& lhs, const Dist& rhs) { return rhs * lhs; }

    Dist operator*(const std::string& rhs) const;
    friend Dist operator*(const std::string& lhs, const Dist& rhs) { return rhs * lhs; }

private:
    Dist incrAssgn(const symbol& v, const ex& term) const;
    symtab makeSymtab() const;



};

#endif //PRODIGY_DIST_H
