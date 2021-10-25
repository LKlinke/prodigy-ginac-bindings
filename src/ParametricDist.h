//
// Created by Tobias Winkler on 22.10.21.
//

#ifndef PRODIGY_PARAMETRICDIST_H
#define PRODIGY_PARAMETRICDIST_H


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
class ParametricDist {

public:

    /* Constructors */

    ParametricDist() = default;

    /*
     * copy constructor
     */
    ParametricDist(const ParametricDist& other) {
        _gf = other._gf;
        params = other.params;
    };

    template<typename Container>
    ParametricDist(const ex &gf, const Container &params) {
        this->_gf = gf;
        this->params = std::vector<symbol>(std::begin(params), std::end(params));
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

    /* operators */

    friend std::ostream& operator<<(std::ostream& os, const ParametricDist& dist) {
        os << dist._gf;
        return os;
    };

    ParametricDist operator*(const numeric& rhs) const { return ParametricDist(_gf * rhs, params); }
    friend ParametricDist operator*(const numeric& lhs, const ParametricDist& rhs) { return rhs * lhs; }

    ParametricDist operator*(const std::string& rhs) const;
    friend ParametricDist operator*(const std::string& lhs, const ParametricDist& rhs) { return rhs * lhs; }


private:
    ex _gf;
    std::vector<symbol> params;

};

#endif //PRODIGY_PARAMETRICDIST_H
