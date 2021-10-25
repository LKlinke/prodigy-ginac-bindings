//
// Created by Tobias Winkler on 22.10.21.
//

#include "DistributionFactory.h"

using namespace GiNaC;

ParametricDist DistributionFactory::geometric(const symbol& x, const ex &p) {

    std::vector<symbol> params;
    // Analyse arguments
    if(p.info(info_flags::rational)) {
        if(!(0 <= p && p < 1)) {
            throw std::invalid_argument("rational parameter of geom distribution must be in [0,1)");
        }
    }
    else if(p.info(info_flags::symbol)) {
        // todo test properly if x and p are same symbol
        const auto& p_as_symbol = ex_to<symbol>(p);
        if(x.is_equal(p_as_symbol)) {
            throw std::invalid_argument("parameter and variable must be different");
        }
        params.push_back(p_as_symbol);
    }
    else {
        // p is neither rational nor a symbol - error
        throw std::invalid_argument("parameter of geom distribution must be a rational number or a symbol");
    }

    ex gf = (1-p) / (1 - p*x);
    return ParametricDist(gf, params);

}
