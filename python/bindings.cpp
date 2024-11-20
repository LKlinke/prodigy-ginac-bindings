#include <pybind11/pybind11.h>
#include <pybind11/operators.h>
#include <pybind11/stl.h> // for STL containers
#include <sstream>

#include "DistributionFactory.h"

#define STRINGIFY(x) #x
#define MACRO_STRINGIFY(x) STRINGIFY(x)

using std::string;
using GiNaC::ex;
using namespace prodigy;

// converts GiNaC::ex to string
string toString(const ex& e) {
    std::stringstream s;
    s << e;
    return s.str();
}

string toString(const Dist& dist){
    std::stringstream s;
    s << dist;
    return s.str();
}

namespace py = pybind11;

PYBIND11_MODULE(pygin, m) {
    m.doc() = "python bindings of pygin - a library for manipulating distributions through generating functions";

    py::class_<Dist>(m, "Dist")
            .def(py::init())
            .def(py::init<const string&>())
            .def(py::init<const string&, const string&>())
            .def(py::init<const string&, const py::args>)
            .def("E", [](const Dist& dist, const string& x) { return toString(dist.E(x)); })
            .def("mass", [](const Dist& dist, const std::optional<std::set<std::string>> vars) { 
                return toString(dist.mass(vars)); }, py::arg("vars") = std::optional<std::set<std::string>>())
            .def("normalize", &Dist::normalize)
            .def("marginal", [](const Dist& dist, const string& x) {return dist.marginal(x);})
            .def("marginal", [](const Dist& dist, const std::vector<string>& vars) {return dist.marginal(vars);})
            .def("is_zero", &Dist::isZero)
            .def("is_trivial", &Dist::isTrivial)
            .def("updateIid", &Dist::updateIid)
            .def("filterLess", &Dist::filterLess)
            .def("filterEq", &Dist::filterEq)
            .def("filterNeq", &Dist::filterNeq)
            .def("filterLeq", &Dist::filterLeq)
            .def("filterGreater", &Dist::filterGreater)
            .def("filterGeq", &Dist::filterGeq)
            .def("coefficient_iterator", &Dist::coefficientIterator)
            .def(py::self == py::self)
            .def(py::self != py::self)
            //.def(string() * py::self)
            //.def(py::self * string())
            .def(py::self * py::self)
            .def(py::self + py::self)
            .def(py::self - py::self)
            .def("__repr__", [](const Dist& dist) { return toString(dist); })
            .def("is_polynomial", &Dist::isPolynomial)
            .def("get_terms", [](const Dist& dist, const std::set<string>& vars) {
                std::vector<std::pair<std::string, std::map<std::string, int>>> res;
                for(auto state : dist.getTerms(vars)) {
                    auto prob = toString(state.prob);
                    res.push_back(std::pair<std::string, std::map<std::string, int>>(prob, state.vals));
                }
                return res;
            })
            .def("arithmetic_progression", &Dist::arithmeticProgression)
            .def("get_symbols", [](const Dist &dist) {
                std::vector<string> res;
                for(auto sym : dist.getSymbols()) {
                    res.push_back(toString(sym));
                }
                return res;
            })
            .def("update_var", &Dist::updateVar)
            .def("update_sum", &Dist::updateSum)
            .def("update_product", &Dist::updateProduct)
            .def("update_subtraction", &Dist::updateSubtraction)
            .def("update_modulo", &Dist::updateModulo)
            .def("update_division", &Dist::updateDivision)
            .def("update_power", &Dist::updatePower)
            .def("update_root", &Dist::updateRoot)
            .def("approximate_unilaterally", &Dist::approximateUnilaterally)
            .def("set_assumptions", &Dist::setAssumptions)
            .def("get_assumptions", &Dist::getAssumptions)
            .def("assumptions_hold", &Dist::assumptionsHold);

    // global functions in python
    m.def("geometric", &DistributionFactory::geometric, "creates a geometric distribution")
    .def("uniform", &DistributionFactory::uniform, "creates a uniform distribution")
    .def("dirac", &DistributionFactory::dirac, "creates a Dirac distribution (point mass)")
    .def("catalan", &DistributionFactory::catalan, "creates a Catalan distribution (runtime distribution of 1D random walk)")
    .def("is_known_as_var", [](const string& x) { return Dist::isKnownAsVar(x); })
    .def("is_known_as_param", [](const string& x) { return Dist::isKnownAsParam(x); })
    .def("reset_symbol_cache", &Dist::resetSymbolCache)
    .def("find_symbols", &Dist::find_symbols)
    .def("get_fresh_variable", &Dist::getFreshVariable)
    .def("evaluate", &Dist::evaluate)
    .def("nth_root", &Dist::nthRoot);

    py::enum_<troolean>(m, "troolean")
            .value("false", troolean::FALSE)
            .value("true", troolean::TRUE)
            .value("unknown", troolean::UNKNOWN);

    py::class_<Dist::CoefficientIterator>(m, "CoefficientIterator")
            .def("next", &Dist::CoefficientIterator::next)
            .def("rest", &Dist::CoefficientIterator::rest);


    // explicitly define that we use shared pointers to reference Assumptions instead of the default unique pointers
    // this means that if an instance goes out of scope in Python, it will keep existing in C++
    py::class_<Assumption, std::shared_ptr<Assumption>>(m, "Assumption")
            .def("implies", &Assumption::implies);

    // explicitly define that Divisible inherits from Assumption
    py::class_<Divisible, std::shared_ptr<Divisible>, Assumption>(m, "Divisible")
            .def("implies", &Divisible::implies)
            .def("get_lhs", &Divisible::getLhs)
            .def("get_rhs", &Divisible::getRhs)
            .def(py::init<const std::string&, const int&>());


    #ifdef VERSION_INFO
    m.attr("__version__") = MACRO_STRINGIFY(VERSION_INFO);
    #else
    m.attr("__version__") = "dev";
    #endif
}
