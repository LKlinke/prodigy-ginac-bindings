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
            .def(py::init<const string&, const std::vector<string>&>())
            .def("E", [](const Dist& dist, const string& x) { return toString(dist.E(x)); })
            .def("mass", [](const Dist& dist) { return toString(dist.mass()); })
            .def("normalize", &Dist::normalize)
            .def("marginal", [](const Dist& dist, const string& x) {return dist.marginal(x);})
            .def("marginal", [](const Dist& dist, const std::vector<string>& vars) {return dist.marginal(vars);})
            .def("is_zero", &Dist::isZero)
            .def("is_trivial", &Dist::isTrivial)
            .def("update", &Dist::update)
            .def("updateIid", &Dist::updateIid)
            .def("filterLess", &Dist::filterLess)
            .def("filterEq", &Dist::filterEq)
            .def("filterNeq", &Dist::filterNeq)
            .def("filterLeq", &Dist::filterLeq)
            .def("filterGreater", &Dist::filterGreater)
            .def("filterGeq", &Dist::filterGeq)
            .def(py::self == py::self)
            .def(py::self != py::self)
            //.def(string() * py::self)
            //.def(py::self * string())
            .def(py::self * py::self)
            .def(py::self + py::self)
            .def(py::self - py::self)
            .def("__repr__", [](const Dist& dist) { return toString(dist); });

    // global functions in python
    m.def("geometric", &DistributionFactory::geometric, "creates a geometric distribution")
    .def("uniform", &DistributionFactory::uniform, "creates a uniform distribution")
    .def("is_known_as_var", [](const string& x) { return Dist::isKnownAsVar(x); })
    .def("is_known_as_param", [](const string& x) { return Dist::isKnownAsParam(x); })
    .def("reset_symbol_cache", &Dist::resetSymbolCache);

    py::enum_<troolean>(m, "troolean")
            .value("false", troolean::FALSE)
            .value("true", troolean::TRUE)
            .value("unknown", troolean::UNKNOWN);


    #ifdef VERSION_INFO
    m.attr("__version__") = MACRO_STRINGIFY(VERSION_INFO);
    #else
    m.attr("__version__") = "dev";
    #endif
}
