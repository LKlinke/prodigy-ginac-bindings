#include <pybind11/pybind11.h>
#include <pybind11/operators.h>
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

PYBIND11_MODULE(prodigy, m) {
    m.doc() = "python bindings of prodigy - a library for manipulating distributions through generating functions";

    py::class_<Dist>(m, "Dist")
            .def(py::init())
            .def(py::init<const string&>())
            .def(py::init<const string&, const string&>())
            .def("E", [](const Dist& dist, const string& x) { return toString(dist.E(x)); })
            .def("update", &Dist::update)
            .def("mass", [](const Dist& dist) { return toString(dist.mass()); })
            .def("normalize", &Dist::normalize)
            .def("isZero", &Dist::isZero)
            .def(py::self == py::self)
            .def(py::self != py::self)
            .def(string() * py::self)
            .def(py::self * string())
            .def(py::self * py::self)
            .def(py::self + py::self)
            .def(py::self - py::self)
            .def("__repr__", [](const Dist& dist) { return toString(dist); });

    m.def("geometric", &DistributionFactory::geometric, "creates a geometric distribution");




    #ifdef VERSION_INFO
    m.attr("__version__") = MACRO_STRINGIFY(VERSION_INFO);
    #else
    m.attr("__version__") = "dev";
    #endif
}
