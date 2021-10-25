#include <pybind11/pybind11.h>
#include <sstream>

#include "DistributionFactory.h"

#define STRINGIFY(x) #x
#define MACRO_STRINGIFY(x) STRINGIFY(x)

// converts GiNaC::ex to string
std::string toString(const ex& e) {
    std::stringstream s;
    s << e;
    return s.str();
}

ParametricDist geometric(const std::string& x, const std::string& p) {
    return DistributionFactory::geometric(symbol(x), symbol(p));
}

//ParametricDist geometric(const std::string& x, double p) {
//    return DistributionFactory::geometric(symbol(x), numeric(p));
//}

namespace py = pybind11;

PYBIND11_MODULE(prodigy, m) {
    m.doc() = "python bindings of prodigy - a library for manipulating distributions through generating functions";

    py::class_<ParametricDist>(m, "ParametricDist")
            .def(py::init())
            .def("E", [](const ParametricDist& dist, const std::string& x) {
                return toString(dist.E(x));
            })
            .def("hasParameter", [](const ParametricDist& dist, const std::string& p) {
                return dist.hasParameter(p);
            })
            .def("__mul__", [](const ParametricDist& lhs, const std::string& rhs) {
                return lhs * rhs;
            }, py::is_operator())
            .def("__mul__", [](const ParametricDist& lhs, float rhs) {
                return lhs * numeric(rhs);
            }, py::is_operator());

    m.def("geometric", &geometric, "creates a geometric distribution");




    #ifdef VERSION_INFO
    m.attr("__version__") = MACRO_STRINGIFY(VERSION_INFO);
    #else
    m.attr("__version__") = "dev";
    #endif
}