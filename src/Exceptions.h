//
// Created by Tobias Winkler on 25.10.21.
//

#ifndef PRODIGY_EXCEPTIONS_H
#define PRODIGY_EXCEPTIONS_H

#include <stdexcept>

namespace prodigy {

    // the python bindings translate it into NotImplementedError
    class not_implemented_error : public std::logic_error {
    public:
        using std::logic_error::logic_error;
    };

    class bad_input_error : public std::invalid_argument {
    public:
        using std::invalid_argument::invalid_argument;
    };

}

#endif //PRODIGY_EXCEPTIONS_H
