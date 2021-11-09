//
// Created by Tobias Winkler on 01.11.21.
//

#ifndef PRODIGY_TROOLEAN_H
#define PRODIGY_TROOLEAN_H

#include <ostream>

namespace prodigy {

    enum class troolean {
        NO,
        YES,
        UNKNOWN
    };

    troolean makeTroolean(bool b);

    std::ostream& operator<<(std::ostream& os, troolean tr);

}

#endif //PRODIGY_TROOLEAN_H
