//
// Created by Tobias Winkler on 09.11.21.
//

#include "troolean.h"

namespace prodigy {

    troolean makeTroolean(bool b) {
        if(b) return troolean::TRUE;
        return troolean::FALSE;
    }

    std::ostream& operator<<(std::ostream& os, troolean tr) {
        if(tr == troolean::TRUE) {
            return os << "yes";
        } else if(tr == troolean::FALSE) {
            return os << "no";
        } else {
            return os << "unknown";
        }
    }
}
