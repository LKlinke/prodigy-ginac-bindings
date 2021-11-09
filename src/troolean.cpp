//
// Created by Tobias Winkler on 09.11.21.
//

#include "troolean.h"

namespace prodigy {

    troolean makeTroolean(bool b) {
        if(b) return troolean::YES;
        return troolean::NO;
    }

    std::ostream& operator<<(std::ostream& os, troolean tr) {
        if(tr == troolean::YES) {
            return os << "yes";
        } else if(tr == troolean::NO) {
            return os << "no";
        } else {
            return os << "unknown";
        }
    }
}
