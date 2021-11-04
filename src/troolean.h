//
// Created by Tobias Winkler on 01.11.21.
//

#ifndef PRODIGY_TROOLEAN_H
#define PRODIGY_TROOLEAN_H

namespace prodigy {

    enum class troolean {
        NO,
        YES,
        UNKNOWN
    };

    constexpr troolean makeTroolean(bool b) {
        if(b) return troolean::YES;
        return troolean::NO;
    }

}

#endif //PRODIGY_TROOLEAN_H
