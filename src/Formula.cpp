#include "Formula.h"


int Formula::bla(int arg1) {
    return arg1 * 2;
}

ex Formula::giveIt2Me() {
    symbol x("x");
    return ex(x + x);
}
