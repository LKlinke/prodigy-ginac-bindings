# pygin python bindings demo

import pygin as pg

FALSE = pg.troolean(0)
TRUE = pg.troolean(1)
UNKNOWN = pg.troolean(2)


def check_if_zero(dist):
    if dist.is_zero() == FALSE:
        print("given distribution does not have zero mass")
    elif dist.is_zero() == TRUE:
        print("given distribution has zero mass")
    elif dist.is_zero() == UNKNOWN:
        print("it is unknown if given distribution has zero mass")


if __name__ == "__main__":
    geom_param = "p"
    var = "x"

    dist = pg.geometric(var, geom_param)
    print(f"expected value of a {geom_param}-geometric distribution is {dist.E(var)}")

    check_if_zero(dist)
    check_if_zero(pg.Dist("0"))
