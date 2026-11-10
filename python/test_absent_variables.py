# A variable that does not occur in the gf is 0 in every state; updates read it as 0, not as a new symbol (issue #22).
import pytest

import prodigy_ginac as pg

GF = "(1+y+y^2)/3"


@pytest.mark.parametrize("update, expected", [
    (lambda d: d.update_var("y", "x"), "1"),
    (lambda d: d.update_var("x", "z"), GF),
    (lambda d: d.update_sum("y", "z", "w"), "1"),
    (lambda d: d.update_sum("y", "z", "1"), "y"),
    (lambda d: d.update_sum("y", "y", "z"), GF),
    (lambda d: d.update_product("y", "z", "y", {"y", "z"}, True, None), "1"),
    (lambda d: d.update_division("y", "z", "2", None), "1"),
    (lambda d: d.update_modulo("y", "z", "2", {"y", "z"}, True, None), "1"),
    (lambda d: d.update_power("y", "2", "z", None), "y"),
])
def test_absent_operand_is_zero(update, expected):
    assert update(pg.Dist(GF)) == pg.Dist(expected)


def test_infinite_support():
    assert pg.Dist("1/(2-y)").update_var("y", "z") == pg.Dist("1")


def test_variable_that_occurs_is_copied():
    assert pg.Dist("z/(2-y)").update_var("y", "z") == pg.Dist("z*y")
