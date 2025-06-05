"""Example tests for example file"""

import pytest
import numpy.testing as npt

from h2ermes_tools.example import dot_product
from h2ermes_tools.example import Cylinder

def test_dot_product1():
    vec1 = [1, 2, 3]
    vec2 = [0.5, 0.2, 1/3]

    expected_output = 1.9

    result = dot_product(vec1, vec2)

    npt.assert_almost_equal(result, expected_output)

def test_dot_product2():
    vec1 = [1, 2, 3]
    vec2 = [0.5, 0.2, 1/3, 5]

    with pytest.raises(ValueError) as exception:
        dot_product(vec1, vec2)
    assert exception.value.args[0] == "vector1 and vector2 must be the same length"

def make_test_cylinder():
    expected_cylinder = Cylinder.__new__(Cylinder)
    expected_cylinder.radius = 0.5
    expected_cylinder.height = 0.2

    return expected_cylinder

def test_cylinder_init():
    expected_cylinder = make_test_cylinder()

    test_cylinder = Cylinder(0.5, 0.2)

    assert expected_cylinder.radius == test_cylinder.radius
    assert expected_cylinder.height == test_cylinder.height

def test_cylinder_area():
    cylinder = make_test_cylinder()
    expected_area = 0.7853981 # pi * radius**2

    cylinder.calculate_base_area()
    npt.assert_almost_equal(cylinder.area, expected_area)

def test_cylinder_volume():
    cylinder = make_test_cylinder()
    cylinder.area = 0.7853981
    expected_volume = 0.157079

    cylinder.calculate_volume()
    npt.assert_almost_equal(cylinder.volume, expected_volume, decimal=3)