import pytest
import numpy.testing as npt
import src.fatigue_calculation_tool as fct


def test_mechanical_stress():
    stress = fct.mechanical_stress(1000, 5, 0.01, 0.1)
    test_stress = 3387.383
    npt.assert_almost_equal(test_stress, stress)