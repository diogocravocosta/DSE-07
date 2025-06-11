import pytest
import numpy.testing as npt
import h2ermes_tools.fatigue_calculation_tool as fct


def test_mechanical_stress():
    stress = fct.mechanical_stress(1000, 5, 0.01, 10)
    test_stress = 3232.2
    npt.assert_almost_equal(test_stress, stress,decimal=1)

def test_t_crit_mechanical():
    t_crit = fct.t_crit_mechanical(1000,5,10e6,10)
    t_act = 0.0000323
    npt.assert_almost_equal(t_act, t_crit,decimal=1)

def test_thermal_stress():
    stress = fct.thermal_stress(200, 1e8, 1e-6,5,0.01, 10)
    test = 20308.53
    npt.assert_almost_equal(stress, test,decimal=1)

