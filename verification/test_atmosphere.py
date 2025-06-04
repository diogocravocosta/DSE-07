import pytest
import numpy.testing as npt

import src.data.constants as ct
from h2ermes_tools import atmosphere as atm

class TestAtmosphere:
    height = 0

    def make_test_atmosphere(self):
        atmosphere = atm.Atmosphere.__new__(atm.Atmosphere)
        atmosphere.g_0 = ct.g_0  # m/s^2
        atmosphere.Re = ct.radius  # m
        atmosphere.R_star = ct.R_star
        atmosphere.M = ct.molecular_mass
        atmosphere.R = atmosphere.R_star / atmosphere.M
        atmosphere.gamma = ct.gamma
        atmosphere.rho_0 = ct.density_sea_level
        atmosphere.h = self.height * 1000

        return atmosphere

    def test_atmosphere_init(self):
        expected_atmosphere = self.make_test_atmosphere()

        test_atmosphere = atm.Atmosphere(self.height)

        assert expected_atmosphere.g_0 == test_atmosphere.g_0
        assert expected_atmosphere.Re == test_atmosphere.Re
        assert expected_atmosphere.R_star == test_atmosphere.R_star
        assert expected_atmosphere.M == test_atmosphere.M
        assert expected_atmosphere.R == test_atmosphere.R
        assert expected_atmosphere.gamma == test_atmosphere.gamma
        assert expected_atmosphere.rho_0 == test_atmosphere.rho_0
        assert expected_atmosphere.h == test_atmosphere.h

    def test_gravitational_acceleration(self):
        atmosphere = self.make_test_atmosphere()

        gravitational_acceleration = atmosphere.gravitational_acceleration()
        expected_gravitational_acceleration = 9.81
        npt.assert_almost_equal(gravitational_acceleration, expected_gravitational_acceleration)


