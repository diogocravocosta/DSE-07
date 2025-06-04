import pytest
import numpy.testing as npt

import data.constants as ct
from h2ermes_tools import atmosphere as atm

class TestAtmosphere:
    height = 400

    def make_test_atmosphere(self):
        atmosphere = atm.Atmosphere.__new__(atm.Atmosphere)
        atmosphere.g_0 = ct.g_0  # m/s^2
        atmosphere.Re = ct.earth_radius  # m
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

        atmosphere.g = atmosphere.gravitational_acceleration()

        # to calculate expected_gravitational_acceleration use: g0  / (1 + h / Re)**2 with h and Re in meters
        expected_gravitational_acceleration = 8.68335742

        npt.assert_almost_equal(atmosphere.g, expected_gravitational_acceleration)

    def test_geopotential_altitude(self):
        atmosphere = self.make_test_atmosphere()

        atmosphere.geop_altitude = atmosphere.geopotential_altitude()

        # to calculate expected_geop_altitude use: h * (1 - h/Re) with h and Re in meters
        expected_geop_altitude = 374914.3049138

        npt.assert_almost_equal(atmosphere.geop_altitude, expected_geop_altitude)

    def test_exponential_atmosphere(self):
        atmosphere = self.make_test_atmosphere()
        atmosphere.T_exp, atmosphere.scale_height, atmosphere.beta, atmosphere.speed_of_sound = atmosphere.exponential_atmosphere(scale_height=7050)

        # following calculations come from the scale_height, which can take values between 7050m and 7200 according to prof. Mooij
        expected_scale_height = 7050 # input value

        # to calculate expected_T_exp use: floor(scale_height * g_0 / R) where g0 = 9.80665 and R = 287
        expected_T_exp = 240
        # to calculate expected_beta use: 1 / scale_height
        expected_beta = 0.0001418439716
        # to calculate expected_speed_of_sounds use: round(sqrt(gamma * temperature * R), 1), with gamma = 1.4 and R = 287
        expected_speed_of_sound = 310.6

        atmosphere.rho_exp = atmosphere.exponential_density()

        # to calculate expected_rho_exp use: rho_0 * exp(-beta * height)
        expected_rho_exp = 2.8010084592781267e-25

        npt.assert_almost_equal(atmosphere.T_exp, expected_T_exp)
        npt.assert_almost_equal(atmosphere.scale_height, expected_scale_height)
        npt.assert_almost_equal(atmosphere.beta, expected_beta)
        npt.assert_almost_equal(atmosphere.speed_of_sound, expected_speed_of_sound)
        npt.assert_almost_equal(atmosphere.rho_exp, expected_rho_exp)


