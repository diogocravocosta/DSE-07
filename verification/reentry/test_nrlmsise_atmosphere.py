import numpy.testing as npt

import data.constants as ct
from h2ermes_tools.reentry import atmosphere as atm

from pyatmos import download_sw_nrlmsise00,read_sw_nrlmsise00
from pyatmos import nrlmsise00

class TestNrlmsiseAtmosphere:
    height = 600
    time = '2014-07-22 22:18:45'
    latitude = 25
    longitude = 102

    def make_test_atmosphere(self):
        atmosphere = atm.NrlmsiseAtmosphere.__new__(atm.NrlmsiseAtmosphere)
        atmosphere.g_0 = ct.g_0  # m/s^2
        atmosphere.Re = ct.earth_radius  # m
        atmosphere.R_star = ct.R_star
        atmosphere.M = ct.molecular_mass
        atmosphere.R = atmosphere.R_star / atmosphere.M
        atmosphere.gamma = ct.gamma
        atmosphere.rho_0 = ct.density_sea_level
        atmosphere.altitude = self.height
        atmosphere.time = self.time  # UTC
        atmosphere.latitude = self.latitude
        atmosphere.longitude = self.longitude

        return atmosphere

    def test_atmosphere_init(self):
        expected_atmosphere = self.make_test_atmosphere()

        test_atmosphere = atm.NrlmsiseAtmosphere(self.height, self.time, self.latitude, self.longitude)

        assert expected_atmosphere.g_0 == test_atmosphere.g_0
        assert expected_atmosphere.Re == test_atmosphere.Re
        assert expected_atmosphere.R_star == test_atmosphere.R_star
        assert expected_atmosphere.M == test_atmosphere.M
        assert expected_atmosphere.R == test_atmosphere.R
        assert expected_atmosphere.gamma == test_atmosphere.gamma
        assert expected_atmosphere.rho_0 == test_atmosphere.rho_0
        assert expected_atmosphere.altitude == test_atmosphere.altitude
        assert expected_atmosphere.time == test_atmosphere.time
        assert expected_atmosphere.latitude == test_atmosphere.latitude
        assert expected_atmosphere.longitude == test_atmosphere.longitude

    def test_atmosphere_function(self):
        expected_atmosphere = self.make_test_atmosphere()
        expected_atmosphere.sw_data = expected_atmosphere.space_weather_data()

        expected_atmosphere.T_nrl, expected_atmosphere.rho_nrl = expected_atmosphere.atmosphere()

        expected_T = 765.8976564552341
        expected_rho = 1.714115212984513e-14

        npt.assert_almost_equal(expected_atmosphere.T_nrl, expected_T)
        npt.assert_almost_equal(expected_atmosphere.rho_nrl, expected_rho)






