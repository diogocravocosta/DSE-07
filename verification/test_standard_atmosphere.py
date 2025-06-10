import numpy.testing as npt
import numpy as np

import data.constants as ct
from h2ermes_tools.reentry import atmosphere as atm

class TestStandardAtmosphere:
    height_test = 400

    def make_test_atmosphere(self,height):
        atmosphere = atm.StandardAtmosphere.__new__(atm.StandardAtmosphere)
        atmosphere.g_0 = ct.g_0  # m/s^2
        atmosphere.Re = ct.earth_radius  # m
        atmosphere.R_star = ct.R_star
        atmosphere.M = ct.molecular_mass
        atmosphere.R = atmosphere.R_star / atmosphere.M
        atmosphere.gamma = ct.gamma
        atmosphere.rho_0 = ct.density_sea_level
        atmosphere.T_0 = ct.temperature_sea_level
        atmosphere.p_0 = ct.pressure_sea_level
        atmosphere.T_inf = 1000
        atmosphere.R_0 = 6356766
        atmosphere.altitude = height * 1000

        return atmosphere

    def test_atmosphere_init(self):
        expected_atmosphere = self.make_test_atmosphere(self.height_test)

        test_atmosphere = atm.StandardAtmosphere(self.height_test)

        assert expected_atmosphere.g_0 == test_atmosphere.g_0
        assert expected_atmosphere.Re == test_atmosphere.Re
        assert expected_atmosphere.R_star == test_atmosphere.R_star
        assert expected_atmosphere.M == test_atmosphere.M
        assert expected_atmosphere.R == test_atmosphere.R
        assert expected_atmosphere.gamma == test_atmosphere.gamma
        assert expected_atmosphere.rho_0 == test_atmosphere.rho_0
        assert expected_atmosphere.altitude == test_atmosphere.altitude

    def test_geopotential_altitude(self):
        atmosphere = self.make_test_atmosphere(self.height_test)

        atmosphere.altitude_gp = atmosphere.geopotential_altitude()

        # to calculate expected_geop_altitude use: h * (1 - h/Re) with h and Re in meters
        expected_geop_altitude = 374914.3049138

        npt.assert_almost_equal(atmosphere.altitude_gp, expected_geop_altitude)

    def test_exp_gravitational_acceleration(self):
        atmosphere = self.make_test_atmosphere(self.height_test)

        atmosphere.g = atmosphere.gravitational_acceleration_std()

        # to calculate expected_gravitational_acceleration use: g0  / (1 + h / Re)**2 with h and Re in meters
        expected_gravitational_acceleration = 9.226097114

        npt.assert_almost_equal(atmosphere.g, expected_gravitational_acceleration)

    def test_all_temperatures(self):
        atm_dict = {10: 223.252,
                    15: 216.65,
                    25: 221.552,
                    35: 236.513,
                    49: 270.65,
                    60: 247.021,
                    75: 208.399,
                    89: 186.87,
                    100: 195.08,
                    115: 300.00,
                    200: 854.56,
                    900: 1000.00}

        for height, expected_temp in atm_dict.items():
            atmosphere = self.make_test_atmosphere(height)
            atmosphere.altitude_gp = atmosphere.geopotential_altitude()
            temperature, pressure, rho = atmosphere.atmosphere()

            expected_temperature = expected_temp

            npt.assert_almost_equal(temperature, expected_temperature, decimal=2)


    def test_all_pressures(self):
        atm_dict = {10: 26499,
                    15: 12111,
                    25: 2549.2,
                    35: 574.59,
                    49: 90.336,
                    60: 21.958,
                    75: 2.3881,
                    89: 0.21919,
                    100: 0.032011,
                    115: 0.0040096,
                    200: 0.000084736,
                    900: 0.000000010873
                    }

        for height, expected_p in atm_dict.items():
            atmosphere = self.make_test_atmosphere(height)
            atmosphere.altitude_gp = atmosphere.geopotential_altitude()
            temperature, pressure, rho = atmosphere.atmosphere()

            expected_pressure = expected_p

            npt.assert_almost_equal(pressure, expected_pressure, decimal=0)




    def test_all_densities(self):
        atm_dict = {10: 4.1351e-1,
                    15: 1.9476e-1,
                    25: 4.0084e-2,
                    35: 8.4634e-3,
                    49: 1.1628e-3,
                    60: 3.0968e-4,
                    75: 3.9921e-5,
                    89: 4.081e-6,
                    100: 5.604e-7,
                    115: 4.289e-8,
                    200: 2.541e-10,
                    900: 5.759e-15}

        for height, expected_rho in atm_dict.items():
            atmosphere = self.make_test_atmosphere(height)
            atmosphere.altitude_gp = atmosphere.geopotential_altitude()
            temperature, pressure, rho = atmosphere.atmosphere()

            expected_density = expected_rho

            npt.assert_almost_equal(rho, expected_density, decimal=2)
