import pytest
import numpy.testing as npt
import numpy as np

import data.constants as ct
from h2ermes_tools.reentry import gliding_trajectory as traj

class TestTrajectory:
    height = 400
    mass = 1000
    surface = 12
    C_D = 1.4
    C_L = 0.7
    latitude = 0
    longitude = 0
    atmosphere_model = "exponential"

    def make_test_trajectory(self):
        trajectory = traj.GlidingEntry.__new__(traj.GlidingEntry)
        trajectory.planet = "earth"
        trajectory.g_0 = ct.g_0
        trajectory.rho_0 = 1.225
        trajectory.R_e = ct.earth_radius
        trajectory.scale_height = 7200
        trajectory.beta = 1 / 7200
        trajectory.c_star = 1.1097e8
        trajectory.m = 3
        trajectory.altitude = self.height * 1000
        trajectory.alt_array = np.arange(0, trajectory.altitude+1)
        trajectory.mass = self.mass
        trajectory.S = self.surface
        trajectory.C_D = self.C_D
        trajectory.C_L = self.C_L
        trajectory.latitude = self.latitude
        trajectory.longitude = self.longitude
        trajectory.atmosphere_model = self.atmosphere_model

        return trajectory

    def test_trajectory_init(self):
        expected_trajectory = self.make_test_trajectory()

        test_trajectory = traj.GlidingEntry("earth",
                                            self.height,
                                            self.mass,
                                            self.surface,
                                            self.C_D,
                                            self.C_L,
                                            self.latitude,
                                            self.longitude,
                                            atmosphere=self.atmosphere_model)

        assert expected_trajectory.planet == test_trajectory.planet
        assert expected_trajectory.g_0 == test_trajectory.g_0
        assert expected_trajectory.rho_0 == test_trajectory.rho_0
        assert expected_trajectory.R_e == test_trajectory.R_e
        assert expected_trajectory.scale_height == test_trajectory.scale_height
        assert expected_trajectory.beta == test_trajectory.beta
        assert expected_trajectory.c_star == test_trajectory.c_star
        assert expected_trajectory.m == test_trajectory.m
        assert expected_trajectory.altitude == test_trajectory.altitude
        assert expected_trajectory.mass == test_trajectory.mass
        assert expected_trajectory.S == test_trajectory.S
        assert expected_trajectory.C_D == test_trajectory.C_D
        assert expected_trajectory.C_L == test_trajectory.C_L
        assert expected_trajectory.latitude == test_trajectory.latitude
        assert expected_trajectory.longitude == test_trajectory.longitude
        assert expected_trajectory.atmosphere_model == test_trajectory.atmosphere_model

    def test_trajectory_calculations(self):
        expected_trajectory = self.make_test_trajectory()
        expected_trajectory.circular_velocity_0 = np.sqrt(expected_trajectory.g_0 * (expected_trajectory.R_e + self.height * 1000))
        expected_trajectory.lift_parameter = expected_trajectory.mass * expected_trajectory.g_0 / (expected_trajectory.S * expected_trajectory.C_L)

        expected_circular_velocity = np.sqrt(ct.g_0 * (ct.earth_radius + self.height * 1000))
        expected_lift_parameter = 1167.458333

        npt.assert_almost_equal(expected_trajectory.circular_velocity_0, expected_circular_velocity)
        npt.assert_almost_equal(expected_trajectory.lift_parameter, expected_lift_parameter, decimal = 3)