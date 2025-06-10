import pytest
import numpy.testing as npt
import numpy as np

import data.constants as ct
from h2ermes_tools.reentry import trajectory_simplified as traj

class TestTrajectory:
    height = 400
    mass = 1000
    surface = 12
    C_D = 1.4
    C_L = 0.7

    def make_test_trajectory(self):
        trajectory = traj.GlidingEntry.__new__(traj.GlidingEntry)
        trajectory.planet = "earth"
        trajectory.g = ct.g_0
        trajectory.rho_0 = 1.225
        trajectory.R = ct.earth_radius
        trajectory.scale_height = 7200
        trajectory.beta = 1 / 7200
        trajectory.c_star = 1.1097e8
        trajectory.m = 3
        trajectory.altitude = self.height * 1000
        trajectory.mass = self.mass
        trajectory.S = self.surface
        trajectory.C_D = self.C_D
        trajectory.C_L = self.C_L

        return trajectory

    def test_trajectory_init(self):
        expected_trajectory = self.make_test_trajectory()

        test_trajectory = traj.GlidingEntry("earth",
                                            self.height,
                                            self.mass,
                                            self.surface,
                                            self.C_D,
                                            self.C_L)

        assert expected_trajectory.planet == test_trajectory.planet
        assert expected_trajectory.g == test_trajectory.g
        assert expected_trajectory.rho_0 == test_trajectory.rho_0
        assert expected_trajectory.R == test_trajectory.R
        assert expected_trajectory.scale_height == test_trajectory.scale_height
        assert expected_trajectory.beta == test_trajectory.beta
        assert expected_trajectory.c_star == test_trajectory.c_star
        assert expected_trajectory.m == test_trajectory.m
        assert expected_trajectory.altitude == test_trajectory.altitude
        assert expected_trajectory.mass == test_trajectory.mass
        assert expected_trajectory.S == test_trajectory.S
        assert expected_trajectory.C_D == test_trajectory.C_D
        assert expected_trajectory.C_L == test_trajectory.C_L

    def test_trajectory_calculations(self):
        expected_trajectory = self.make_test_trajectory()
        expected_trajectory.circular_velocity_0 = np.sqrt(expected_trajectory.g * (expected_trajectory.R + self.height * 1000))

        expected_circular_velocity = np.sqrt(ct.g_0 * (ct.earth_radius + self.height * 1000))

        npt.assert_almost_equal(expected_trajectory.circular_velocity_0, expected_circular_velocity)