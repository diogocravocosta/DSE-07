import pytest
import numpy.testing as npt
import numpy as np

import data.constants as ct
from h2ermes_tools.reentry import trajectory_simplified as traj


class TestRange:
    lift_drag_ratio = 1.84
    V_c_0 = 7800
    entry_speed = 0.8 * V_c_0

    def make_test_trajectory(self):
        trajectory = traj.GlidingEntry.__new__(traj.GlidingEntry)
        trajectory.lift_drag_ratio = self.lift_drag_ratio
        trajectory.V_c_0 = self.V_c_0
        trajectory.entry_speed = self.entry_speed

        return trajectory

    def test_range_ratio(self):
        # need: entry_speed, V_c_0, lift_drag_ratio
        range_trajectory = self.make_test_trajectory()
        range_trajectory.flight_range_ratio = range_trajectory.get_flight_range_ratio()

        expected_range_ratio = 0.940
        npt.assert_almost_equal(range_trajectory.flight_range_ratio, expected_range_ratio, decimal = 2)


class TestEqAngle:
    scale_height = 7166.44
    beta = 1 / scale_height
    R_e = 6378137
    V_c_0 = 7800
    V = 0.1 * V_c_0
    lift_drag_ratio = 2

    def make_test_trajectory(self):
        trajectory = traj.GlidingEntry.__new__(traj.GlidingEntry)
        trajectory.lift_drag_ratio = self.lift_drag_ratio
        trajectory.V_c_0 = self.V_c_0
        trajectory.V = self.V
        trajectory.beta = self.beta
        trajectory.R_e = self.R_e
        trajectory.v_ratio = trajectory.V / trajectory.V_c_0

        return trajectory

    def test_eq_angle(self):
        equilibrium_angle_trajectory = self.make_test_trajectory()

        equilibrium_angle_trajectory.flight_path_angle_equilibrium_rad = equilibrium_angle_trajectory.equilibrium_glide_angle()
        equilibrium_angle_trajectory.flight_path_angle_equilibrium_deg = equilibrium_angle_trajectory.flight_path_angle_equilibrium_rad * 180/np.pi
        expected_radians = -0.112
        expected_degrees = -6.44

        npt.assert_almost_equal(equilibrium_angle_trajectory.flight_path_angle_equilibrium_deg, expected_degrees, decimal = 2)
        npt.assert_almost_equal(equilibrium_angle_trajectory.flight_path_angle_equilibrium_rad, expected_radians, decimal = 2)

class TestFlightTime:
    V_c = 7920
    V_entry = 0.8 * V_c
    lift_drag_ratio = 2

    def make_test_trajectory(self):
        trajectory = traj.GlidingEntry.__new__(traj.GlidingEntry)
        trajectory.lift_drag_ratio = self.lift_drag_ratio
        trajectory.V_c_0 = self.V_c
        trajectory.entry_speed = self.V_entry
        trajectory.v_ratio = trajectory.entry_speed / trajectory.V_c_0
        trajectory.g_0 = ct.g_0

        return trajectory

    def test_flight_time(self):
        flight_time_trajectory = self.make_test_trajectory()
        flight_time_trajectory.t_h, flight_time_trajectory.t_min, flight_time_trajectory.t_sec = flight_time_trajectory.get_flight_time()

        expected_time_min = 29.6

        npt.assert_almost_equal(np.round(flight_time_trajectory.t_min, decimals = 1), expected_time_min, decimal = 2)


class TestDeceleration:
    scale_height = 7166.44
    beta = 1 / scale_height
    R_e = 6378137
    V_c_0 = 7800
    V = 0.1 * V_c_0
    lift_drag_ratio = 2
    lift_parameter = 6000

    def make_test_trajectory(self):
        trajectory = traj.GlidingEntry.__new__(traj.GlidingEntry)
        trajectory.lift_drag_ratio = self.lift_drag_ratio
        trajectory.V_c_0 = self.V_c_0
        trajectory.entry_speed = self.V
        trajectory.v_ratio = trajectory.entry_speed / trajectory.V_c_0
        trajectory.g_0 = ct.g_0
        trajectory.beta = self.beta
        trajectory.R_e = self.R_e
        trajectory.lift_parameter = self.lift_parameter
        trajectory.rho_0 = 1.225

        return trajectory

    def test_deceleration(self):
        deceleration_trajectory = self.make_test_trajectory()

        deceleration_trajectory.a, deceleration_trajectory.a_max = deceleration_trajectory.get_gliding_deceleration()

        deceleration_trajectory.v_ratio_a_max, deceleration_trajectory.V_a_max = deceleration_trajectory.get_max_deceleration_v_ratio()

        expected_V_a_max = 1698.26

        deceleration_trajectory.rho_a_max, deceleration_trajectory.h_a_max = deceleration_trajectory.get_max_deceleration_altitude()

        expected_a_max = 0.4525

        expected_v_ratio_a_max = 0.2177

        expected_rho =0.003963
        expected_h = 41089

        npt.assert_almost_equal(deceleration_trajectory.v_ratio_a_max, expected_v_ratio_a_max, decimal=2)
        npt.assert_almost_equal(deceleration_trajectory.V_a_max, expected_V_a_max, decimal=2)
        npt.assert_almost_equal(deceleration_trajectory.a_max, expected_a_max, decimal = 2)
        npt.assert_almost_equal(deceleration_trajectory.rho_a_max, expected_rho, decimal=2)
        npt.assert_almost_equal(deceleration_trajectory.h_a_max, expected_h, decimal = 0)













