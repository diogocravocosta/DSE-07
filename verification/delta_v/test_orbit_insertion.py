import numpy as np

import data.constants as cn
from h2ermes_tools.delta_v.orbit_insertion import simulate_ascent

def test_orbit_insertion_circular_orbit():
    altitude = 2e5
    velocity = np.sqrt(cn.gravitational_parameter/(altitude + cn.earth_radius))

    trajectory_df = simulate_ascent(0, 1, velocity, 0, altitude, timestep=10, simulation_time=140 * 60)

    np.testing.assert_allclose(trajectory_df['r'], altitude+cn.earth_radius, atol=1)

def test_vertical_ballistic_trajectory():
    altitude = 0
    vertical_velocity = 1000

    final_altitude = -1 / (vertical_velocity ** 2 / 2 / cn.gravitational_parameter - 1 / cn.earth_radius) - cn.earth_radius # expected altitude calculated using mechanical energy balance

    trajectory_df = simulate_ascent(0, 1, 0, vertical_velocity, altitude, timestep=0.1, simulation_time=110)

    assert abs(max(trajectory_df['r']) - cn.earth_radius - final_altitude) < 50

def test_vertical_powered_flight():
    time = 20

    thrust_to_weight = 1.1
    specific_impulse = 300
    mass_flow = thrust_to_weight / specific_impulse

    final_mass = 1 - mass_flow * time

    expected_burnout_altitude = -specific_impulse*cn.g_0 * (final_mass*np.log(1/final_mass) - mass_flow*time)/mass_flow - 1/2*cn.g_0*time**2
    expected_burnout_velocity = specific_impulse*cn.g_0*np.log(1/final_mass) - cn.g_0*time

    trajectory_df = simulate_ascent(thrust_to_weight,
                                    specific_impulse,
                                    0,
                                    0,
                                    0,
                                    0.01,
                                    time,
                                    guidance = 'vertical')

    burnout_altitude = trajectory_df.loc[len(trajectory_df)-1, 'r'] - cn.earth_radius
    burnout_velocity = trajectory_df.loc[len(trajectory_df)-1, 'r_dot']

    assert burnout_altitude > expected_burnout_altitude # simulation should be higher as analytic calculation assumes constant gravity
    assert abs(burnout_altitude - expected_burnout_altitude) < 5
    assert burnout_velocity > expected_burnout_velocity
    assert abs(burnout_velocity - expected_burnout_velocity) < 1

