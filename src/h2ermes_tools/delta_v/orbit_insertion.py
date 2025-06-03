import pandas as pd
import numpy as np
import time

import data.constants as cn
from h2ermes_tools.delta_v.helpers import delta_v_from_final_mass
from helpers import calculate_circular_orbit_energy, plot_polar, plot_rectangular

def simulate_ascent(initial_thrust_to_weight_ratio: float,
                    specific_impulse: float,
                    initial_horizontal_velocity: float,
                    initial_vertical_velocity: float,
                    initial_altitude: float,
                    timestep: float = 1.,
                    simulation_time:float = 600,
                    desired_orbital_altitude: float|None = None,
                    guidance: str = 'gravity turn',):
    """

    Args:
        guidance:
        desired_orbital_altitude:
        initial_thrust_to_weight_ratio:
        specific_impulse:
        initial_horizontal_velocity:
        initial_vertical_velocity:
        simulation_time:
        timestep:
        initial_altitude:

    Returns:

    """
    states = ['time', # s
              'T/W', # -, thrust to weight ratio
              'mass_ratio',  # -, Mass ratio is defined as current total mass over initial mass
              'r', # m, radial distance from Earth center
              'theta', # rad, angle between original Earth center-vehicle line to current Earth center-vehicle line
              'r_dot', # m/s, change of r
              'theta_dot', # rad/s, change of theta
              'throttle' # -, throttle setting, between 0 and 1
              ]

    data = pd.DataFrame(columns=states)

    # Initialize states
    i = 0
    mass_flow_non_dimensional = initial_thrust_to_weight_ratio / specific_impulse
    throttle = 1

    # Save initial values
    data.loc[0] = [0,
                   initial_thrust_to_weight_ratio,
                   1.,
                   initial_altitude + cn.earth_radius,
                   0.,
                   initial_vertical_velocity,
                   initial_horizontal_velocity / (initial_altitude + cn.earth_radius),
                   throttle]

    while (data.loc[i, 'time'] < simulation_time
           and data.loc[i, 'r'] >= cn.earth_radius/2
           and data.loc[i, 'mass_ratio'] > 0):
        past = data.loc[i]
        i += 1
        curr = pd.Series(index = states)

        curr['time'] = past['time'] + timestep
        curr['T/W'] = initial_thrust_to_weight_ratio / past['mass_ratio'] * throttle
        curr['mass_ratio'] = past['mass_ratio'] - mass_flow_non_dimensional * timestep * throttle
        curr['r'] = past['r'] + past['r_dot'] * timestep
        curr['theta'] = past['theta'] + past['theta_dot'] * timestep

        if guidance == 'gravity turn':
            pitch_angle = np.atan2(past['r_dot'], past['r'] * past['theta_dot']) # pitch angle is equal to flightpath angle for a gravity turn
        elif guidance == 'vertical':
            pitch_angle = np.pi/2

        vertical_acceleration = -cn.gravitational_parameter / past['r'] ** 2 + past['T/W'] * cn.g_0 * np.sin(pitch_angle)
        r_double_dot = vertical_acceleration + past['r'] * past['theta_dot'] ** 2

        horizontal_acceleration = past['T/W'] * cn.g_0 * np.cos(pitch_angle)
        theta_double_dot = (horizontal_acceleration - 2*past['r_dot']*past['theta_dot'])/past['r']

        curr['r_dot'] = past['r_dot'] + r_double_dot * timestep
        curr['theta_dot'] = past['theta_dot'] + theta_double_dot * timestep
        curr['throttle'] = throttle

        data.loc[i] = curr

        velocity = np.sqrt(past['r_dot']**2 + (past['r'] * past['theta_dot'])**2)
        energy = velocity**2/2 - cn.gravitational_parameter/past['r']

        if desired_orbital_altitude and energy > calculate_circular_orbit_energy(desired_orbital_altitude) and throttle > 0:
            throttle = 0
            timestep = 0.1
            print(data.loc[len(data)-1])
            print('energy', energy)
    print('energy', energy)
    return data

if __name__ == '__main__':

    start = time.time()
    trajectory = simulate_ascent(1.1,
                                 420,
                                 1820,
                                 1050,
                                 81700,
                                 timestep=0.1,
                                 simulation_time=4e2,
                                 desired_orbital_altitude=2e5)

    end = time.time()

    plot_polar(trajectory)
    plot_rectangular(trajectory)

    print(trajectory.loc[len(trajectory)-1])

    print("delta_v:", delta_v_from_final_mass(trajectory.loc[len(trajectory)-1, 'mass_ratio'], 420))
    print('time:', end-start)
