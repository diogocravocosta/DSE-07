import pandas as pd
import numpy as np
import time

import data.constants as cn
from h2ermes_tools.delta_v.helpers import delta_v_from_final_mass, calculate_circular_orbit_energy, plot_polar, plot_rectangular

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

    data = np.empty((int(simulation_time//timestep + 3), len(states)))

    # Initialize states
    time = 0
    past_tw = initial_thrust_to_weight_ratio
    past_mass = 1.
    past_r = initial_altitude + cn.earth_radius
    past_theta = 0.
    past_r_dot = initial_vertical_velocity
    past_theta_dot = initial_horizontal_velocity/past_r
    throttle = 1

    mass_flow_non_dimensional = initial_thrust_to_weight_ratio / specific_impulse

    # Save initial values
    data[0] = [time,
               past_tw,
               past_mass,
               past_r,
               past_theta,
               past_r_dot,
               past_theta_dot,
               throttle]
    if desired_orbital_altitude:
        desired_orbit_energy = calculate_circular_orbit_energy(desired_orbital_altitude)

    for i in range(len(data)):
        time += timestep
        curr_tw = initial_thrust_to_weight_ratio / past_mass * throttle
        curr_mass = past_mass - mass_flow_non_dimensional * timestep * throttle
        curr_r = past_r + past_r_dot * timestep
        curr_theta = past_theta + past_theta_dot * timestep

        if guidance == 'gravity turn':
            pitch_angle = np.atan2(past_r_dot, past_r * past_theta_dot) # pitch angle is equal to flightpath angle for a gravity turn
        elif guidance == 'vertical':
            pitch_angle = np.pi/2

        vertical_acceleration = -cn.gravitational_parameter / past_r ** 2 + past_tw * cn.g_0 * np.sin(pitch_angle)
        r_double_dot = vertical_acceleration + past_r * past_theta_dot ** 2

        horizontal_acceleration = past_tw * cn.g_0 * np.cos(pitch_angle)
        theta_double_dot = (horizontal_acceleration - 2*past_r_dot*past_theta_dot)/past_r

        curr_r_dot = past_r_dot + r_double_dot * timestep
        curr_theta_dot = past_theta_dot + theta_double_dot * timestep

        data[i] = (
            time,
            past_tw,
            past_mass,
            past_r,
            past_theta,
            past_r_dot,
            past_theta_dot,
            throttle
         ) = (
            time,
            curr_tw,
            curr_mass,
            curr_r,
            curr_theta,
            curr_r_dot,
            curr_theta_dot,
            throttle
        )


        velocity = np.sqrt(past_r_dot**2 + (past_r * past_theta_dot)**2)
        energy = velocity**2/2 - cn.gravitational_parameter/past_r

        if desired_orbital_altitude and energy > desired_orbit_energy and throttle > 0:
            throttle = 0
            print(pd.Series(data[i], index = states))
            print('energy', energy)

        if past_r < cn.earth_radius/2 or past_mass < 0:
            break

    print('energy', energy)
    return pd.DataFrame(data, columns=states)

if __name__ == '__main__':

    start = time.time()
    trajectory = simulate_ascent(1.18,
                                 420,
                                 2280,
                                 1040,
                                 81700,
                                 timestep=0.01,
                                 simulation_time=5e3,
                                 desired_orbital_altitude=2e5)

    end = time.time()

    plot_polar(trajectory)
    plot_rectangular(trajectory)

    print(trajectory.loc[len(trajectory)-1])

    print("delta_v:", delta_v_from_final_mass(trajectory.loc[len(trajectory)-1, 'mass_ratio'], 420))
    print('time:', end-start)
