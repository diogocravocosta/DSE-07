import pandas as pd
import numpy as np
import time
import matplotlib.pyplot as plt

import data.constants as cn
import h2ermes_tools.delta_v.helpers as hv

np.seterr(all='raise')  # Raise exceptions for all numpy errors to catch them

def simulate_ascent(initial_thrust_to_weight_ratio: float,
                    specific_impulse: float,
                    initial_horizontal_velocity: float,
                    initial_vertical_velocity: float,
                    initial_altitude: float,
                    timestep: float = 1.,
                    simulation_time:float = 600,
                    target_orbital_altitude: float|None = None,
                    guidance: str = 'gravity turn',
                    gravity_turn_offset: float = 0,
                    enable_terminal_guidance: bool = False) -> pd.DataFrame:
    """
    Simulates the ascent from arbitrary initial conditions until an energy equivalent
    to the energy of a circular orbit at the target orbital altitude.

    Args:
        enable_terminal_guidance(bool): if True, the ascent will use terminal guidance to adjust the pitch angle at the end of orbit insertion
        gravity_turn_offset: offset from flightpath angle in degrees for gravity turn guidance in radians
        guidance (str): defines the guidance scheme used to control pitch angle. currently only gravity turn with offset and vertical burn are implemented
        target_orbital_altitude:
        initial_thrust_to_weight_ratio:
        specific_impulse:
        initial_horizontal_velocity:
        initial_vertical_velocity:
        simulation_time:
        timestep:
        initial_altitude:

    Returns:
    pd.DataFrame logging the states throughout the ascent
    """
    states = ['time', # s
              'T/W', # -, thrust to weight ratio
              'mass_ratio',  # -, Mass ratio is defined as current total mass over initial mass
              'r', # m, radial distance from Earth center
              'theta', # rad, angle between original Earth center-vehicle line to current Earth center-vehicle line
              'r_dot', # m/s, change of r
              'theta_dot', # rad/s, change of theta
              'throttle', # -, throttle setting, between 0 and 1
              'pitch_angle', # rad
              'flightpath_angle', # rad
              'target_flightpath_angle' #rad
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
    flightpath_angle = np.atan2(past_r_dot, past_r * past_theta_dot)
    pitch_angle = initial_pitch_angle = flightpath_angle
    target_flightpath_angle = flightpath_angle

    mass_flow_non_dimensional = initial_thrust_to_weight_ratio / specific_impulse

    if target_orbital_altitude:
        target_orbit_energy = hv.calculate_circular_orbit_energy(target_orbital_altitude)
        target_orbit_velocity = np.sqrt(cn.gravitational_parameter / (cn.earth_radius + target_orbital_altitude))
    else:
        target_orbit_energy = None

    # Save initial values
    data[0] = [time,
               past_tw,
               past_mass,
               past_r,
               past_theta,
               past_r_dot,
               past_theta_dot,
               throttle,
               pitch_angle,
               flightpath_angle,
               target_flightpath_angle]

    for i in range(len(data)):
        time += timestep
        curr_tw = initial_thrust_to_weight_ratio / past_mass * throttle
        curr_mass = past_mass - mass_flow_non_dimensional * timestep * throttle
        curr_r = past_r + past_r_dot * timestep
        curr_theta = past_theta + past_theta_dot * timestep

        flightpath_angle = np.atan2(past_r_dot, past_r * past_theta_dot)

        velocity = np.sqrt(past_r_dot**2 + (past_r * past_theta_dot)**2)
        energy = velocity**2/2 - cn.gravitational_parameter/past_r

        if (enable_terminal_guidance
                and target_orbital_altitude
                and abs(past_r - (cn.earth_radius + target_orbital_altitude))/(cn.earth_radius + target_orbital_altitude) < 0.05
                and velocity > target_orbit_velocity * 0.9
                and throttle > 0):
            # Terminal guidance
            target_r_double_dot = -past_r_dot / 10
            try:
                target_pitch_angle = np.arcsin((target_r_double_dot - past_r * past_theta_dot ** 2 + cn.gravitational_parameter/past_r**2) / (cn.g_0 * past_tw))
            except FloatingPointError:
                target_pitch_angle = np.pi/2

            pitch_angle = max(target_pitch_angle, 0) # Don't pitch below horizontal
        elif guidance == 'gravity turn':
            pitch_angle = flightpath_angle + gravity_turn_offset # offset from flightpath angle needs to be manually adjusted
            pitch_angle = max(pitch_angle, 0)  # Don't pitch below horizontal
        elif guidance == 'vertical':
            pitch_angle = np.pi/2

        vertical_acceleration = -cn.gravitational_parameter / past_r ** 2 + past_tw * cn.g_0 * np.sin(pitch_angle)
        r_double_dot = vertical_acceleration + past_r * past_theta_dot ** 2

        horizontal_acceleration = past_tw * cn.g_0 * np.cos(pitch_angle)
        theta_double_dot = (horizontal_acceleration - 2*past_r_dot*past_theta_dot)/past_r

        curr_r_dot = past_r_dot + r_double_dot * timestep
        curr_theta_dot = past_theta_dot + theta_double_dot * timestep

        data[i] = (
            _,
            past_tw,
            past_mass,
            past_r,
            past_theta,
            past_r_dot,
            past_theta_dot,
            throttle,
            _,
            _,
            _
         ) = (
            time,
            curr_tw,
            curr_mass,
            curr_r,
            curr_theta,
            curr_r_dot,
            curr_theta_dot,
            throttle,
            pitch_angle,
            flightpath_angle,
            target_flightpath_angle
        )

        if target_orbital_altitude and throttle > 0 and energy > target_orbit_energy:
            throttle = 0
            print(pd.Series(data[i], index = states))

        if past_r < cn.earth_radius/2 or past_mass < 0:
            data = data[:i+1]  # Truncate data to the current length
            break

    return pd.DataFrame(data, columns=states)

if __name__ == '__main__':

    start = time.time()
    specific_impulse = 450

    trajectory = simulate_ascent(1.4,
                                 specific_impulse,
                                 2280,
                                 1040,
                                 81700,
                                 timestep=0.01,
                                 simulation_time=5e3,
                                 target_orbital_altitude=2e5,
                                 guidance='gravity turn',
                                 gravity_turn_offset=np.deg2rad(-4),
                                 enable_terminal_guidance=True
                                 )

    end = time.time()

    hv.plot_polar(trajectory)
    hv.plot_rectangular(trajectory)

    thrust_df = trajectory[trajectory['throttle'] > 0]

    plt.plot(thrust_df['time'], np.rad2deg(thrust_df['pitch_angle']), label='pitch', color='red')
    plt.plot(thrust_df['time'], np.rad2deg(thrust_df['flightpath_angle']), label='flightpath', color='orange')
    # plt.plot(thrust_df['time'], np.rad2deg(thrust_df['target_flightpath_angle']), label='target flightpath', color='green')
    # plt.plot(thrust_df['time'],
    #          np.rad2deg(-thrust_df['flightpath_angle'] + thrust_df['target_flightpath_angle']), label='diff', color='purple')
    plt.plot(thrust_df['time'], (thrust_df['r'] - cn.earth_radius)/10000, label='altitude [10km]', color='blue')
    plt.xlabel('Time [s]')
    plt.ylabel('Angle [deg]')
    plt.legend()
    plt.grid(True)
    plt.show()



    print(trajectory.loc[len(trajectory)-1])

    print("delta_v:", hv.delta_v_from_final_mass(trajectory.loc[len(trajectory) - 1, 'mass_ratio'], specific_impulse))
    print('time:', end-start)

    """
    Given
    specific_impulse = 450
    initial_horizontal_velocity = 2280
    initial_vertical_velocity = 1040
    initial_altitude = 81700
    
    All of the following use terminal guidance enabled
    
    The delta V is 6210 m/s with thrust to weight ratio of 0.6 and 26.5 degrees offset from gravity turn guidance.
    The delta V is 5937 m/s with thrust to weight ratio of 0.7 and 17 degrees offset from gravity turn guidance.
    The delta V is 5810 m/s with thrust to weight ratio of 0.8 and 11 degrees offset from gravity turn guidance.
    The delta V is 5739 m/s with thrust to weight ratio of 0.9 and 6.5 degrees offset from gravity turn guidance.
    The delta V is 5673 m/s with thrust to weight ratio of 1.0 and 3.7 degrees offset from gravity turn guidance.
    The delta V is 5645 m/s with thrust to weight ratio of 1.1 and 1.5 degrees offset from gravity turn guidance.
    The delta V is 5629 m/s with thrust to weight ratio of 1.2 and 0 degrees offset from gravity turn guidance.
    The delta V is 5618 m/s with thrust to weight ratio of 1.3 and -1.5 degrees offset from gravity turn guidance.
    The delta V is 5610 m/s with thrust to weight ratio of 1.4 and -3.5 degrees offset from gravity turn guidance.
    """

    known_data = {'delta_v': [6210, 5937, 5810, 5739, 5673, 5645, 5629, 5618, 5610],
                  'thrust_to_weight_ratio': [0.6, 0.7, 0.8, 0.9, 1.0, 1.1, 1.2, 1.3, 1.4],}