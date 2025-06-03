import pandas as pd
import numpy as np
import time
import matplotlib.pyplot as plt

import data.constants as cn
import h2ermes_tools.delta_v.helpers as hv

def simulate_ascent(initial_thrust_to_weight_ratio: float,
                    specific_impulse: float,
                    initial_horizontal_velocity: float,
                    initial_vertical_velocity: float,
                    initial_altitude: float,
                    timestep: float = 1.,
                    simulation_time:float = 600,
                    target_orbital_altitude: float|None = None,
                    guidance: str = 'gravity turn',):
    """

    Args:
        guidance:
        target_orbital_altitude:
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
    else:
        target_orbit_energy = None

    if target_orbital_altitude is None and (guidance != 'gravity turn' or guidance != 'vertical'):
        raise ValueError('this guidance requires target orbital altitude')

    if guidance == 'linear tangent' or guidance == 'another tangent':
        target_orbit_velocity = np.sqrt(cn.gravitational_parameter / (target_orbital_altitude + cn.earth_radius))

        delta_v_estimate = 5700 #m/s, calculate based on inputs

        final_mass_estimate = np.e**(-delta_v_estimate/(specific_impulse * cn.g_0))
        burn_time_estimate = (1-final_mass_estimate)/mass_flow_non_dimensional
    elif guidance == 'altitude':
        slope = flightpath_angle / (target_orbital_altitude - initial_altitude)

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

        if guidance == 'gravity turn':
            pitch_angle = flightpath_angle + np.deg2rad(7.81) # offset from flightpath angle needs to be manually adjusted
        elif guidance == 'vertical':
            pitch_angle = np.pi/2
        # elif guidance == 'linear tangent':
        #     target_flightpath_angle = np.atan2(cn.g_0 * (burn_time_estimate - time), target_orbit_velocity)
        #
        #     pitch_gain = 0.5
        #     pitch_angle = target_flightpath_angle + (target_flightpath_angle-flightpath_angle)*pitch_gain
        # elif guidance == 'altitude':
        #     target_flightpath_angle = slope * (target_orbital_altitude - past_r + cn.earth_radius)
        #
        #     pitch_gain = 15
        #     pitch_angle = target_flightpath_angle + (target_flightpath_angle-flightpath_angle)*pitch_gain
        # elif guidance == 'another tangent':
        #     pitch_angle = np.atan((1-time/burn_time_estimate) * np.tan(initial_pitch_angle*1.095)) # The additional guidance profiles don't work currently

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


        velocity = np.sqrt(past_r_dot**2 + (past_r * past_theta_dot)**2)
        energy = velocity**2/2 - cn.gravitational_parameter/past_r

        if target_orbital_altitude and throttle > 0 and energy > target_orbit_energy:
            throttle = 0
            print(pd.Series(data[i], index = states))

        if past_r < cn.earth_radius/2 or past_mass < 0:
            break

    return pd.DataFrame(data, columns=states)

if __name__ == '__main__':

    start = time.time()
    specific_impulse = 450

    trajectory = simulate_ascent(0.88,
                                 specific_impulse,
                                 2280,
                                 1040,
                                 81700,
                                 timestep=0.01,
                                 simulation_time=5e3,
                                 target_orbital_altitude=2e5,
                                 guidance='gravity turn',
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

    print("delta_v:", hv.delta_v_from_final_mass(trajectory.loc[len(trajectory)-1, 'mass_ratio'], specific_impulse))
    print('time:', end-start)
