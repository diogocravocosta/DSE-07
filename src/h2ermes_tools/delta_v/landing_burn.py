import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

import data.constants as cn

def landing_burn_no_drag(initial_velocity: float,
                         thrust_to_weight_ratio: float,
                         ) -> tuple[float, float]:
    """
    Calculate the landing burn distance and delta v without considering drag assuming purely vertical motion.

    Args:
        initial_velocity (float): Initial velocity in m/s.
        thrust_to_weight_ratio (float): Thrust-to-weight ratio (dimensionless).

    Returns:
        tuple: (landing_burn_distance, delta_v) where:
            landing_burn_distance (float): Distance required for landing burn in meters.
            delta_v (float): Change in velocity during the landing burn in m/s.
    """
    if thrust_to_weight_ratio <= 1:
        raise ValueError("Thrust-to-weight ratio must be greater than 1 for a landing burn.")

    time = initial_velocity / ((thrust_to_weight_ratio - 1) * cn.g_0)
    landing_burn_distance = initial_velocity * time - 0.5 * (thrust_to_weight_ratio - 1) * cn.g_0 * time**2
    delta_v = thrust_to_weight_ratio * cn.g_0 * time

    return landing_burn_distance, delta_v

def landing_burn_with_drag(initial_velocity: float,
                           thrust_to_weight_ratio: float,
                           ballistic_coefficient: float,
                           time_step: float = 0.1,
                           max_time: float = 100.0
                           ):
    """
    Simulate the landing burn including drag.

    Assumptions:
    - The landing burn is vertical.
    - The mass is constant during the burn.
    - The thrust is constant.
    - The drag coefficient is constant.
    - The density of air is constant (using sea level standard density for simplicity).

    Args:
        initial_velocity(float): Initial velocity in m/s.
        thrust_to_weight_ratio(float): Thrust-to-weight ratio (dimensionless).
        ballistic_coefficient(float): Ballistic coefficient in kg/m^2.
        Defined as mass divided by cross-sectional area times drag coefficient.
        time_step(float): Time step for the simulation in seconds.
        max_time(float): Maximum time for the simulation in seconds.

    Returns:

    """
    if thrust_to_weight_ratio <= 1:
        raise ValueError("Thrust-to-weight ratio must be greater than 1 for a landing burn.")

    states = ['time', 'velocity', 'acceleration', 'distance']

    data = np.empty((int(max_time / time_step) + 1, len(states)))

    # Initialize the initial conditions
    i = 0
    time = 0.
    past_velocity = -initial_velocity
    past_acceleration = 0.
    past_distance = 0.

    density = cn.rho_0  # Assuming sea level standard density for simplicity

    data[i] = [time, past_velocity, past_acceleration, past_distance]

    for i in range(1, data.shape[0]):
        time += time_step

        # Update velocity, and distance
        curr_velocity = past_velocity + past_acceleration * time_step
        curr_distance = past_distance + past_velocity * time_step
        curr_acceleration = 1 / 2 * past_velocity ** 2 * density * ballistic_coefficient + (
                    thrust_to_weight_ratio - 1) * cn.g_0

        data[i] = [time, past_velocity, past_acceleration, past_distance] = [time, curr_velocity, curr_acceleration, curr_distance]

        if past_velocity <= 0:
            # trim the data to the last valid point
            data = data[:i + 1]
            break

    # Convert to a pandas dataframe for easier access
    return pd.DataFrame(data, columns=states)

