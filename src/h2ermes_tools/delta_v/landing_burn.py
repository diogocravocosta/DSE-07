import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

import data.constants as cn
from h2ermes_tools.delta_v.helpers import delta_v_from_final_mass

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
                           initial_thrust_to_weight_ratio: float,
                           ballistic_coefficient: float,
                           specific_impulse: float,
                           time_step: float = 0.01,
                           max_time: float = 100.0
                           ) -> pd.DataFrame:
    """
    Simulate the landing burn including drag.

    Assumptions:
    - The landing burn is vertical.
    - The thrust is constant.
    - The drag coefficient is constant.
    - The density of air is constant (using sea level standard density for simplicity).

    Args:
        initial_velocity(float): Initial velocity in m/s.
        initial_thrust_to_weight_ratio(float): Thrust-to-weight ratio (dimensionless).
        ballistic_coefficient(float): Ballistic coefficient in kg/m^2.
        Defined as mass divided by cross-sectional area times drag coefficient.
        time_step(float): Time step for the simulation in seconds.
        max_time(float): Maximum time for the simulation in seconds.

    Returns:

    """
    if initial_thrust_to_weight_ratio <= 1:
        raise ValueError("Thrust-to-weight ratio must be greater than 1 for a landing burn.")

    states = ['time', 'velocity', 'acceleration', 'distance', 'mass_fraction', 'thrust_to_weight_ratio']

    data = np.empty((int(max_time / time_step) + 1, len(states)))

    density = cn.density_sea_level  # Assuming sea level standard density for simplicity
    mass_flow_non_dimensional = initial_thrust_to_weight_ratio / specific_impulse  # Non-dimensional mass flow rate

    # Initialize the initial conditions
    i = 0
    time = 0.
    past_velocity = -initial_velocity # up is defined as positive, so initial velocity is negative
    past_acceleration = 1 / 2 * past_velocity ** 2 * density / ballistic_coefficient + (
            initial_thrust_to_weight_ratio - 1) * cn.g_0
    past_distance = 0.
    past_mass_fraction = 1.0  # Mass fraction is defined as current mass over initial mass, so it starts at 1
    past_thrust_to_weight_ratio = initial_thrust_to_weight_ratio

    data[i] = [time, past_velocity, past_acceleration, past_distance, past_mass_fraction, past_thrust_to_weight_ratio]

    for i in range(1, data.shape[0]):
        time += time_step

        # Update velocity, and distance
        curr_velocity = past_velocity + past_acceleration * time_step
        curr_distance = past_distance + past_velocity * time_step
        curr_acceleration = 1 / 2 * past_velocity ** 2 * density / ballistic_coefficient + (
                past_thrust_to_weight_ratio - 1) * cn.g_0
        # Update mass fraction and thrust-to-weight ratio
        curr_mass_fraction = past_mass_fraction - mass_flow_non_dimensional * time_step
        curr_thrust_to_weight_ratio = initial_thrust_to_weight_ratio / past_mass_fraction

        data[i] = [
            time,
            past_velocity,
            past_acceleration,
            past_distance,
            past_mass_fraction,
            past_thrust_to_weight_ratio
        ] = [
            time,
            curr_velocity,
            curr_acceleration,
            curr_distance,
            curr_mass_fraction,
            curr_thrust_to_weight_ratio
        ]

        if past_velocity >= 0:
            # trim the data to the last valid point
            data = data[:i + 1]
            break

    # Convert to a pandas dataframe for easier access
    return pd.DataFrame(data, columns=states)

if __name__ == "__main__":
    # Example usage
    initial_velocity = 400  # m/s
    thrust_to_weight_ratio = 2.9  # dimensionless
    ballistic_coefficient = 2600  # kg/m^2, example value, float('inf') for no drag
    specific_impulse = 361  # seconds, example value

    landing_burn_distance, delta_v = landing_burn_no_drag(initial_velocity, thrust_to_weight_ratio)
    print(f"Landing burn distance (no drag): {landing_burn_distance:.2f} m, Delta V: {delta_v:.2f} m/s")

    df = landing_burn_with_drag(initial_velocity, thrust_to_weight_ratio, ballistic_coefficient, specific_impulse)
    delta_v_with_drag = delta_v_from_final_mass(df['mass_fraction'].iloc[-1], specific_impulse)
    print(f"Landing burn distance (with drag): {-df['distance'].iloc[-1]:.2f} m, "
          f"Delta V: {delta_v_with_drag:.2f} m/s, "
          f"Time: {df['time'].iloc[-1]:.2f} s")

    # Plot velocity, acceleration, and distance over time in one plot with scaled y-axes
    plt.plot(df['time'], df['velocity'], label='Velocity (m/s)', color='blue')
    plt.plot(df['time'], df['acceleration']*10, label='Acceleration (0.1m/s²)', color='orange')
    plt.plot(df['time'], (df['distance'] - min(df['distance']))/10, label='Distance (10m)', color='green')
    plt.xlabel('Time (s)')
    plt.title('Landing Burn Profile')
    plt.legend()
    plt.grid()
    plt.show()
