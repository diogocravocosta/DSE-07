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
                           time_step: float = 0.01,
                           max_time: float = 100.0
                           ) -> tuple[pd.DataFrame, float]:
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

    density = cn.rho_0  # Assuming sea level standard density for simplicity

    # Initialize the initial conditions
    i = 0
    time = 0.
    past_velocity = -initial_velocity # up is defined as positive, so initial velocity is negative
    past_acceleration = 1 / 2 * past_velocity ** 2 * density / ballistic_coefficient + (
                    thrust_to_weight_ratio - 1) * cn.g_0
    past_distance = 0.

    data[i] = [time, past_velocity, past_acceleration, past_distance]

    for i in range(1, data.shape[0]):
        time += time_step

        # Update velocity, and distance
        curr_velocity = past_velocity + past_acceleration * time_step
        curr_distance = past_distance + past_velocity * time_step
        curr_acceleration = 1 / 2 * past_velocity ** 2 * density / ballistic_coefficient + (
                    thrust_to_weight_ratio - 1) * cn.g_0

        data[i] = [time, past_velocity, past_acceleration, past_distance] = [time, curr_velocity, curr_acceleration, curr_distance]

        if past_velocity >= 0:
            # trim the data to the last valid point
            data = data[:i + 1]
            break

    delta_v = data[-1, 0] * thrust_to_weight_ratio * cn.g_0  # Final delta_v is thrust-to-weight ratio times g_0 times time
    # Convert to a pandas dataframe for easier access
    return pd.DataFrame(data, columns=states), delta_v

def plot_landing_burn_heatmaps(ballistic_coefficients=None,
                               thrust_to_weight_ratios=None,
                               initial_velocities=None) -> None:
    """
    Plot heatmaps of landing burn delta v for different ballistic coefficients, thrust-to-weight ratios, and initial velocities.

    Args:
        ballistic_coefficients:
        thrust_to_weight_ratios:
        initial_velocities:
    """
    for bc in ballistic_coefficients:
        delta_vs = []
        for ttw in thrust_to_weight_ratios:
            row = []
            for iv in initial_velocities:
                _, delta_v = landing_burn_with_drag(iv, ttw, bc)
                row.append(delta_v)
            delta_vs.append(row)

        df = pd.DataFrame(delta_vs, index=thrust_to_weight_ratios, columns=initial_velocities)
        plt.figure(figsize=(10, 6))
        plt.title(f'Landing Burn Delta V Heatmap (BC={bc})')
        plt.xlabel('Initial Velocity (m/s)')
        plt.ylabel('Thrust-to-Weight Ratio')
        plt.imshow(df, aspect='auto', cmap='viridis', origin='lower',
                   extent=(initial_velocities[0], initial_velocities[-1],
                           thrust_to_weight_ratios[0], thrust_to_weight_ratios[-1]))
        plt.colorbar(label='Delta V (m/s)')
        plt.xticks(initial_velocities)
        plt.yticks(thrust_to_weight_ratios)
        plt.grid(False)
        plt.show()

if __name__ == "__main__":
    # Example usage
    initial_velocity = 100.0  # m/s
    thrust_to_weight_ratio = 5  # dimensionless
    ballistic_coefficient = 100  # kg/m^2, example value, float('inf') for no drag

    landing_burn_distance, delta_v = landing_burn_no_drag(initial_velocity, thrust_to_weight_ratio)
    print(f"Landing burn distance (no drag): {landing_burn_distance:.2f} m, Delta V: {delta_v:.2f} m/s")

    df, delta_v_with_drag = landing_burn_with_drag(initial_velocity, thrust_to_weight_ratio, ballistic_coefficient)
    print(f"Landing burn distance (with drag): {-df['distance'].iloc[-1]:.2f} m, Delta V: {delta_v_with_drag:.2f} m/s")

    # Plot velocity, acceleration, and distance over time in one plot with scaled y-axes
    plt.plot(df['time'], df['velocity'], label='Velocity (m/s)', color='blue')
    plt.plot(df['time'], df['acceleration'], label='Acceleration (m/s²)', color='orange')
    plt.plot(df['time'], df['distance'] - min(df['distance']), label='Distance (m)', color='green')
    plt.xlabel('Time (s)')
    plt.title('Landing Burn Profile')
    plt.legend()
    plt.grid()
    plt.show()

    # Uncomment the following line to generate delta v heatmaps
    # ballistic_coefficients = [float('inf'), 1000, 500, 200, 100]
    # thrust_to_weight_ratios = [1.5, 2.0, 2.5, 3.0, 3.5]
    # initial_velocities = [100, 150, 200, 250, 300]
    #
    # plot_landing_burn_heatmaps(ballistic_coefficients,
    #                            thrust_to_weight_ratios,
    #                            initial_velocities)
