import pandas as pd
import matplotlib.pyplot as plt
import numpy as np

import data.constants as cn

def calculate_circular_orbit_energy(altitude: float) -> float:
    """

    Args:
        altitude:

    Returns:

    """
    return -cn.gravitational_parameter/(2*(altitude + cn.earth_radius))

def delta_v_from_final_mass(final_mass: float, specific_impulse: float) -> float:
    """

    Args:
        specific_impulse:
        final_mass:

    Returns:

    """
    return cn.g_0 * specific_impulse * np.log(1/final_mass)

def plot_polar(df: pd.DataFrame):
    """

    Args:
        df: dataframe containing the trajectory data
    """
    thrust_df = df[df['throttle'] > 0]
    coast_df = df[df['throttle'] == 0]

    fig, ax = plt.subplots(subplot_kw={'projection': 'polar'})
    ax.plot(thrust_df['theta'], thrust_df['r']/1000, color='blue')
    ax.plot(coast_df['theta'], coast_df['r']/1000, color='red')
    ax.set_rmax(max(df['r']/1000) * 1.1)
    ax.grid(True)


    plt.show()

def plot_rectangular(df: pd.DataFrame):
    thrust_df = df[df['throttle'] > 0]
    coast_df = df[df['throttle'] == 0]
    plt.plot(thrust_df['theta']*cn.earth_radius/1000, (thrust_df['r']-cn.earth_radius)/1000, color='blue')
    plt.plot(coast_df['theta']*cn.earth_radius/1000, (coast_df['r']-cn.earth_radius)/1000, color='red')
    plt.xlabel('Downrange [km]')
    plt.ylabel('Altitude [km]')
    plt.grid(True)
    plt.show()

    plt.plot(thrust_df['theta'] * cn.earth_radius / 1000, (thrust_df['r'] - cn.earth_radius) / 1000, color='blue')
    plt.xlabel('Downrange [km]')
    plt.ylabel('Altitude [km]')
    plt.grid(True)
    plt.show()