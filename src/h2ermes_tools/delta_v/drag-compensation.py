import numpy as np

from h2ermes_tools.atmosphere import Atmosphere
import data.constants as cn

def drag_compensation_delta_v(altitude: float,
                              ballistic_coefficient: float,
                              time_in_orbit: float) -> float:
    """
    Computes the delta-v required for drag compensation in orbit.
    Args:
        altitude: orbital altitude in meters
        ballistic_coefficient: ballistic coefficient in kg/m^2
        time_in_orbit: time in orbit in seconds
    Returns:
        delta_v: required delta-v in m/s
    """
    # Create an instance of the Atmosphere class
    density = 2.91E-10  # Get the atmospheric density at the given altitude
    # todo update density calculation

    circular_velocity = np.sqrt(cn.gravitational_parameter / (cn.earth_radius + altitude))
    acceleration = 1 / 2 * density * circular_velocity ** 2 / ballistic_coefficient

    return acceleration * time_in_orbit

if __name__ == "__main__":
    # Example usage
    altitude = 200000  # 400 km
    ballistic_coefficient = 100  # kg/m^2
    time_in_orbit = 3600 * 24 # 1 hour in seconds

    delta_v = drag_compensation_delta_v(altitude, ballistic_coefficient, time_in_orbit)
    print(f"Required delta-v for drag compensation: {delta_v:.4f} m/s")