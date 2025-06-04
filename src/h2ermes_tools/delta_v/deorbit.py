import numpy as np

import data.constants as cn

def calculate_orbit_given_flightpath_angle_and_velocity(flightpath_angle: float,
                                                        velocity: float,
                                                        altitude: float) -> float:
     """
     Calculate the orbital radius given the flight path angle and velocity.

     Args:
          flightpath_angle (float): Flight path angle in radians.
          velocity (float): Velocity in m/s.

     Returns:
          float: Orbital radius in meters.
     """
     if flightpath_angle == 0:
          raise ValueError("Flight path angle cannot be zero for this calculation.")

     radius = altitude + cn.earth_radius

     radial_velocity = velocity * np.sin(flightpath_angle)
     tangential_velocity = velocity * np.cos(flightpath_angle)

     semi_major_axis = -1/(velocity**2 / cn.gravitational_parameter - 2/radius)

     eccentricity = 1/cn.gravitational_parameter * np.sqrt((radius * tangential_velocity**2 - cn.gravitational_parameter)**2 + radial_velocity**2 * tangential_velocity**2 * radius**2)

     periapsis_radius = semi_major_axis * (1 - eccentricity)
     apoapsis_radius = semi_major_axis * (1 + eccentricity)

     return semi_major_axis, periapsis_radius, apoapsis_radius

def calculate_velocity_given_flightpath_angle_and_apogee(flightpath_angle: float,
                                                         apoapsis_radius: float,
                                                         altitude: float) -> float:
     """
     Calculate the velocity at a given flight path angle and apoapsis radius.

     Args:
          flightpath_angle (float): Flight path angle in radians.
          apoapsis_radius (float): Apoapsis radius in meters.

     Returns:
          float: Velocity in m/s.
     """


     return velocity

if __name__ == '__main__':
    # Example usage
    flightpath_angle = np.radians(-8) # Flight path angle in radians
    velocity = 7200  # m/s
    altitude = 100_000  # meters (400 km)

    semi_major_axis, periapsis_radius, apoapsis_radius = calculate_orbit_given_flightpath_angle_and_velocity(flightpath_angle, velocity, altitude)

    print(f"Semi-major axis: {semi_major_axis} m")
    print(f"Periapsis radius: {periapsis_radius} m")
    print(f"Apoapsis radius: {apoapsis_radius} m")

    print(f"Periapsis altitude: {periapsis_radius - cn.earth_radius} m")
    print(f"Apoapsis altittude: {apoapsis_radius - cn.earth_radius} m")
