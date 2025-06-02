import numpy as np
vehicle_mass = 1000  # kg
vehicle_shape = 'rectangular_prism'  # Shape of the spacecraft
vehicle_dimensions = np.array([5.0, 3.0, 2.0])  # Length, Width, Height in meters
c_light = 3e8  # Speed of light in m/s
altitude = 600 # Altitude of the spacecraft in km
COM = 0  # Center of mass of the spacecraft in meters from the geometric center
def space_craft_properties(vehicle_mass):
    """
    Function to define the spacecraft properties.
    Returns:
        MMOI_vehicle (float): Mass moment of inertia of the spacecraft in kg.
    """
    # Calculate the mass moment of inertia (MMOI) for a rectangular prism
    length, width, height = vehicle_dimensions
    Ixx = (1/12) * vehicle_mass * (width**2 + height**2)
    Iyy = (1/12) * vehicle_mass * (length**2 + height**2)
    Izz = (1/12) * vehicle_mass * (length**2 + width**2)
    MMOI_vehicle = np.array([Ixx, Iyy, Izz])
    return MMOI_vehicle

MMOI_vehicle = space_craft_properties(vehicle_mass)


# print("Mass Moment of Inertia (MMOI) of the spacecraft:", MMOI_vehicle)
# print("MMOI in the x-direction:", MMOI_vehicle[0], "kg*m^2")


def disturbance_loads(altitude, MMOI_vehicle):
    """
    Function to calculate the disturbance loads on the spacecraft.
    Args:
        altitude (float): Altitude of the spacecraft in meters.
        MMOI_vehicle (np.array): Mass moment of inertia of the spacecraft in kg*m^2.
    Returns:
        disturbance_loads (list): List of all disturbance loads on the spacecraft.
    """
    altitude = (altitude + 6371) * 1000  # in meters
    graviational_constant = 3.986e14  # m^3/s^2
    vehicle_velocity = np.sqrt(graviational_constant) / altitude  # in m/s
    theta = np.radians(20)  # Random angle in radians
    gravity_gradient_torque = 1.5 * (graviational_constant / altitude**3) * np.abs(MMOI_vehicle[2] - MMOI_vehicle[0]) * np.sin(2 * theta)
    vehicle_dipole = 0.1  # in radians
    magnetic_constant = 7.8e15  # T*m^3
    vehicle_lambda = 2
    magnetic_torque = vehicle_dipole * (magnetic_constant / altitude**3) * vehicle_lambda
    drag_coefficient = 2.2  # Dimensionless
    atmospheric_density = 1.225 * np.exp(-altitude / 8500)  # kg/m^3
    aerodynamic_drag = 0.5 * atmospheric_density * drag_coefficient * vehicle_velocity**2 * vehicle_dimensions[0] * vehicle_dimensions[1]
    geo_midpoint = (vehicle_dimensions[0] / 2)**2 + (vehicle_dimensions[1] / 2)**2  # Geometric midpoint of the spacecraft
    phi_1, phi_2, solar_flux = 0.1, 0.2, 1361  # Solar flux in W/m^2
    solar_area = vehicle_dimensions[0] * vehicle_dimensions[1]  # Area exposed to the sun
    solar_torque = solar_area * (phi_1/c_light) * (1+solar_flux) * np.cos(phi_2) * (geo_midpoint - COM)  # Solar radiation pressure torque


    disturbance_loads = [gravity_gradient_torque, magnetic_torque, aerodynamic_drag, solar_torque]
    return disturbance_loads

gravity_gradient_torque, magnetic_torque, aerodynamic_drag, solar_torque = disturbance_loads(altitude, MMOI_vehicle)
# print(disturbance_loads(altitude, MMOI_vehicle))
# Calculate the required torque to counteract the disturbances
required_torque = max(disturbance_loads(altitude, MMOI_vehicle))
# print("Torque required to counteract disturbances:", required_torque, "Nm")





