import numpy as np

# Constants and parameters for the spacecraft
vehicle_mass = 1000  # kg
vehicle_shape = 'rectangular_prism'  # Shape of the spacecraft
vehicle_dimensions = np.array([10, 15, 5])  # Length, Width, Height in meters

altitude = 600 # Altitude of the spacecraft in km
COM = np.array([vehicle_dimensions[0] / 2,vehicle_dimensions[1] / 2, vehicle_dimensions[2] / 2])  # Center of mass of the spacecraft in meters from the geometric center

redundancy_factor = 2 # Redundancy factor for thrusters



# velocity = np.sqrt(gravitational_constant / ((altitude + 6371) * 1000))  # in m/s
# print("Velocity of the spacecraft at altitude", altitude, "km:", velocity, "m/s")

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
print("Mass Moment of Inertia (MMOI) of the spacecraft:", MMOI_vehicle)

geo_midpoint = np.array([vehicle_dimensions[0] / 2,vehicle_dimensions[1] / 2,0])

# print("Mass Moment of Inertia (MMOI) of the spacecraft:", MMOI_vehicle)
# print("MMOI in the x-direction:", MMOI_vehicle[0], "kg*m^2")


def disturbance_loads(altitude, MMOI_vehicle):
    """
    Function to calculate the disturbance loads on the spacecraft.
    Args:
        altitude (float): Altitude of the spacecraft in kilometers.
        MMOI_vehicle (np.array): Mass moment of inertia of the spacecraft in kg*m^2
        gravitational_constant (float): Gravitational constant in m^3/s^2.
    Returns:
        disturbance_loads (list): List of all disturbance loads on the spacecraft.
    """
    altitude_r  = (altitude + 6371) * 1000  # Convert altitude to meters
    gravitational_constant = 3.986 * 10**14  # m^3/s^2
    theta = 20 # Random angle in degrees
    c_light = 3 * 10**8  # Speed of light in m/s
    vehicle_velocity = np.sqrt(gravitational_constant / altitude_r)  # in m/s
    theta_r = np.radians(theta)  # Random angle in radians


    vehicle_dipole_moment = 0.01  # in Am^2, assuming a small dipole moment for the spacecraft
    magnetic_constant = 7.8 * 10**15  # T*m^3
    vehicle_lambda = 2
    
    drag_coefficient = 2.2  # Dimensionless
    atmospheric_density = 1.225 * np.exp(-altitude_r / 8500)  # kg/m^3

    geo_midpoint = np.array([vehicle_dimensions[0] / 2,vehicle_dimensions[1] / 2,0])  # Geometric midpoint of the spacecraft
    phi_1 = 0.1 # in radians
    phi_2 = 0.2 # in radians
    solar_flux = 1361  # W/m^2, solar constant
    solar_area = vehicle_dimensions[0] * vehicle_dimensions[1]  # Area exposed to the sun

    magnetic_torque = vehicle_dipole_moment * (magnetic_constant / altitude_r**3) * vehicle_lambda
    gravity_gradient_torque = 1.5 * (gravitational_constant / altitude_r**3) * np.abs(MMOI_vehicle[2] - MMOI_vehicle[0]) * np.sin(2 * theta_r)
    aerodynamic_drag = 0.5 * atmospheric_density * drag_coefficient * vehicle_velocity**2 * vehicle_dimensions[0] * vehicle_dimensions[1]
    solar_torque = solar_area * (phi_1/c_light) * (1 + solar_flux) * np.cos(phi_2) * (geo_midpoint - COM)  # Solar radiation pressure torque
    solar_torque = np.linalg.norm(solar_torque)  # Convert to magnitude

    disturbance_loads = [gravity_gradient_torque, magnetic_torque, aerodynamic_drag, solar_torque]
    return disturbance_loads

gravity_gradient_torque, magnetic_torque, aerodynamic_drag, solar_torque = disturbance_loads(altitude, MMOI_vehicle)
# print(disturbance_loads(altitude, MMOI_vehicle))
# Calculate the required torque to counteract the disturbances
required_torque = sum(np.abs([gravity_gradient_torque, magnetic_torque, aerodynamic_drag, solar_torque]))
# print("Torque required to counteract disturbances:", required_torque, "Nm")

slew_rate_x = 0.5  # in rad/s^2
slew_rate_y = 0.5  # in rad/s^2
slew_rate_z = 0.57  # in rad/s^2

slew_rate_max = max(slew_rate_x,slew_rate_y,slew_rate_z)  # in rad/s^2

maneuver_time = 120  # in seconds



thrusters = {
    'alphard_20': {
        'thrust': 10,  # in N
        'Isp': 160,  # seconds
        'power': 18,  # in W
        'moment_arm': 5  # in meters
    },
    'nammo_220': {
        'thrust': 100,  # in N
        'Isp': 160,  # seconds
        'power': 50,  # in W
        'moment_arm': 5  # in meters
    }
}


def thruster_sizing(thrusters, slew_rate_maximum, MMOI_vehicle):
    """
    Function to size the thrusters based on the required torque and slew rate.
    Args:
        thrusters (dict): Dictionary containing thruster specifications.
        slew_rate_maximum (float): Maximum slew rate in rad/s^2.
        MMOI_vehicle (np.array): Mass moment of inertia of the spacecraft in kg*m^2.
    Returns:
        number_of_thrusters (dict): Dictionary containing the number of thrusters required for each type.
        thruster_moment_arm (dict): Dictionary containing the moment arm for each thruster type.
    """

    number_of_thrusters = {}
    thruster_moment_arm = {}

    torque_produced = np.array([MMOI_vehicle[0] * slew_rate_maximum, MMOI_vehicle[1] * slew_rate_maximum, MMOI_vehicle[2] * slew_rate_maximum])  # Torque produced by each thruster in Nm
    
    thrust_x_direction_one = torque_produced[0] / (2 * thrusters['alphard_20']['moment_arm']) + required_torque  # Thrust required in the x-direction
    thrust_x_direction_two = torque_produced[0]/ (2 * thrusters['nammo_220']['moment_arm']) + required_torque  # Thrust required in the x-direction
    thrust_y_direction_one = torque_produced[1] / (2 * thrusters['alphard_20']['moment_arm']) + required_torque  # Thrust required in the y-direction
    thrust_y_direction_two = torque_produced[1] / (2 * thrusters['nammo_220']['moment_arm']) + required_torque  # Thrust required in the y-direction
    thrust_z_direction_one = torque_produced[2] / (2 * thrusters['alphard_20']['moment_arm']) + required_torque  # Thrust required in the z-direction
    thrust_z_direction_two = torque_produced[2] / (2 * thrusters['nammo_220']['moment_arm']) + required_torque    # Calculate the number of thrusters required for each type  
    number_of_thrusters['alphard_20'] = {
        'x': np.ceil(2 * (thrust_x_direction_two / thrusters['alphard_20']['thrust']) * redundancy_factor),
        'y': np.ceil(2 * (thrust_y_direction_two / thrusters['alphard_20']['thrust']) * redundancy_factor),
        'z': np.ceil(2 * (thrust_z_direction_two / thrusters['alphard_20']['thrust']) * redundancy_factor)
    }

    number_of_thrusters['nammo_220'] = {
        'x': np.ceil(2 * (thrust_x_direction_two / thrusters['nammo_220']['thrust']) * redundancy_factor),
        'y': np.ceil(2 * (thrust_y_direction_two / thrusters['nammo_220']['thrust']) * redundancy_factor),
        'z': np.ceil(2 * (thrust_z_direction_two / thrusters['nammo_220']['thrust']) * redundancy_factor)
    }

    # Calculate the moment arm for each thruster type
    thruster_moment_arm['alphard_20'] = {
        'x': torque_produced[0] / (number_of_thrusters['alphard_20']['x'] * thrusters['alphard_20']['thrust']),
        'y': torque_produced[1] / (number_of_thrusters['alphard_20']['y'] * thrusters['alphard_20']['thrust']),
        'z': torque_produced[2] / (number_of_thrusters['alphard_20']['z'] * thrusters['alphard_20']['thrust'])
    }

    thruster_moment_arm['nammo_220'] = {
        'x': torque_produced[0] / (number_of_thrusters['nammo_220']['x'] * thrusters['nammo_220']['thrust']),
        'y': torque_produced[1] / (number_of_thrusters['nammo_220']['y'] * thrusters['nammo_220']['thrust']),
        'z': torque_produced[2] / (number_of_thrusters['nammo_220']['z'] * thrusters['nammo_220']['thrust'])
    }

    power_thrusters = {
        'alphard_20': number_of_thrusters['alphard_20']['x'] * thrusters['alphard_20']['power'] + \
                      number_of_thrusters['alphard_20']['y'] * thrusters['alphard_20']['power'] + \
                      number_of_thrusters['alphard_20']['z'] * thrusters['alphard_20']['power'],
        'nammo_220': number_of_thrusters['nammo_220']['x'] * thrusters['nammo_220']['power'] + \
                     number_of_thrusters['nammo_220']['y'] * thrusters['nammo_220']['power'] + \
                     number_of_thrusters['nammo_220']['z'] * thrusters['nammo_220']['power']
    }

    return number_of_thrusters, thruster_moment_arm, power_thrusters

number_of_thrusters, thruster_moment_arm, power_thrusters = thruster_sizing(thrusters, slew_rate_max, MMOI_vehicle)
print("Number of thrusters required:", number_of_thrusters)
print("Moment arm for each thruster type:", thruster_moment_arm)
print("Power required for each thruster type:", power_thrusters)




# moment_arm = required_torque/np.array([thrusters['nammo_220']['thrust']])
# print(moment_arm)















