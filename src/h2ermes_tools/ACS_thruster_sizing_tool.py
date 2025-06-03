import numpy as np

# Constants and parameters for the spacecraft
vehicle_mass = 60000  # kg
vehicle_shape = 'rectangular_prism'  # Shape of the spacecraft
vehicle_dimensions = np.array([10, 10, 5])  # Length, Width, Height in meters
gravitational_constant = 3.986 * 10**14  # m^3/s^2
theta = 50 # Angle between the spacecraft and the local vertical in degrees
altitude = 600 # Altitude of the spacecraft in km
COM = np.array([vehicle_dimensions[0] / 2,vehicle_dimensions[1] / 2, vehicle_dimensions[2] / 2])  # Center of mass of the spacecraft in meters from the geometric center
re_entry_moment = 166000/3 # Maximum moment during re-entry in Nm
redundancy_factor = 2 # Redundancy factor for thrusters

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
    },

    'nammo_220_3': {
        'thrust': 215,  # in N
        'Isp': 160,  # seconds
        'power': 50,  # in W
        'moment_arm': 5  # in meters
        },
    'nammo_220_4': {
        'thrust': 250,  # in N
        'Isp': 160,  # seconds
        'power': 50,  # in W
        'moment_arm': 5  # in meters
    },

}



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

def solar_radiation_pressure_torque (altitude, MMOI_vehicle):
    """
    Function to calculate the solar radiation pressure torque on the spacecraft.
    Args:
        altitude (float): Altitude of the spacecraft in kilometers.
        MMOI_vehicle (np.array): Mass moment of inertia of the spacecraft in kg*m^2
    Returns:
        solar_torque (float): Solar radiation pressure torque on the spacecraft in Nm.
    """
    altitude_r  = (altitude + 6371) * 1000  # Convert altitude to meters
    c_light = 3 * 10**8  # Speed of light in m/s
    vehicle_dipole_moment = 0.01  # in Am^2, assuming a small dipole moment for the spacecraft
    magnetic_constant = 7.8 * 10**15  # T*m^3
    vehicle_lambda = 2
    
    q_solar = np.average([0.15,0.8]) # Reflectance factor for solar radiation pressure torque
    phi_2 = 60  # in degrees
    solar_flux = 1361  # W/m^2, solar constant
    solar_area = vehicle_dimensions[0] * vehicle_dimensions[1]  # Area exposed to the sun

    geo_midpoint = np.array([vehicle_dimensions[0] / 2,vehicle_dimensions[1] / 2,0])  # Geometric midpoint of the spacecraft

    solar_torque = (solar_flux/c_light) * solar_area * (1 + q_solar) * (geo_midpoint - COM) * np.cos(np.deg2rad(phi_2))  # Solar radiation pressure torque
    solar_torque = np.linalg.norm(solar_torque)  # Convert to magnitude

    return solar_torque




# Uncomment the lines below to print the solar radiation pressure torque
solar_torque = solar_radiation_pressure_torque(altitude, MMOI_vehicle)
# print("Solar radiation pressure torque on the spacecraft:", solar_torque, "Nm")

def gravity_gradient_torque(altitude, theta, gravitational_constant, MMOI_vehicle):

    """
    Function to calculate the gravity gradient torque on the spacecraft.
    args:
        altitude (float): Altitude of the spacecraft in kilometers.
        gravitational_constant (float): Gravitational constant in m^3/s^2.
        MMOI_vehicle (np.array): Mass moment of inertia of the spacecraft in kg*m^2
    Returns:
        gravity_gradient_torque (float): Gravity gradient torque on the spacecraft in Nm.
    """
    altitude_r = (altitude + 6371) * 1000  # Convert altitude to meters
    theta_r = np.radians(theta)  # Convert angle to radians
    gravity_gradient_torque = 1.5 * (gravitational_constant / altitude_r**3) * np.abs(MMOI_vehicle[2] - MMOI_vehicle[1]) * np.sin(2*theta_r)  # Gravity gradient torque
    return gravity_gradient_torque

# Uncomment the lines below to print the gravity gradient torque
gravity_gradient_torque = gravity_gradient_torque(altitude, theta, gravitational_constant, MMOI_vehicle)
# print("Gravity gradient torque on the spacecraft:", gravity_gradient_torque, "Nm")

def aerodynamic_drag_torque(altitude, vehicle_dimensions, gravitational_constant):
    """    Function to calculate the aerodynamic drag on the spacecraft.
    Args:
        altitude (float): Altitude of the spacecraft in kilometers.
        vehicle_dimensions (np.array): Dimensions of the spacecraft in meters.
        gravitational_constant (float): Gravitational constant in m^3/s^2.
    Returns:
        aerodynamic_drag (float): Torque due to aerodynamic drag on the spacecraft in Nm.
    """
    altitude_r = (altitude + 6371) * 1000  # Convert altitude to meters
    vehicle_velocity_squared = gravitational_constant / altitude_r  # in m/s

    drag_coefficient = 2.2  # Dimensionless; typical value for spacecraft
    atmospheric_density = 1.03 * 10**-14  # kg/m^3
    aerodynamic_surface_area = vehicle_dimensions[0] * vehicle_dimensions[1]  # Area exposed to the flow

    aerodynamic_drag = 0.5 * atmospheric_density * drag_coefficient * aerodynamic_surface_area * vehicle_velocity_squared * (np.linalg.norm(COM - geo_midpoint))  # Aerodynamic drag torque
    return aerodynamic_drag

# Uncomment the lines below to print the aerodynamic drag torque
aerodynamic_drag = aerodynamic_drag_torque(altitude, vehicle_dimensions, gravitational_constant)
# print("Aerodynamic drag torque on the spacecraft:", aerodynamic_drag, "Nm")



def magnetic_torque(altitude):
    """
    Function to calculate the magnetic torque on the spacecraft.
    Args:
        altitude (float): Altitude of the spacecraft in kilometers.
    Returns:
        magnetic_torque (float): Magnetic torque on the spacecraft in Nm.
    """
    altitude_r = (altitude + 6371) * 1000  # Convert altitude to meters
    vehicle_dipole_moment = 0.01  # in Am^2, assuming a small dipole moment for the spacecraft
    magnetic_constant = 7.8 * 10**15  # T*m^3
    vehicle_lambda = 1.1  # Dimensionless factor for magnetic torque; ranges from 1 at the equator to 2 at the poles
    magnetic_torque = vehicle_dipole_moment * (magnetic_constant / altitude_r**3) * vehicle_lambda  # Magnetic torque
    return magnetic_torque
# Uncomment the lines below to print the magnetic torque
magnetic_torque = magnetic_torque(altitude)
# print("Magnetic torque on the spacecraft:", magnetic_torque, "Nm")

torque_list = [solar_torque, gravity_gradient_torque, aerodynamic_drag, magnetic_torque]  # List of all torques
print(torque_list)
disturbance_load = np.sum([magnetic_torque, gravity_gradient_torque, aerodynamic_drag, solar_torque])  # Maximum disturbance load


# print("Total disturbance torque on the spacecraft:", disturbance_load, "Nm")

# print("Torque required to counteract disturbances:", required_torque, "Nm")

ang_acc_x = 0.0025  # in rad/s^2
ang_acc_y = 0.003  # in rad/s^2
ang_acc_z = 0.0008  # in rad/s^2

ang_acc_max = max(ang_acc_x,ang_acc_y,ang_acc_z)  # in rad/s^2








def thruster_sizing(thrusters, ang_acc_max, MMOI_vehicle):
    """
    Function to size the thrusters based on the required torque and slew rate.
    Args:
        thrusters (dict): Dictionary containing thruster specifications.
        ang_acc_max (float): Maximum slew rate in rad/s^2.
        MMOI_vehicle (np.array): Mass moment of inertia of the spacecraft in kg*m^2.
    Returns:
        number_of_thrusters (dict): Dictionary containing the number of thrusters required for each type.
        thruster_moment_arm (dict): Dictionary containing the moment arm for each thruster type.
    """

    number_of_thrusters = {}
    thruster_moment_arm = {}

    torque_produced = np.array([MMOI_vehicle[0] * ang_acc_max, MMOI_vehicle[1] * ang_acc_max, MMOI_vehicle[2] * ang_acc_max]) + disturbance_load + re_entry_moment  # Torque produced by each thruster in Nm
    
    # thrust_x_direction_one = torque_produced[0] / (2 * thrusters['alphard_20']['moment_arm']) + disturbance_load  # Thrust required in the x-direction
    thrust_x_direction_two = torque_produced[0]/ (2 * thrusters['nammo_220']['moment_arm'])  # Thrust required in the x-direction
    thrust_x_direction_three = torque_produced[0]/ (2 * thrusters['nammo_220_3']['moment_arm'])  # Thrust required in the x-direction
    thrust_x_direction_four = torque_produced[0]/ (2 * thrusters['nammo_220_4']['moment_arm'])  # Thrust required in the x-direction
    # thrust_y_direction_one = torque_produced[1] / (2 * thrusters['alphard_20']['moment_arm']) + disturbance_load  # Thrust required in the y-direction
    thrust_y_direction_two = torque_produced[1] / (2 * thrusters['nammo_220']['moment_arm'])  # Thrust required in the y-direction
    thrust_y_direction_three = torque_produced[1] / (2 * thrusters['nammo_220_3']['moment_arm'])  # Thrust required in the y-direction
    thrust_y_direction_four = torque_produced[1] / (2 * thrusters['nammo_220_4']['moment_arm'])  # Thrust required in the y-direction
    # thrust_z_direction_one = torque_produced[2] / (2 * thrusters['alphard_20']['moment_arm']) + disturbance_load  # Thrust required in the z-direction
    thrust_z_direction_two = torque_produced[2] / (2 * thrusters['nammo_220']['moment_arm'])    # Calculate the number of thrusters required for each type
    thrust_z_direction_three = torque_produced[2] / (2 * thrusters['nammo_220_3']['moment_arm'])  # Thrust required in the z-direction
    thrust_z_direction_four = torque_produced[2] / (2 * thrusters['nammo_220_4']['moment_arm'])  # Thrust required in the z-direction  
    # number_of_thrusters['alphard_20'] = {
    #     'x': np.ceil((thrust_x_direction_two / thrusters['alphard_20']['thrust']) * redundancy_factor),
    #     'y': np.ceil((thrust_y_direction_two / thrusters['alphard_20']['thrust']) * redundancy_factor),
    #     'z': np.ceil((thrust_z_direction_two / thrusters['alphard_20']['thrust']) * redundancy_factor)
    # }

    number_of_thrusters['nammo_220'] = {
        'x': np.ceil((thrust_x_direction_two / thrusters['nammo_220']['thrust']) * redundancy_factor),
        'y': np.ceil((thrust_y_direction_two / thrusters['nammo_220']['thrust']) * redundancy_factor),
        'z': np.ceil((thrust_z_direction_two / thrusters['nammo_220']['thrust']) * redundancy_factor)
    }

    number_of_thrusters['nammo_220_3'] = {
    'x': np.ceil((thrust_x_direction_three / thrusters['nammo_220_3']['thrust']) * redundancy_factor),
    'y': np.ceil((thrust_y_direction_three / thrusters['nammo_220_3']['thrust']) * redundancy_factor),
    'z': np.ceil((thrust_z_direction_three / thrusters['nammo_220_3']['thrust']) * redundancy_factor)
    }

    number_of_thrusters['nammo_220_4'] = {
        'x': np.ceil((thrust_x_direction_four / thrusters['nammo_220_4']['thrust']) * redundancy_factor),
        'y': np.ceil((thrust_y_direction_four / thrusters['nammo_220_4']['thrust']) * redundancy_factor),
        'z': np.ceil((thrust_z_direction_four / thrusters['nammo_220_4']['thrust']) * redundancy_factor)
}

    # Calculate the moment arm for each thruster type
    # thruster_moment_arm['alphard_20'] = {
    #     'x': torque_produced[0] / (number_of_thrusters['alphard_20']['x'] * thrusters['alphard_20']['thrust']),
    #     'y': torque_produced[1] / (number_of_thrusters['alphard_20']['y'] * thrusters['alphard_20']['thrust']),
    #     'z': torque_produced[2] / (number_of_thrusters['alphard_20']['z'] * thrusters['alphard_20']['thrust'])
    # }

    thruster_moment_arm['nammo_220'] = {
        'x': torque_produced[0] / (number_of_thrusters['nammo_220']['x'] * thrusters['nammo_220']['thrust']), 
        'y': torque_produced[1] / (number_of_thrusters['nammo_220']['y'] * thrusters['nammo_220']['thrust']),
        'z': torque_produced[2] / (number_of_thrusters['nammo_220']['z'] * thrusters['nammo_220']['thrust'])
    }
    thruster_moment_arm['nammo_220_3'] = {
    'x': torque_produced[0] / (number_of_thrusters['nammo_220_3']['x'] * thrusters['nammo_220_3']['thrust']),
    'y': torque_produced[1] / (number_of_thrusters['nammo_220_3']['y'] * thrusters['nammo_220_3']['thrust']),
    'z': torque_produced[2] / (number_of_thrusters['nammo_220_3']['z'] * thrusters['nammo_220_3']['thrust'])
    }
    thruster_moment_arm['nammo_220_4'] = {
        'x': torque_produced[0] / (number_of_thrusters['nammo_220_4']['x'] * thrusters['nammo_220_4']['thrust']),       
        'y': torque_produced[1] / (number_of_thrusters['nammo_220_4']['y'] * thrusters['nammo_220_4']['thrust']),
        'z': torque_produced[2] / (number_of_thrusters['nammo_220_4']['z'] * thrusters['nammo_220_4']['thrust'])
    }
    power_thrusters = {
        # 'alphard_20': number_of_thrusters['alphard_20']['x'] * thrusters['alphard_20']['power'] + \
        #               number_of_thrusters['alphard_20']['y'] * thrusters['alphard_20']['power'] + \
        #               number_of_thrusters['alphard_20']['z'] * thrusters['alphard_20']['power'],
        'nammo_220': number_of_thrusters['nammo_220']['x'] * thrusters['nammo_220']['power'] + \
                     number_of_thrusters['nammo_220']['y'] * thrusters['nammo_220']['power'] + \
                     number_of_thrusters['nammo_220']['z'] * thrusters['nammo_220']['power']
    }

    return number_of_thrusters, thruster_moment_arm, power_thrusters

number_of_thrusters, thruster_moment_arm, power_thrusters = thruster_sizing(thrusters, ang_acc_max, MMOI_vehicle)
print("Number of thrusters required:", number_of_thrusters)
print("Moment arm for each thruster type:", thruster_moment_arm['nammo_220_4'])
# print("Power required for each thruster type:", power_thrusters)

print("Coordinates of the center of mass of H2ermes: ", COM)



















