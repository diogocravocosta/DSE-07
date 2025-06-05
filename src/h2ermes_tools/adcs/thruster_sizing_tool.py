import numpy as np
import data.constants as cn
from h2ermes_tools.variables import mass_vehicle, launch_vehicle_dimensions

# Constants and parameters for the spacecraft
vehicle_mass = mass_vehicle.value  # kg
vehicle_shape = "rectangular_prism"  # Shape of the spacecraft (Supposed to be a cone)
# Spacecraft_dimensions, needs to be imported from constants file
vehicle_dimensions = launch_vehicle_dimensions.value  # Length, Width, Height in meters
gravitational_constant = cn.gravitational_parameter  # m^3/s^2
# MMOI_vehicle # To be imported from constants file
theta = 50  # Angle between the spacecraft and the local vertical in degrees
altitude = 600  # Altitude of the spacecraft in km
COM = np.array(
    [vehicle_dimensions[0] / 2, vehicle_dimensions[1] / 2, vehicle_dimensions[2] / 2]
)  # Center of mass of the spacecraft in meters
re_entry_moment = 166000  # Maximum moment during re-entry in Nm
redundancy_factor = 2  # Redundancy factor for thrusters

thrusters = {
    "alphard_20": {
        "thrust": 10,  # in N
        "Isp": 160,  # seconds
        "power": 18,  # in W
        "moment_arm": 5,  # in meters
    },
    "nammo_220": {  # sea level thrust
        "thrust": 148,  # in N
        "Isp": 160,  # seconds
        "power": 50,  # in W
        "moment_arm": 5,  # in meters
    },
    "nammo_220_3": {  # nominal vacuum thrust
        "thrust": 220,  # in N
        "Isp": 160,  # seconds
        "power": 50,  # in W
        "moment_arm": 5,  # in meters
    },
    "nammo_220_4": {  # max vacuum thrust
        "thrust": 250,  # in N
        "Isp": 160,  # seconds
        "power": 50,  # in W
        "moment_arm": 5,  # in meters
    },
}


# velocity = np.sqrt(gravitational_constant / ((altitude + 6371) * 1000))  # in m/s
# print("Velocity of the spacecraft at altitude", altitude, "km:", velocity, "m/s")


def space_craft_properties(
    vehicle_mass,
):  # MMOI Values from constants file will be used for actual calculations
    """
    Function to define the spacecraft properties.

    Returns:
        MMOI_vehicle (float): Mass moment of inertia of the spacecraft in kg.
    """
    # Calculate the mass moment of inertia (MMOI) for a rectangular prism
    length, width, height = vehicle_dimensions
    Ixx = (1 / 12) * vehicle_mass * (width**2 + height**2)
    Iyy = (1 / 12) * vehicle_mass * (length**2 + height**2)
    Izz = (1 / 12) * vehicle_mass * (length**2 + width**2)
    MMOI_vehicle = np.array([Ixx, Iyy, Izz])
    return MMOI_vehicle


MMOI_vehicle = space_craft_properties(vehicle_mass)
print("Mass Moment of Inertia (MMOI) of the spacecraft:", MMOI_vehicle)

geo_midpoint = np.array(
    [vehicle_dimensions[0] / 2, vehicle_dimensions[1] / 2, 0]
)  # Actual value to be determined from the constants file

# print("Mass Moment of Inertia (MMOI) of the spacecraft:", MMOI_vehicle)
# print("MMOI in the x-direction:", MMOI_vehicle[0], "kg*m^2")


def solar_radiation_pressure_torque(altitude):
    """
    Function to calculate the solar radiation pressure torque on the spacecraft.
    Args:
        altitude (float): Altitude of the spacecraft in kilometers.
    Returns:
        solar_torque (float): Solar radiation pressure torque on the spacecraft in Nm.
    """
    c_light = cn.speed_of_light  # Speed of light in m/s

    q_solar = np.average(
        [0.15, 0.8]
    )  # Reflectance factor for solar radiation pressure torque
    phi_2 = 60  # in degrees
    solar_flux = cn.solar_constant  # W/m^2, solar constant
    solar_area = (
        vehicle_dimensions[0] * vehicle_dimensions[1]
    )  # Area exposed to the sun

    geo_midpoint = np.array(
        [vehicle_dimensions[0] / 2, vehicle_dimensions[1] / 2, 0]
    )  # Geometric midpoint of the spacecraft

    solar_torque = (
        (solar_flux / c_light)
        * solar_area
        * (1 + q_solar)
        * (geo_midpoint - COM)
        * np.cos(np.deg2rad(phi_2))
    )  # Solar radiation pressure torque

    solar_torque = np.linalg.norm(solar_torque)  # Convert to magnitude

    return solar_torque


# Uncomment the lines below to print the solar radiation pressure torque
solar_torque = solar_radiation_pressure_torque(altitude)
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
    gravity_gradient_torque = (
        1.5
        * (gravitational_constant / altitude_r**3)
        * np.abs(MMOI_vehicle[2] - MMOI_vehicle[1])
        * np.sin(2 * theta_r)
    )  # Gravity gradient torque
    return gravity_gradient_torque


# Uncomment the lines below to print the gravity gradient torque
gravity_gradient_torque = gravity_gradient_torque(
    altitude, theta, gravitational_constant, MMOI_vehicle
)
# print("Gravity gradient torque on the spacecraft:", gravity_gradient_torque, "Nm")


def aerodynamic_drag_torque(altitude, vehicle_dimensions, gravitational_constant):
    """Function to calculate the aerodynamic drag on the spacecraft.
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
    atmospheric_density = cn.orbit_altitude_density  # kg/m^3
    aerodynamic_surface_area = (
        vehicle_dimensions[0] * vehicle_dimensions[1]
    )  # Area exposed to the flow

    aerodynamic_drag = (
        0.5
        * atmospheric_density
        * drag_coefficient
        * aerodynamic_surface_area
        * vehicle_velocity_squared
        * (np.linalg.norm(COM - geo_midpoint))
    )  # Aerodynamic drag torque
    return aerodynamic_drag


# Uncomment the lines below to print the aerodynamic drag torque
aerodynamic_drag = aerodynamic_drag_torque(
    altitude, vehicle_dimensions, gravitational_constant
)
# print("Aerodynamic drag torque on the spacecraft:", aerodynamic_drag, "Nm")


def magnetic_torque(altitude):
    """
    Function to calculate the magnetic torque on the spacecraft.
    Args:
        altitude (float): Altitude of the spacecraft in kilometers.
    Returns:
        magnetic_torque (float): Magnetic torque on the spacecraft in Nm.
    """
    altitude_r = cn.earth_radius + (altitude * 10**3)  # Convert altitude to meters
    vehicle_dipole_moment = (
        0.01  # in Am^2, assuming a small dipole moment for the spacecraft
    )
    magnetic_constant = cn.magnetic_permeability_constant  # T*m^3
    vehicle_lambda = 1.1  # Dimensionless factor for magnetic torque; ranges from 1 at the equator to 2 at the poles
    magnetic_torque = (
        vehicle_dipole_moment * (magnetic_constant / altitude_r**3) * vehicle_lambda
    )  # Magnetic torque
    return magnetic_torque


# Uncomment the lines below to print the magnetic torque
magnetic_torque = magnetic_torque(altitude)
# print("Magnetic torque on the spacecraft:", magnetic_torque, "Nm")

torque_list = [
    solar_torque,
    gravity_gradient_torque,
    aerodynamic_drag,
    magnetic_torque,
]  # List of all torques
# print(torque_list)


def plot_torques(torque_list):
    """
    Function to plot the torques on the spacecraft.
    Args:
        torque_list (list): List of torques on the spacecraft.
    """
    import matplotlib.pyplot as plt

    labels = [
        "Solar Torque",
        "Gravity Gradient Torque",
        "Aerodynamic Drag",
        "Magnetic Torque",
    ]
    fig, ax = plt.subplots(figsize=(10, 6))
    ax.bar(labels, torque_list, color=["orange", "blue", "green", "red"])
    ax.set_ylabel("Torque (Nm)")
    ax.set_title("Torques on the Spacecraft")
    ax.set_ylim(0, max(torque_list) * 1.2)  # Set y-axis limit for better visibility
    ax.grid(axis="y", linestyle="--", alpha=0.7)
    plt.xticks(rotation=45)
    plt.tight_layout()
    plt.show()


# plot_torques(torque_list)  # Plot the torques on the spacecraft

disturbance_load = np.sum(
    [magnetic_torque, gravity_gradient_torque, aerodynamic_drag, solar_torque]
)  # Maximum disturbance load

# print("Total disturbance torque on the spacecraft:", disturbance_load, "Nm")
# print("Torque required to counteract disturbances:", required_torque, "Nm")

# slew rate = 3 deg/s
slew_rate = np.deg2rad(3)  # in rad/s
# average slew rate for low earth orbit spacecraft with similar maneuverability requirements

burn_time = 120  # in seconds, Max burn time of the Nammo thrusters

ang_acc_max = slew_rate / burn_time  # in rad/s^2


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
    thruster_rotation_factor = 3  # Nammo thruster is capable of rotation of 90 degrees in the x, y, and z directions

    torque_produced = (
        np.array(
            [
                MMOI_vehicle[0] * ang_acc_max,
                MMOI_vehicle[1] * ang_acc_max,
                MMOI_vehicle[2] * ang_acc_max,
            ]
        )
        + disturbance_load
    )  # Torque produced by each thruster in Nm

    # thrust_x_direction_one = torque_produced[0] / (2 * thrusters['alphard_20']['moment_arm']) + disturbance_load  # Thrust required in the x-direction
    thrust_x_direction_two = (torque_produced[0] + re_entry_moment / 3) / (
        2 * thrusters["nammo_220"]["moment_arm"]
    )  # Thrust required in the x-direction
    thrust_x_direction_three = (torque_produced[0] + re_entry_moment / 3) / (
        2 * thrusters["nammo_220_3"]["moment_arm"]
    )  # Thrust required in the x-direction
    thrust_x_direction_four = (torque_produced[0] + re_entry_moment / 3) / (
        2 * thrusters["nammo_220_4"]["moment_arm"]
    )  # Thrust required in the x-direction
    # thrust_y_direction_one = torque_produced[1] / (2 * thrusters['alphard_20']['moment_arm']) + disturbance_load  # Thrust required in the y-direction
    thrust_y_direction_two = (torque_produced[1] + re_entry_moment / 3) / (
        2 * thrusters["nammo_220"]["moment_arm"]
    )  # Thrust required in the y-direction
    thrust_y_direction_three = (torque_produced[1] + re_entry_moment / 3) / (
        2 * thrusters["nammo_220_3"]["moment_arm"]
    )  # Thrust required in the y-direction
    thrust_y_direction_four = (torque_produced[1] + re_entry_moment / 3) / (
        2 * thrusters["nammo_220_4"]["moment_arm"]
    )  # Thrust required in the y-direction
    # thrust_z_direction_one = torque_produced[2] / (2 * thrusters['alphard_20']['moment_arm']) + disturbance_load  # Thrust required in the z-direction
    thrust_z_direction_two = (torque_produced[2] + re_entry_moment / 3) / (
        2 * thrusters["nammo_220"]["moment_arm"]
    )  # Calculate the number of thrusters required for each type
    thrust_z_direction_three = (torque_produced[2] + re_entry_moment / 3) / (
        2 * thrusters["nammo_220_3"]["moment_arm"]
    )  # Thrust required in the z-direction
    thrust_z_direction_four = (torque_produced[2] + re_entry_moment / 3) / (
        2 * thrusters["nammo_220_4"]["moment_arm"]
    )  # Thrust required in the z-direction
    # number_of_thrusters['alphard_20'] = {
    #     'x': np.ceil((thrust_x_direction_two / thrusters['alphard_20']['thrust']) * redundancy_factor),
    #     'y': np.ceil((thrust_y_direction_two / thrusters['alphard_20']['thrust']) * redundancy_factor),
    #     'z': np.ceil((thrust_z_direction_two / thrusters['alphard_20']['thrust']) * redundancy_factor)
    # }

    number_of_thrusters["nammo_220"] = {
        "x": np.ceil(
            (
                thrust_x_direction_two
                / (thruster_rotation_factor * thrusters["nammo_220"]["thrust"])
            )
            * redundancy_factor
        ),
        "y": np.ceil(
            (
                thrust_y_direction_two
                / (thruster_rotation_factor * thrusters["nammo_220"]["thrust"])
            )
            * redundancy_factor
        ),
        "z": np.ceil(
            (
                thrust_z_direction_two
                / (thruster_rotation_factor * thrusters["nammo_220"]["thrust"])
            )
            * redundancy_factor
        ),
    }

    number_of_thrusters["nammo_220_3"] = {
        "x": np.ceil(
            (
                thrust_x_direction_three
                / (thruster_rotation_factor * thrusters["nammo_220_3"]["thrust"])
                * redundancy_factor
            )
        ),
        "y": np.ceil(
            (
                thrust_y_direction_three
                / (thruster_rotation_factor * thrusters["nammo_220_3"]["thrust"])
            )
            * redundancy_factor
        ),
        "z": np.ceil(
            (
                thrust_z_direction_three
                / (thruster_rotation_factor * thrusters["nammo_220_3"]["thrust"])
            )
            * redundancy_factor
        ),
    }

    number_of_thrusters["nammo_220_4"] = {
        "x": np.ceil(
            (
                thrust_x_direction_four
                / (thruster_rotation_factor * thrusters["nammo_220_4"]["thrust"])
            )
            * redundancy_factor
        ),
        "y": np.ceil(
            (
                thrust_y_direction_four
                / (thruster_rotation_factor * thrusters["nammo_220_4"]["thrust"])
            )
            * redundancy_factor
        ),
        "z": np.ceil(
            (
                thrust_z_direction_four
                / (thruster_rotation_factor * thrusters["nammo_220_4"]["thrust"])
            )
            * redundancy_factor
        ),
    }

    # Calculate the moment arm for each thruster type
    # thruster_moment_arm['alphard_20'] = {
    #     'x': torque_produced[0] / (number_of_thrusters['alphard_20']['x'] * thrusters['alphard_20']['thrust']),
    #     'y': torque_produced[1] / (number_of_thrusters['alphard_20']['y'] * thrusters['alphard_20']['thrust']),
    #     'z': torque_produced[2] / (number_of_thrusters['alphard_20']['z'] * thrusters['alphard_20']['thrust'])
    # }

    thruster_moment_arm["nammo_220"] = {
        "x": torque_produced[0]
        / (number_of_thrusters["nammo_220"]["x"] * thrusters["nammo_220"]["thrust"]),
        "y": torque_produced[1]
        / (number_of_thrusters["nammo_220"]["y"] * thrusters["nammo_220"]["thrust"]),
        "z": torque_produced[2]
        / (number_of_thrusters["nammo_220"]["z"] * thrusters["nammo_220"]["thrust"]),
    }
    thruster_moment_arm["nammo_220_3"] = {
        "x": torque_produced[0]
        / (
            number_of_thrusters["nammo_220_3"]["x"] * thrusters["nammo_220_3"]["thrust"]
        ),
        "y": torque_produced[1]
        / (
            number_of_thrusters["nammo_220_3"]["y"] * thrusters["nammo_220_3"]["thrust"]
        ),
        "z": torque_produced[2]
        / (
            number_of_thrusters["nammo_220_3"]["z"] * thrusters["nammo_220_3"]["thrust"]
        ),
    }
    thruster_moment_arm["nammo_220_4"] = {
        "x": torque_produced[0]
        / (
            number_of_thrusters["nammo_220_4"]["x"] * thrusters["nammo_220_4"]["thrust"]
        ),
        "y": torque_produced[1]
        / (
            number_of_thrusters["nammo_220_4"]["y"] * thrusters["nammo_220_4"]["thrust"]
        ),
        "z": torque_produced[2]
        / (
            number_of_thrusters["nammo_220_4"]["z"] * thrusters["nammo_220_4"]["thrust"]
        ),
    }
    power_thrusters = {
        # 'alphard_20': number_of_thrusters['alphard_20']['x'] * thrusters['alphard_20']['power'] + \
        #               number_of_thrusters['alphard_20']['y'] * thrusters['alphard_20']['power'] + \
        #               number_of_thrusters['alphard_20']['z'] * thrusters['alphard_20']['power'],
        "nammo_220": number_of_thrusters["nammo_220"]["x"]
        * thrusters["nammo_220"]["power"]
        + number_of_thrusters["nammo_220"]["y"] * thrusters["nammo_220"]["power"]
        + number_of_thrusters["nammo_220"]["z"] * thrusters["nammo_220"]["power"],
        "nammo_220_3": number_of_thrusters["nammo_220_3"]["x"]
        * thrusters["nammo_220_3"]["power"]
        + number_of_thrusters["nammo_220_3"]["y"] * thrusters["nammo_220_3"]["power"]
        + number_of_thrusters["nammo_220_3"]["z"] * thrusters["nammo_220_3"]["power"],
        "nammo_220_4": number_of_thrusters["nammo_220_4"]["x"]
        * thrusters["nammo_220_4"]["power"]
        + number_of_thrusters["nammo_220_4"]["y"] * thrusters["nammo_220_4"]["power"]
        + number_of_thrusters["nammo_220_4"]["z"] * thrusters["nammo_220_4"]["power"],
    }

    return number_of_thrusters, thruster_moment_arm, power_thrusters


number_of_thrusters, thruster_moment_arm, power_thrusters = thruster_sizing(
    thrusters, ang_acc_max, MMOI_vehicle
)
print("Number of thrusters required:", number_of_thrusters)
print("Moment arm for each thruster type:", thruster_moment_arm["nammo_220_4"])
print("Power required for each thruster type:", power_thrusters)

# print("Coordinates of the center of mass of H2ermes: ", COM)


def get_thruster_positions(COM, thrusters, number_of_thrusters):
    """
    Returns positions for all thrusters distributed around the vehicle geometry.
    Args:
        COM (np.array): Center of mass of the vehicle in meters
        thrusters (dict): Dictionary containing thruster specifications
        number_of_thrusters (dict): Dictionary containing number of thrusters per axis
    Returns:
        dict: Dictionary containing thruster positions and directions for each type and axis
    """
    all_positions = {}

    for thruster_type in number_of_thrusters.keys():
        all_positions[thruster_type] = {"x": [], "y": [], "z": []}
        moment_arm = thrusters[thruster_type]["moment_arm"]
        # dims = vehicle_dimensions

        # # Helper to clamp position within geometry
        # def clamp_position(pos):
        #     return np.array([
        #         np.clip(pos[0], -dims[0]/2, dims[0]/2),
        #         np.clip(pos[1], -dims[1]/2, dims[1]/2),
        #         np.clip(pos[2], -dims[2]/2, dims[2]/2)
        #     ])

        # X-axis thrusters
        n_thrusters_x = int(number_of_thrusters[thruster_type]["x"])
        if n_thrusters_x > 0:
            for i in range(n_thrusters_x):
                if i == 0:
                    sign = 1
                    y_offset = sign * moment_arm
                    z_offset = sign * moment_arm
                else:
                    # Alternate offsets for symmetry
                    sign = 1 if i % 2 == 0 else -1
                    y_offset = sign * moment_arm
                    z_offset = -sign * moment_arm
                pos = np.array([0, y_offset, z_offset])
                direction = np.array([sign, 0, 0])
                all_positions[thruster_type]["x"].append(
                    {"position": pos, "direction": direction}
                )

        # Y-axis thrusters
        n_thrusters_y = int(number_of_thrusters[thruster_type]["y"])
        if n_thrusters_y > 0:
            for i in range(n_thrusters_y):
                if i == 0:
                    sign = 1
                    x_offset = sign * moment_arm
                    z_offset = sign * moment_arm
                else:
                    # Alternate offsets for symmetry
                    sign = 1 if i % 2 == 0 else -1
                    x_offset = sign * moment_arm
                    z_offset = -sign * moment_arm
                pos = np.array([x_offset, 0, z_offset])
                direction = np.array([0, sign, 0])
                all_positions[thruster_type]["y"].append(
                    {"position": pos, "direction": direction}
                )

        # Z-axis thrusters
        n_thrusters_z = int(number_of_thrusters[thruster_type]["z"])
        if n_thrusters_z > 0:
            for i in range(n_thrusters_z):
                if i == 0:
                    sign = 1
                    x_offset = sign * moment_arm
                    y_offset = sign * moment_arm
                else:
                    # Alternate offsets for symmetry
                    sign = 1 if i % 2 == 0 else -1
                    x_offset = sign * moment_arm
                    y_offset = -sign * moment_arm
                pos = np.array([x_offset, y_offset, 0])
                direction = np.array([0, 0, sign])
                all_positions[thruster_type]["z"].append(
                    {"position": pos, "direction": direction}
                )

    return all_positions


all_thruster_positions = get_thruster_positions(COM, thrusters, number_of_thrusters)
# print("\n thruster positions for the two best options:")
# for thruster_type in all_thruster_positions:
#     print(f"\n{thruster_type}:")
#     for axis in ['x', 'y', 'z']:
#         print(f"\n  {axis}-axis thrusters:")
#         for i, thruster in enumerate(all_thruster_positions[thruster_type][axis]):
#             print(f"    Thruster {i+1}:")
#             print(f"      Position: {thruster['position']}")
#             print(f"      Direction: {thruster['direction']}")
