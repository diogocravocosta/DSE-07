import numpy as np
import data.constants as cn
from h2ermes_tools.variables import mass_vehicle, launch_vehicle_dimensions, MMOI_vehicle, vehicle_surface_area

# Constants and parameters for the spacecraft
vehicle_mass = mass_vehicle.value  # kg
rcs_dry_mass = 1.48 # kg, Dry mass of the RCS thrusters including valves
# vehicle_shape = "rectangular_prism"
vehicle_dimensions = launch_vehicle_dimensions.value  # Length, Width, Height in meters
gravitational_constant = cn.gravitational_parameter  # m^3/s^2
theta = 50  # Angle between the spacecraft and the local vertical in degrees
altitude = 600000  # Altitude of the spacecraft in meters (600 km)
radius_earth = cn.earth_radius  # Radius of the Earth in meters
COM = np.array(
    [vehicle_dimensions[0] / 2, vehicle_dimensions[1] / 2, vehicle_dimensions[2] / 2]
)  # Center of mass of the spacecraft in meters
re_entry_moment = 166000  # Maximum moment during re-entry in Nm
redundancy_factor = 2  # Redundancy factor for thrusters
I_sp_thrusters_nammo = 160 # s, Specific impulse of the Nammo thrusters
g_0 = cn.g_0  # m/s^2, Standard 
burn_time = 120  # in seconds, Max burn time of the Nammo thrusters



# velocity = np.sqrt(gravitational_constant / ((altitude + 6371) * 1000))  # in m/s
# print("Velocity of the spacecraft at altitude", altitude, "km:", velocity, "m/s")





MMOI_vehicle = MMOI_vehicle.value  # Mass Moment of Inertia (MMOI) of the spacecraft in kg*m^2
# print("Mass Moment of Inertia (MMOI) of the spacecraft:", MMOI_vehicle)

thrusters = {
    "nammo_220": {  # sea level thrust
        "thrust": 148,  # in N
        "power": 50,  # in W
        "moment_arm": np.linalg.norm(COM),  # in meters
    },
    "nammo_220_3": {  # nominal vacuum thrust
        "thrust": 220,  # in N
        "power": 50,  # in W
        "moment_arm": np.linalg.norm(COM),  # in meters
    },
    "nammo_220_4": {  # max vacuum thrust
        "thrust": 250,  # in N
        "power": 50,  # in W
        "moment_arm": np.linalg.norm(COM),  # in meters
    },
}

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
    solar_area = vehicle_surface_area.value  # m^2, surface area of the spacecraft
      # Area exposed to the sun

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
    altitude_r = (altitude + radius_earth)  # Convert altitude to meters
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
    altitude_r = (altitude + radius_earth)  # Convert altitude to meters
    vehicle_velocity_squared = gravitational_constant / altitude_r  # in m/s

    drag_coefficient = 2.2  # Dimensionless; typical value for spacecraft
    atmospheric_density = cn.orbit_altitude_density  # kg/m^3
    aerodynamic_surface_area = vehicle_surface_area.value # Area exposed to the flow

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
    altitude_r = (altitude + radius_earth)  # Convert altitude to meters
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

ang_acc_max = slew_rate / (burn_time)  # in rad/s^2





def thruster_sizing(thrusters, ang_acc_max, MMOI_vehicle):
    """
    Function to size the thrusters based on the required torque and angular acceleration.
    Args:
        thrusters (dict): Dictionary containing thruster specifications.
        ang_acc_max (float): Maximum angular acceleration in rad/s^2.
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


    # Number of thrusters required for the nammo_220 engine with nominal sea level thrust
    thrust_x_direction_two = (torque_produced[0] + re_entry_moment / 3) / (
        2 * thrusters["nammo_220"]["moment_arm"]
    )  # Thrust required in the x-direction
    thrust_y_direction_two = (torque_produced[1] + re_entry_moment / 3) / (
        2 * thrusters["nammo_220"]["moment_arm"]
    )  # Thrust required in the y-direction
    thrust_z_direction_two = (torque_produced[2] + re_entry_moment / 3) / (
    2 * thrusters["nammo_220"]["moment_arm"]
    )  # Thrust required in the z-direction
    print("Thrust in x, y, and z directions for nammo_220 are: ", thrust_x_direction_two, thrust_y_direction_two, thrust_z_direction_two, "N")
    
    
    # Number of thrusters required for the nammo_220_3 engine with nominal vacuum thrust    
    thrust_x_direction_three = (torque_produced[0] + re_entry_moment / 3) / (
        2 * thrusters["nammo_220_3"]["moment_arm"]
    )  # Thrust required in the x-direction
    thrust_y_direction_three = (torque_produced[1] + re_entry_moment / 3) / (
        2 * thrusters["nammo_220_3"]["moment_arm"]
    )  # Thrust required in the y-direction
    thrust_z_direction_three = (torque_produced[2] + re_entry_moment / 3) / (
        2 * thrusters["nammo_220_3"]["moment_arm"]
    )  # Thrust required in the z-direction
    # print("Thrust in x, y, and z directions for nammo_220_3 are: ", thrust_x_direction_three, thrust_y_direction_three, thrust_z_direction_three, "N")

    # Number of thrusters required for the nammo_220_4 engine with maximum vacuum thrust
    thrust_x_direction_four = (torque_produced[0] + re_entry_moment / 3) / (
        2 * thrusters["nammo_220_4"]["moment_arm"]
    )  # Thrust required in the x-direction
    thrust_y_direction_four = (torque_produced[1] + re_entry_moment / 3) / (
        2 * thrusters["nammo_220_4"]["moment_arm"]
    )  # Thrust required in the y-direction
    thrust_z_direction_four = (torque_produced[2] + re_entry_moment / 3) / (
        2 * thrusters["nammo_220_4"]["moment_arm"]
    )  # Thrust required in the z-direction
    # print("Thrust in x, y, and z directions for nammo_220_4 are: ", thrust_x_direction_four, thrust_y_direction_four, thrust_z_direction_four, "N")

    # Calculate the number of thrusters required for each direction and thruster type
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

    return number_of_thrusters


number_of_thrusters = thruster_sizing(
    thrusters, ang_acc_max, MMOI_vehicle
)
# print("Number of thrusters required:", number_of_thrusters)

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

        if thruster_type == "nammo_220_3":
            
            number_of_thrusters[thruster_type]["x"] = (
                number_of_thrusters[thruster_type]["x"] // 2
            )
            number_of_thrusters[thruster_type]["y"] = (
                number_of_thrusters[thruster_type]["y"] // 2
            )
            number_of_thrusters[thruster_type]["z"] = (
                number_of_thrusters[thruster_type]["z"] // 2 
            )
            
            
        if thruster_type == "nammo_220_4":
            
            number_of_thrusters[thruster_type]["x"] = (
                number_of_thrusters[thruster_type]["x"] // 2
            )
            number_of_thrusters[thruster_type]["y"] = ( 
                number_of_thrusters[thruster_type]["y"] // 2
            )
            number_of_thrusters[thruster_type]["z"] = ( 
                number_of_thrusters[thruster_type]["z"] // 2
            )
            
        elif thruster_type == "nammo_220":
           
            number_of_thrusters[thruster_type]["x"] = (
                number_of_thrusters[thruster_type]["x"] // 2
            )
            number_of_thrusters[thruster_type]["y"] = (
                number_of_thrusters[thruster_type]["y"] // 2
            )
            number_of_thrusters[thruster_type]["z"] = (
                number_of_thrusters[thruster_type]["z"] // 2
            )
        
        
        

        moment_arm = thrusters[thruster_type]["moment_arm"]
        # X-axis thrusters
        n_thrusters_x = int(number_of_thrusters[thruster_type]["x"])         
        if n_thrusters_x > 0:
            for i in range(n_thrusters_x):
                if i == 0:
                    sign = 1
                    y_offset = sign * moment_arm
                    z_offset = - sign * moment_arm
                else:
                    # Alternate offsets for symmetry
                    sign = 1 if i % 2 == 0 else -1
                    y_offset = sign * moment_arm
                    z_offset = sign * moment_arm
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
                    x_offset = - sign * moment_arm
                    z_offset = sign * moment_arm
                else:
                    # Alternate offsets for symmetry
                    sign = 1 if i % 2 == 0 else -1
                    x_offset = sign * moment_arm
                    z_offset = sign * moment_arm
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
                    y_offset = - sign * moment_arm
                else:
                    # Alternate offsets for symmetry
                    sign = 1 if i % 2 == 0 else -1
                    x_offset = sign * moment_arm
                    y_offset =  sign * moment_arm
                pos = np.array([x_offset, y_offset, 0])
                direction = np.array([0, 0, sign])
                all_positions[thruster_type]["z"].append(
                    {"position": pos, "direction": direction}
                )
    updated_positions = number_of_thrusters


    return all_positions, updated_positions


all_thruster_positions, updated_positions = get_thruster_positions(COM, thrusters, number_of_thrusters)
print(updated_positions)

# print("\nAll thruster positions:")
# for thruster_type in all_thruster_positions:
#     print(f"\n{thruster_type}:")
#     for axis in ['x', 'y', 'z']:
#         print(f"\n  {axis}-axis thrusters:")
#         for i, thruster in enumerate(all_thruster_positions[thruster_type][axis]):
#             print(f"    Thruster {i+1}:")
#             print(f"      Position: {thruster['position']}")
#             print(f"      Direction: {thruster['direction']}")

# print("\nUpdated number of thrusters after reduction:", updated_positions)




def mass_and_power_estimation(updated_positions, thrusters, burn_time = burn_time):
    """
    Function to compute total mass and power consumption of the thrusters.
    Args:
        updated_positions (dict): Dictionary containing updated thruster positions and directions.
        thrusters (dict): Dictionary containing thruster specifications.
        burn_time (int): Burn time for the thrusters in seconds.
    Returns:
        total_power(float): Total power consumption of the thrusters in watts.
        
   """
    #Calculate total dry mass of the thrusters
    total_dry_mass = {
        "nammo_220": updated_positions["nammo_220"]["x"] * rcs_dry_mass
        + updated_positions["nammo_220"]["y"] * rcs_dry_mass
        + updated_positions["nammo_220"]["z"] * rcs_dry_mass,
        "nammo_220_3": updated_positions["nammo_220_3"]["x"] * rcs_dry_mass
        + updated_positions["nammo_220_3"]["y"] * rcs_dry_mass
        + updated_positions["nammo_220_3"]["z"] * rcs_dry_mass,
        "nammo_220_4": updated_positions["nammo_220_4"]["x"] * rcs_dry_mass
        + updated_positions["nammo_220_4"]["y"] * rcs_dry_mass
        + updated_positions["nammo_220_4"]["z"] * rcs_dry_mass,
    }
    #Calculate total propellant mass for each thrust level
    thrust_nammo_220 = thrusters["nammo_220"]["thrust"]
    thrust_nammo_220_3 = thrusters["nammo_220_3"]["thrust"]
    thrust_nammo_220_4 = thrusters["nammo_220_4"]["thrust"]
    # Propellant mass for each thruster type
    propellant_mass = {
        "nammo_220": (updated_positions["nammo_220"]["x"] * thrust_nammo_220
        + updated_positions["nammo_220"]["y"] * thrust_nammo_220
        + updated_positions["nammo_220"]["z"] * thrust_nammo_220
    ) * burn_time / (I_sp_thrusters_nammo * g_0),  # in kg

       "nammo_220_3":  (updated_positions["nammo_220_3"]["x"] * thrust_nammo_220_3
        + updated_positions["nammo_220_3"]["y"] * thrust_nammo_220_3
        + updated_positions["nammo_220_3"]["z"] * thrust_nammo_220_3
        ) * burn_time / (I_sp_thrusters_nammo * g_0),  # in kg

        "nammo_220_4": (updated_positions["nammo_220_4"]["x"] * thrust_nammo_220_4
        + updated_positions["nammo_220_4"]["y"] * thrust_nammo_220_4
        + updated_positions["nammo_220_4"]["z"] * thrust_nammo_220_4
         ) * burn_time / (I_sp_thrusters_nammo * g_0)  # in kg

    }

    # Mass calculations for thrusters
    total_propellant_mass = sum(propellant_mass.values())  
    # print(total_propellant_mass) # Total propellant mass in kg
    total_dry_mass = sum(total_dry_mass.values())  # Total dry mass in kg
    total_mass = total_dry_mass + total_propellant_mass  # Total mass in kg
    #Calculate total power consumption of the thrusters
    power_thrusters = {
    "nammo_220": updated_positions["nammo_220"]["x"]
    * thrusters["nammo_220"]["power"]
    + updated_positions["nammo_220"]["y"] * thrusters["nammo_220"]["power"]
    + updated_positions["nammo_220"]["z"] * thrusters["nammo_220"]["power"],
    "nammo_220_3": updated_positions["nammo_220_3"]["x"]
    * thrusters["nammo_220_3"]["power"]
    + updated_positions["nammo_220_3"]["y"] * thrusters["nammo_220_3"]["power"]
    + updated_positions["nammo_220_3"]["z"] * thrusters["nammo_220_3"]["power"],
    "nammo_220_4": updated_positions["nammo_220_4"]["x"]
    * thrusters["nammo_220_4"]["power"]
    + updated_positions["nammo_220_4"]["y"] * thrusters["nammo_220_4"]["power"]
    + updated_positions["nammo_220_4"]["z"] * thrusters["nammo_220_4"]["power"],
    }
    total_power = sum(power_thrusters.values())  # Total power consumption in W
    return total_power, total_mass, total_dry_mass, total_propellant_mass

total_power_thrusters, total_mass_thrusters, total_dry_mass_thrusters, total_prop_mass_thrusters = mass_and_power_estimation(updated_positions, thrusters)
print("The total power requirement of the RCS thrusters is: ", total_power_thrusters, "W")
# print("The total dry mass of the RCS thrusters is: ", total_dry_mass_thrusters, "kg")
# print("The total propellant mass of the RCS thrusters is: ", total_prop_mass_thrusters, "kg")
print("The total mass of the RCS thrusters is: ", total_mass_thrusters, "kg")

htp_density = 1134.5 #kg/m^3, limiting density at highest inlet temperature of 80 degrees celsius

def acs_tank_design(htp_density, total_prop_mass):
    tank_pressure = 8e6  # in Pa
    Materials = {
        "Stainless Steel SS304L": 
                     { 'yield_strength': 215 * 10**6, # in Pa
                      'density': 7930 # kg/m**3 
        },

                "Titanium Ti6Al4V Alloy": 
                    { 'yield_strength': 310 * 10**6, # in Pa
                    'density': 4430 # kg/m**3
                }  
    }
    
    total_prop_mass = total_prop_mass_thrusters
    tank_volume = total_prop_mass / htp_density  # in m^3
    tank_radius = ((tank_volume / np.pi) * 0.75) ** (1/3)  # Assuming a spherical tank for simplicity
    tank_thickness_SS = (tank_pressure * tank_radius)/(Materials['Stainless Steel SS304L']['yield_strength'])  # Using the formula for thin-walled pressure vessels
    tank_thickness_Ti = (tank_pressure * tank_radius)/(Materials["Titanium Ti6Al4V Alloy"]['yield_strength'])  # Using the formula for thin-walled pressure vessels
    print(tank_thickness_SS, tank_thickness_Ti)
    tank_mass_SS = 4 * np.pi * tank_radius**2 * tank_thickness_SS * Materials['Stainless Steel SS304L']['density']  # Mass of the tank in kg
    tank_mass_Ti = 4 * np.pi * tank_radius**2 * tank_thickness_Ti * Materials["Titanium Ti6Al4V Alloy"]['density']  # Mass of the tank in kg
    tank_mass = min(tank_mass_SS, tank_mass_Ti)  # Taking the maximum mass of the tank
    tank_material = "Titanium Ti6Al4V Alloy" if tank_mass_Ti < tank_mass_SS else "Stainless Steel SS304L"
    return tank_mass, tank_material

tank_mass, tank_material = acs_tank_design(htp_density, total_prop_mass_thrusters)
# print("The mass of the ACS tank is: ", tank_mass, "kg")
# print("The material of the ACS tank is: ", tank_material)
    

