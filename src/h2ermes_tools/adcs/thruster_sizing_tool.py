import numpy as np
import data.constants as cn

# vehicle_mass = 60000  # kg, Mass of the launch vehicle
# rcs_dry_mass = 1.48  # kg, Dry mass of the RCS thrusters including valves
# I_sp_thrusters_nammo = 160  # s, Specific impulse of the Nammo thrusters
# vehicle_dimensions = np.array([5,2.5,15])  # Lower radius, upper radius, height in meters
# altitude = 600000  # Altitude of the spacecraft in meters (600 km)
# COM = np.array([vehicle_dimensions[0]/2, vehicle_dimensions[1]/2, vehicle_dimensions[2]/2])  # Center of mass of the vehicle in meters
# burn_time = 120  # in seconds, Max burn time of the Nammo thrusters
# MMOI_vehicle = np.array([2723220, 4117860, 4117860])  # Mass Moment of Inertia (MMOI) of the spacecraft in kg*m^2
# redundancy_factor = 2  # Redundancy factor for thrusters
# solar_surface_area = 25 #input("Enter the surface area of the spacecraft exposed to sunlight in m^2: ")
# aerodynamic_surface_area = 50 #input("Enter the aerodynamic surface area of the spacecraft in m^2: ")
# # geo_midpoint = np.array([vehicle_dimensions[0]/2, vehicle_dimensions[1]/2, vehicle_dimensions[2]/4]) # Centroid in terms of base radius, top radius and height
# thrusters = {
#     # "nammo_220": {  # maximum sea level thrust
#     #     "thrust": 180,  # in N
#     #     "power": 50,  # in W
        
#     # },
#     # "nammo_220_3": {  # nominal vacuum thrust
#     #     "thrust": 220,  # in N
#     #     "power": 50,  # in W
        
#     # },
#     "nammo_220_4": {  # max vacuum thrust
#         "thrust": 250,  # in N
#         "power": 50,  # in W
#     }     
# }
def centroid(vehicle_dimensions):
    """
    Function to calculate the centroid of the spacecraft.
    Args:
        vehicle_dimensions (np.array): Dimensions of the spacecraft in meters.
    Returns:
        COM (np.array): Center of mass of the spacecraft in meters.
    """
    centroid = np.array([vehicle_dimensions[0]/2, vehicle_dimensions[1]/2, vehicle_dimensions[2]/4])
    return centroid
# geo_midpoint = centroid(vehicle_dimensions)  # Geometric midpoint of the spacecraft in meters


# # print(thrusters["nammo_220"]["moment_arm"], thrusters["nammo_220_3"]["moment_arm"], thrusters["nammo_220_4"]["moment_arm"])
# geo_midpoint = COM  # Geometric midpoint of the spacecraft in meters
#  # Actual value to be determined from the constants file

# # print("Mass Moment of Inertia (MMOI) of the spacecraft:", MMOI_vehicle)
# # print("MMOI in the x-direction:", MMOI_vehicle[0], "kg*m^2")


def solar_radiation_pressure_torque(geo_midpoint):
    """
    Function to calculate the solar radiation pressure torque on the spacecraft.
    Args:
        geo_midpoint (np.array): Centroid of the spacecraft in meters.
    Returns:
        solar_torque (float): Solar radiation pressure torque on the spacecraft in Nm.
    """

    solar_area = solar_surface_area  # m^2, surface area of the spacecraft
      # Area exposed to the sun
    q_solar = np.average(
        [0.15, 0.8]
    )  # Reflectance factor for solar radiation pressure torque
    phi_2 = 60  # in degrees

    solar_torque = (
        float(cn.solar_constant/ cn.speed_of_light)
        * solar_area
        * (1 + q_solar)
        * (geo_midpoint - COM)
        * np.cos(np.deg2rad(phi_2))
    )  # Solar radiation pressure torque

    solar_torque = np.linalg.norm(solar_torque)  # Convert to magnitude

    return solar_torque


# # # Uncomment the lines below to print the solar radiation pressure torque
# solar_torque = solar_radiation_pressure_torque(geo_midpoint)
# print("Solar radiation pressure torque on the spacecraft:", solar_torque, "Nm")


def gravity_gradient_torque(altitude, MMOI_vehicle):
    """
    Function to calculate the gravity gradient torque on the spacecraft.
    args:
        altitude (float): Altitude of the spacecraft in kilometers.
        gravitational_constant (float): Gravitational constant in m^3/s^2.
        MMOI_vehicle (np.array): Mass moment of inertia of the spacecraft in kg*m^2
    Returns:
        gravity_gradient_torque (float): Gravity gradient torque on the spacecraft in Nm.
    """
    theta = 50  # Angle between the spacecraft and the local vertical in degrees
    altitude_r = (altitude + cn.earth_radius)  # Convert altitude to meters
    theta_r = np.radians(theta)  # Convert angle to radians
    gravity_gradient_torque = (
        1.5
        * (cn.gravitational_constant / altitude_r**3)
        * np.abs(MMOI_vehicle[0] - MMOI_vehicle[1])
        * np.sin(2 * theta_r)
    )  # Gravity gradient torque
    return gravity_gradient_torque


# # # Uncomment the lines below to print the gravity gradient torque
# gravity_gradient = gravity_gradient_torque(
#     altitude, MMOI_vehicle
# )
# print("Gravity gradient torque on the spacecraft:", gravity_gradient_torque, "Nm")


def aerodynamic_drag_torque(altitude, geo_midpoint, COM):
    """Function to calculate the aerodynamic drag on the spacecraft.
    Args:
        altitude (float): Altitude of the spacecraft in meters.
        geo_midpoint (np.array): Geometric midpoint of the spacecraft in meters.
    Returns:
        aerodynamic_drag (float): Torque due to aerodynamic drag on the spacecraft in Nm.
    """
    altitude_r = (altitude + cn.earth_radius)  # Convert altitude to meters
    vehicle_velocity_squared = cn.gravitational_constant / altitude_r  # in m/s

    drag_coefficient = 2.2  # Dimensionless; typical value for spacecraft

    aerodynamic_drag = (
        0.5
        * cn.orbit_altitude_density
        * drag_coefficient
        * aerodynamic_surface_area
        * vehicle_velocity_squared
        * (np.linalg.norm(COM - geo_midpoint))
    )  # Aerodynamic drag torque
    return aerodynamic_drag


# # # Uncomment the lines below to print the aerodynamic drag torque
# aerodynamic_drag = aerodynamic_drag_torque(
#      altitude, geo_midpoint, COM
#  )
# print("Aerodynamic drag torque on the spacecraft:", aerodynamic_drag, "Nm")


def magnetic_torque(altitude):
    """
    Function to calculate the magnetic torque on the spacecraft.
    Args:
        altitude (float): Altitude of the spacecraft in meters.
        radius_earth (float): Radius of the Earth in meters.
        magnetic_constant (float): Magnetic constant in T*m^3.
    Returns:
        magnetic_torque (float): Magnetic torque on the spacecraft in Nm.
    """
    altitude_r = (altitude + cn.earth_radius)  # Convert altitude to meters
    vehicle_dipole_moment = (
        0.01  # in Am^2, assuming a small dipole moment for the spacecraft
    )

    vehicle_lambda = 1.1  # Dimensionless factor for magnetic torque; ranges from 1 at the equator to 2 at the poles
    magnetic_torque = (
        vehicle_dipole_moment * (cn.magnetic_permeability_constant / altitude_r**3) * vehicle_lambda
    )  # Magnetic torque
    return magnetic_torque


# # # Uncomment the lines below to print the magnetic torque
# magnetic_t = magnetic_torque(altitude)
# print("Magnetic torque on the spacecraft:", magnetic_torque, "Nm")

# torque_list = [
#     solar_torque,
#     gravity_gradient,
#     aerodynamic_drag,
#     magnetic_t, ]
# List of all torques
# print(torque_list)



def plot_disturbance_torques(torque_list):
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




# plot_disturbance_torques(torque_list)  # Plot the disturbance torques on the spacecraft

# disturbance_load = max(torque_list)  # Maximum disturbance load

# print("Total disturbance torque on the spacecraft:", disturbance_load, "Nm")


# slew rate = 3 deg/s; Average slew rate for low earth orbit spacecraft with similar maneuverability requirements
# slew_rate = np.deg2rad(3)  # in rad/s
# average slew rate for low earth orbit spacecraft with similar maneuverability requirements
# maneuver_time = burn_time  # in seconds, assuming a third of the maximum burn time for maneuvering
# ang_acc_max = slew_rate / (maneuver_time)  # in rad/s^2

# moment_arm_thrusters = np.array([vehicle_dimensions[0]/2, vehicle_dimensions[2]/2, vehicle_dimensions[1]/2])





def thruster_sizing(thrusters, ang_acc_max, MMOI_vehicle, torque_list, moment_arm_thrusters):
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
#     thruster_rotation_factor = 3  # Nammo thruster is capable of rotation of 90 degrees in the x, y, and z directions
    disturbance_load = max(torque_list)  # Maximum disturbance load
    torque_produced = (
        np.array(
            [
                MMOI_vehicle[0] * ang_acc_max,
                MMOI_vehicle[1] * ang_acc_max,
                MMOI_vehicle[2] * ang_acc_max,
            ]
        )
        + disturbance_load
    ) 
#     # print(torque_produced)
    # Number of thrusters required for the nammo_220 engine with nominal sea level thrust
    thrust_x_direction_two = (torque_produced[0]) / (
        moment_arm_thrusters[2]/2
        )  # Thrust required in the x-direction
    thrust_y_direction_two = (torque_produced[1]) / (
        moment_arm_thrusters[0]/2
        )  # Thrust required in the y-direction
    thrust_z_direction_two = (torque_produced[2]) / (
        moment_arm_thrusters[0]/2
        )  # Thrust required in the z-direction
    # print("Thrust in x, y, and z directions for nammo_220 are: ", thrust_x_direction_two, thrust_y_direction_two, thrust_z_direction_two, "N")
    
    
    # Number of thrusters required for the nammo_220_3 engine with nominal vacuum thrust    
    thrust_x_direction_three = (torque_produced[0]) / (
        moment_arm_thrusters[2]/2
    )  # Thrust required in the x-direction
    thrust_y_direction_three = (torque_produced[1]) / (
        moment_arm_thrusters[0]/2
    )  # Thrust required in the y-direction
    thrust_z_direction_three = (torque_produced[2]) / (
        moment_arm_thrusters[0]/2
    )  # Thrust required in the z-direction
    # print("Thrust in x, y, and z directions for nammo_220_3 are: ", thrust_x_direction_three, thrust_y_direction_three, thrust_z_direction_three, "N")
    
    # Number of thrusters required for the nammo_220_4 engine with maximum vacuum thrust
    thrust_x_direction_four = (torque_produced[0]) / (
        moment_arm_thrusters[2]/2
    )  # Thrust required in the x-direction
    thrust_y_direction_four = (torque_produced[1]) / (
        moment_arm_thrusters[0]/2
    )  # Thrust required in the y-direction
    thrust_z_direction_four = (torque_produced[2]) / (
        moment_arm_thrusters[0]/2
    )  # Thrust required in the z-direction
    # print("Thrust in x, y, and z directions for nammo_220_4 are: ", thrust_x_direction_four, thrust_y_direction_four, thrust_z_direction_four, "N")
    # print("Maximum thrust required: ", max(thrust_x_direction_four, thrust_y_direction_four, thrust_z_direction_four), "N")
#     

   
    # number_of_thrusters["nammo_220"] = {
    #     "x": np.ceil(
    #         (
    #             thrust_x_direction_two
    #             / (thrusters["nammo_220"]["thrust"])
    #         )
    #         * redundancy_factor
    #     ),
    #     "y": np.ceil(
    #         (
    #             thrust_y_direction_two
    #             / (thrusters["nammo_220"]["thrust"])
    #         )
    #         * redundancy_factor
    #     ),
    #     "z": np.ceil(
    #         (
    #             thrust_z_direction_two
    #             / (thrusters["nammo_220"]["thrust"])
    #         )
    #         * redundancy_factor
    #     ),
    # }

    # number_of_thrusters["nammo_220_3"] = {
    #     "x": np.ceil(
    #         (
    #             thrust_x_direction_three
    #             / (thrusters["nammo_220_3"]["thrust"])
    #             * redundancy_factor
    #         )
    #     ),
    #     "y": np.ceil(
    #         (
    #             thrust_y_direction_three
    #             / (thrusters["nammo_220_3"]["thrust"])
    #         )
    #         * redundancy_factor
    #     ),
    #     "z": np.ceil(
    #         (
    #             thrust_z_direction_three
    #             / (thrusters["nammo_220_3"]["thrust"])
    #         )
    #         * redundancy_factor
    #     ),
    # }

    number_of_thrusters["nammo_220_4"] = {
        "x": np.ceil(
            (
                thrust_x_direction_four
                / (thrusters["nammo_220_4"]["thrust"])
            )
            * redundancy_factor
        ),
        "y": np.ceil(
            (
                thrust_y_direction_four
                / (thrusters["nammo_220_4"]["thrust"])
            )
            * redundancy_factor
        ),
        "z": np.ceil(
            (
                thrust_z_direction_four
                / (thrusters["nammo_220_4"]["thrust"])
            )
            * redundancy_factor
        ),
    }

    return number_of_thrusters


# number_of_thrusters = thruster_sizing(
#     thrusters, ang_acc_max, MMOI_vehicle
# )
# print("Number of thrusters required:", number_of_thrusters["nammo_220_4"])

#### IGNORE this
# def get_thruster_positions(COM, thrusters, number_of_thrusters):
#     """
#     Returns positions for all thrusters distributed around the vehicle geometry.
#     Args:
#         COM (np.array): Center of mass of the vehicle in meters
#         thrusters (dict): Dictionary containing thruster specifications
#         number_of_thrusters (dict): Dictionary containing number of thrusters per axis
#     Returns:
#         dict: Dictionary containing thruster positions and directions for each type and axis
#     """
#     all_positions = {}

     
        
        

#         moment_arm = thrusters[thruster_type]["moment_arm"]
#         # X-axis thrusters
#         n_thrusters_x = int(number_of_thrusters[thruster_type]["x"])         
#         if n_thrusters_x > 0:
#             for i in range(n_thrusters_x):
#                 if i == 0:
#                     sign = 1
#                     y_offset = sign * moment_arm
#                     z_offset = - sign * moment_arm
#                 else:
#                     # Alternate offsets for symmetry
#                     sign = 1 if i % 2 == 0 else -1
#                     y_offset = sign * moment_arm
#                     z_offset = sign * moment_arm
#                 pos = np.array([0, y_offset, z_offset])
#                 direction = np.array([sign, 0, 0])
#                 all_positions[thruster_type]["x"].append(
#                     {"position": pos, "direction": direction}
#                 )

#         # Y-axis thrusters
#         n_thrusters_y = int(number_of_thrusters[thruster_type]["y"])
#         if n_thrusters_y > 0:
#             for i in range(n_thrusters_y):
#                 if i == 0:
#                     sign = 1
#                     x_offset = - sign * moment_arm
#                     z_offset = sign * moment_arm
#                 else:
#                     # Alternate offsets for symmetry
#                     sign = 1 if i % 2 == 0 else -1
#                     x_offset = sign * moment_arm
#                     z_offset = sign * moment_arm
#                 pos = np.array([x_offset, 0, z_offset])
#                 direction = np.array([0, sign, 0])
#                 all_positions[thruster_type]["y"].append(
#                     {"position": pos, "direction": direction}
#                 )

#         # Z-axis thrusters
#         n_thrusters_z = int(number_of_thrusters[thruster_type]["z"])
#         if n_thrusters_z > 0:
#             for i in range(n_thrusters_z):
#                 if i == 0:
#                     sign = 1
#                     x_offset = sign * moment_arm
#                     y_offset = - sign * moment_arm
#                 else:
#                     # Alternate offsets for symmetry
#                     sign = 1 if i % 2 == 0 else -1
#                     x_offset = sign * moment_arm
#                     y_offset =  sign * moment_arm
#                 pos = np.array([x_offset, y_offset, 0])
#                 direction = np.array([0, 0, sign])
#                 all_positions[thruster_type]["z"].append(
#                     {"position": pos, "direction": direction}
#                 )


#     return all_positions


# all_thruster_positions = get_thruster_positions(COM, thrusters, number_of_thrusters)


# # print("\nAll thruster positions:")
# # for thruster_type in all_thruster_positions:
# #     print(f"\n{thruster_type}:")
# #     for axis in ['x', 'y', 'z']:
# #         print(f"\n  {axis}-axis thrusters:")
# #         for i, thruster in enumerate(all_thruster_positions[thruster_type][axis]):
# #             print(f"    Thruster {i+1}:")
# #             print(f"      Position: {thruster['position']}")
# #             print(f"      Direction: {thruster['direction']}")






def mass_and_power_estimation(number_of_thrusters, rcs_dry_mass, thrusters, burn_time):
    """
    Function to compute total mass and power consumption of the thrusters.
    Args:
        number_of_thrusters (dict): Dictionary containing updated thruster positions and directions.
        thrusters (dict): Dictionary containing thruster specifications.
        burn_time (int): Burn time for the thrusters in seconds.
    Returns:
        total_power(float): Total power consumption of the thrusters in watts.
        
   """
    #Calculate total dry mass for each thrust level
    dry_mass = { 
        
        # "nammo_220": (
        #         number_of_thrusters["nammo_220"]["x"] +
        #          number_of_thrusters["nammo_220"]["y"] +
        #          number_of_thrusters["nammo_220"]["z"]
        #  ) * rcs_dry_mass,  # in kg
        # "nammo_220_3": (
        #         number_of_thrusters["nammo_220_3"]["x"] +   
        #         number_of_thrusters["nammo_220_3"]["y"] +
        #         number_of_thrusters["nammo_220_3"]["z"]
        #  ) * rcs_dry_mass,  # in kg
        "nammo_220_4": (
                number_of_thrusters["nammo_220_4"]["x"] +
                number_of_thrusters["nammo_220_4"]["y"] +
                number_of_thrusters["nammo_220_4"]["z"]
         ) * rcs_dry_mass  # in kg

    }

    
        
    #Calculate total propellant mass for each thrust level
    # thrust_nammo_220 = thrusters["nammo_220"]["thrust"]
    # thrust_nammo_220_3 = thrusters["nammo_220_3"]["thrust"]
    thrust_nammo_220_4 = thrusters["nammo_220_4"]["thrust"]
    # Propellant mass for each thruster type
    propellant_mass = {
    #     "nammo_220": (number_of_thrusters["nammo_220"][ "x"] + 
    #                       number_of_thrusters["nammo_220"]["y"] + 
    #                       number_of_thrusters["nammo_220"]["z"]) * ((thrust_nammo_220 * burn_time) / (I_sp_thrusters_nammo * g_0)),  # in kg

    #    "nammo_220_3":  ( number_of_thrusters["nammo_220_3"][ "x"] + 
    #                       number_of_thrusters["nammo_220_3"]["y"] + 
    #                       number_of_thrusters["nammo_220_3"]["z"]) * ((thrust_nammo_220_3 * burn_time) / (I_sp_thrusters_nammo * g_0)),  # in kg

        "nammo_220_4": ( number_of_thrusters["nammo_220_4"][ "x"] + 
                          number_of_thrusters["nammo_220_4"]["y"] + 
                          number_of_thrusters["nammo_220_4"]["z"]) * ((thrust_nammo_220_4 * burn_time) / (I_sp_thrusters_nammo * cn.g_0))  # in kg

    }


    # Mass calculations for thrusters
    total_propellant_mass = propellant_mass["nammo_220_4"]  
    # print(total_propellant_mass) # Total propellant mass in kg
    total_dry_mass = dry_mass["nammo_220_4"]  # Total dry mass in kg
    total_mass = total_dry_mass + total_propellant_mass  # Total mass in kg
    #Calculate total power consumption of the thrusters
    power_thrusters = {
    # "nammo_220": ( number_of_thrusters["nammo_220"][ "x"] +
    #                       number_of_thrusters["nammo_220"]["y"] +
    #                       number_of_thrusters["nammo_220"]["z"]) * thrusters["nammo_220"]["power"],
                          
    # "nammo_220_3": ( number_of_thrusters["nammo_220_3"][ "x"] + 
    #                       number_of_thrusters["nammo_220_3"]["y"] + 
    #                       number_of_thrusters["nammo_220_3"]["z"]) * thrusters["nammo_220_3"]["power"],

    "nammo_220_4": ( number_of_thrusters["nammo_220_4"][ "x"] + 
                          number_of_thrusters["nammo_220_4"]["y"] + 
                          number_of_thrusters["nammo_220_4"]["z"]) * thrusters["nammo_220_4"]["power"],
    }
    total_power = power_thrusters["nammo_220_4"]  # Peak power requirement in W
    return total_power, total_mass, total_dry_mass, total_propellant_mass
# total_power_thrusters, total_mass_thrusters, total_dry_mass_thrusters, total_prop_mass_thrusters = mass_and_power_estimation(number_of_thrusters, rcs_dry_mass, thrusters, burn_time)
# print("The total power requirement of the RCS thrusters is: ", total_power_thrusters, "W")
# print("The total dry mass of the RCS thrusters is: ", total_dry_mass_thrusters, "kg")
# print("The total propellant mass of the RCS thrusters is: ", total_prop_mass_thrusters, "kg")
# print("The total mass of the RCS thrusters is: ", total_mass_thrusters, "kg")

# htp_density = 1134.5 #kg/m^3, limiting density at highest inlet temperature of 80 degrees celsius

def acs_tank_design(htp_density, total_prop_mass):
    tank_pressure = 5e6  # in Pa
    burst_safety_factor = 2 # Safety factor for burst pressure
    tank_pressure *= burst_safety_factor  # Adjusting the tank pressure for safety factor
    Materials = {
                "Stainless Steel 316":
                    { 'yield_strength': 205*10**6, # in Pa    
                    'density': 7870 # kg/m**3
                },

                "Aluminum AA6000 T6 Series": 
                    { 'yield_strength': 276 * 10**6, # in Pa
                     'density': 2700 # kg/m**3}
    }
    }
    

    ullage_factor = 1.1  # Ullage factor to account for propellant expansion and sloshing
    tank_volume = total_prop_mass / htp_density  # in m^3
    # print(tank_volume)
    new_tank_volume = tank_volume * ullage_factor  # Adjusting the tank volume for ullage factor
    # print("The volume of the ACS tank is: ", new_tank_volume, "m^3")
    tank_radius = ((tank_volume / np.pi) * 0.75) ** (1/3)  # Assuming a spherical tank for simplicity
    # print(tank_radius)
    tank_thickness_SS_316 = (tank_pressure * tank_radius)/(2 * Materials['Stainless Steel 316']['yield_strength'])  # Using the formula for thin-walled pressure vessels
    tank_thickness_AA_6000 = (tank_pressure * tank_radius)/(2 * Materials['Aluminum AA6000 T6 Series']['yield_strength'])  # Using the formula for thin-walled pressure vessels
    # print("aluminum thickness", tank_thickness_AA_6000)
    tank_mass_SS_316 = 4 * np.pi * tank_radius**2 * tank_thickness_SS_316 * Materials['Stainless Steel 316']['density']  # Mass of the tank in kg
    tank_mass_AA_6000 = 4 * np.pi * tank_radius**2 * tank_thickness_AA_6000 * Materials['Aluminum AA6000 T6 Series']['density']  # Mass of the tank in kg
    tank_mass = min(tank_mass_SS_316, tank_mass_AA_6000)  # Taking the minimum mass of the tank
    tank_material = 'Aluminum AA6000 T6 Series' if tank_mass_AA_6000 < tank_mass_SS_316 else "Stainless Steel 316"  # Choosing the material with the lower mass
    return tank_mass, tank_material

# tank_mass, tank_material = acs_tank_design(htp_density, total_prop_mass_thrusters)
# print("The mass of the ACS tank is: ", tank_mass, "kg")
# print("The material of the ACS tank is: ", tank_material)

# print("The total mass of the ACS system is: ", total_mass_thrusters + tank_mass, "kg")
# print("The total power requirement of the ACS system is: ", total_power_thrusters, "W")

def aocs_mass_function(MMOI_vehicle, COM, vehicle_dimensions, htp_density, slew_rate, burn_time, moment_arm_thrusters):
    """
    Function to compute the total mass of the AOCS for integration 
    """

    ang_acc_max = slew_rate / (burn_time)  # in rad/s^2, Maximum angular acceleration required for the spacecraft
    #Calculate the geometric midpoint of the spacecraft
    geo_midpoint = centroid(vehicle_dimensions)  # Geometric midpoint of the spacecraft in meters
    # Calculate the disturbance torques
    solar_torque = solar_radiation_pressure_torque(geo_midpoint)
    gravity_gradient = gravity_gradient_torque(altitude, MMOI_vehicle)
    aerodynamic_drag = aerodynamic_drag_torque(altitude, geo_midpoint, COM)
    magnetic_t = magnetic_torque(altitude)

    # List of all torques
    torque_list = [solar_torque, gravity_gradient, aerodynamic_drag, magnetic_t]
    
    # Calculate the maximum disturbance load
    disturbance_load = max(torque_list)

    # Calculate the number of thrusters required
    number_of_thrusters = thruster_sizing(thrusters, ang_acc_max, MMOI_vehicle, torque_list, moment_arm_thrusters)

    # Estimate the mass consumption of the thrusters
    _, total_mass_thrusters, _, total_prop_mass_thrusters = mass_and_power_estimation(number_of_thrusters, rcs_dry_mass, thrusters, burn_time)

    # Design the ACS tank
    tank_mass, _ = acs_tank_design(htp_density, total_prop_mass_thrusters)

    # Calculate the total mass of the ACS system
    total_mass_acs_system = total_mass_thrusters + tank_mass

    return total_mass_acs_system


# total_mass_aocs = aocs_mass_function(MMOI_vehicle, COM, vehicle_dimensions, htp_density, slew_rate)
# print("The total mass of the Attitude and Orbit Control System (AOCS) is: ", total_mass_aocs, "kg")



if __name__ == "__main__":
    vehicle_mass = 60000  # kg, Mass of the launch vehicle
    rcs_dry_mass = 1.48  # kg, Dry mass of the RCS thrusters including valves
    I_sp_thrusters_nammo = 160  # s, Specific impulse of the Nammo thrusters
    vehicle_dimensions = np.array([5,2.5,15])  # Lower radius, upper radius, height in meters
    altitude = 600000  # Altitude of the spacecraft in meters (600 km)
    COM = np.array([vehicle_dimensions[0]/2, vehicle_dimensions[1]/2, vehicle_dimensions[2]/2])  # Center of mass of the vehicle in meters
    burn_time = 120  # in seconds, Max burn time of the Nammo thrusters
    MMOI_vehicle = np.array([2723220, 4117860, 4117860])  # Mass Moment of Inertia (MMOI) of the spacecraft in kg*m^2
    redundancy_factor = 2  # Redundancy factor for thrusters
    solar_surface_area = 25 #input("Enter the surface area of the spacecraft exposed to sunlight in m^2: ")
    aerodynamic_surface_area = 50 #input("Enter the aerodynamic surface area of the spacecraft in m^2: ")
    htp_density = 1134.5 #kg/m^3, limiting density at highest inlet temperature of 80 degrees celsius
    moment_arm_thrusters = np.array([vehicle_dimensions[0]/2, vehicle_dimensions[2]/2, vehicle_dimensions[1]/2])  # Moment arm for the thrusters in meters
      # in rad/s^2, Average slew rate for low earth orbit spacecraft with similar maneuverability requirements
    thrusters = {
        # "nammo_220": {  # maximum sea level thrust
        #     "thrust": 180,  # in N
        #     "power": 50,  # in W
            
        # },
        # "nammo_220_3": {  # nominal vacuum thrust
        #     "thrust": 220,  # in N
        #     "power": 50,  # in W
            
        # },
        "nammo_220_4": {  # max vacuum thrust
            "thrust": 250,  # in N
            "power": 50,  # in W
        }     
    }
    # slew rate = 3 deg/s; Average slew rate for low earth orbit spacecraft with similar maneuverability requirements
    slew_rate = np.deg2rad(3)  # in rad/s, average slew rate for low earth orbit spacecraft with similar maneuverability requirements
    # maneuver_time = burn_time  # in seconds, assuming a third of the maximum burn time for maneuvering
    # ang_acc_max = slew_rate / (maneuver_time)  # in rad/s^2

    total_mass_aocs = aocs_mass_function(MMOI_vehicle, COM, vehicle_dimensions, htp_density, slew_rate, burn_time, moment_arm_thrusters)
    print("The total mass of the Attitude and Orbit Control System (AOCS) is: ", total_mass_aocs, "kg")


# List of all torques