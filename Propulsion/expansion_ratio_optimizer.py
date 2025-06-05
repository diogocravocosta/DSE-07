import numpy as np
import matplotlib.pyplot as plt
import math

# Data
g0 = 9.80665
expansion_ratios = np.array([40, 50, 60, 70, 80, 85, 90, 95, 100, 110, 120, 130, 140])
vacuum_Isp = np.array([435.5235, 439.7381, 443.0285, 445.669, 447.8579, 448.8227, 449.7163, 450.5474, 451.3234, 452.7236, 453.9782, 455.0922, 456.0993])
sea_level_Isp = np.array([274.4114, 238.3587, 201.3354, 163.7011, 125.6168, 106.4455, 87.2033, 67.8989, 48.5397, 9.6704,-29.7567, -68.9614, -108.2737])
sigma = 0.121029372  # Structural mass ratio
structural_mass = 20642  # Structural mass in kg
delta_V_vacuum = 7264.29 # Delta V in vacuum in m/s
delta_V_sea_level = 250  # Delta V in m/s for landing burn
structural_mass = 21 #assuming structural mass is constant throughout the mission
payload_mass = 15000
#REVISE THIS VALUE


def optimize_expansion_ratio(expansion_ratios, vacuum_Isp, sea_level_Isp, g0, structural_mass, delta_V_vacuum, delta_V_sea_level, payload_mass):
    # Filter out negative sea level Isp values
    mask = sea_level_Isp >= 0
    filtered_expansion_ratios = expansion_ratios[mask]
    filtered_vacuum_Isp = vacuum_Isp[mask]
    filtered_sea_level_Isp = sea_level_Isp[mask]
    
    total_prop_mass_list = []

    # Loop over filtered values and compute total propellant mass
    for i in range(len(filtered_expansion_ratios)):
        #mass ratio for vacuum burn
        mass_ratio_vacuum = math.exp(delta_V_vacuum / (filtered_vacuum_Isp[i] * g0))
        propellant_mass_vacuum = mass_ratio_vacuum * (payload_mass + structural_mass) - structural_mass - payload_mass
        print(propellant_mass_vacuum)
    #mass ratio for landing burn
        mass_ratio_landing= math.exp(delta_V_sea_level / (filtered_sea_level_Isp[i] * g0))
        propellant_mass_sea_level = structural_mass * (mass_ratio_landing - 1)
        print(propellant_mass_sea_level)
        total_mass = propellant_mass_vacuum + propellant_mass_sea_level
        total_prop_mass_list.append(total_mass)

    # Plotting
    plt.plot(filtered_expansion_ratios, total_prop_mass_list, marker='o', label='Total Propellant Mass')
    plt.xlabel('Expansion Ratio')
    plt.ylabel('Total Propellant Mass (kg)')
    plt.title('Total Propellant Mass vs Expansion Ratio')
    plt.grid(True)
    plt.legend()
    plt.tight_layout()
    plt.show()

    return filtered_expansion_ratios, total_prop_mass_list

print(optimize_expansion_ratio(expansion_ratios, vacuum_Isp, sea_level_Isp, g0, structural_mass, delta_V_vacuum, delta_V_sea_level, payload_mass))

