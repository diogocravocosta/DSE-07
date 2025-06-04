import numpy as np
from Mixture_Ratio_Optimizer import calculate_propmass
import matplotlib.pyplot as plt

# Data
expansion_ratios = np.linspace(40, 140, 11)
vacuum_Isp = np.array([435.5235, 439.7381, 443.0285, 445.669, 447.8579, 466.1641, 451.3234, 452.7236, 453.9782, 455.0922, 456.0993])
sea_level_Isp = np.array([274.4114, 238.3587, 201.3354, 163.7011, 125.6168, 103.6511, 48.5397, 9.6704, -29.7567, -68.9614, -108.2737])

def optimize_expansion_ratio(expansion_ratios, vacuum_Isp, sea_level_Isp):
    # Filter out negative sea level Isp values
    mask = sea_level_Isp >= 0
    filtered_expansion_ratios = expansion_ratios[mask]
    filtered_vacuum_Isp = vacuum_Isp[mask]
    filtered_sea_level_Isp = sea_level_Isp[mask]

    total_prop_mass_list = []

    # Loop over filtered values and compute total propellant mass
    for i in range(len(filtered_expansion_ratios)):
        propellant_mass_vacuum = calculate_propmass(
            struct_ratio=0.121029372268477, 
            Isp=filtered_vacuum_Isp[i], 
            deltaV=7264.29, 
            payload=15000
        )
        propellant_mass_sea_level = calculate_propmass(
            struct_ratio=0.121029372268477, 
            Isp=filtered_sea_level_Isp[i], 
            deltaV=200, 
            payload=15000
        )
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

print(optimize_expansion_ratio(expansion_ratios, vacuum_Isp, sea_level_Isp))

