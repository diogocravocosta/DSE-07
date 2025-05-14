import numpy as np

#material density [kg/m3]
#material strength in cryogenic conditions [MPa]
materials_properties = {
    "Al-Li": {"density": 3000, "strength": 560},
    "Annealed 304L Stainless Steel": {"density": 7800, "strength": 190},
    "Annealed Ti-6Al-4V": {"density": 4400, "strength": 910},
    "Inconel 718": {"density": 8300, "strength": 1190}
}
#Constraints
tank_diameter = 7
LH2_LOX_O_F_ratio = 6/1
LH2_volume_shrinkage_coefficient = 1.2 #%
LOX_volume_shrinkage_coefficient = 1.43 #%
CH4_volume_shrinkage_coefficient = 1.28 #%

LH2_pressure = 1035 #kPa (10 bar)
LH2_density = 0.07085 * 10^3 #kg/m3
LOX_pressure = 2330,47 #kPa (23 bar)
Ch4_pressure = 0

#HydroLOx Jarvis
propellant_mass_jarvis = 14.7 *10^3 #14.7 tons
structural_mass_jarvis = 2.03 *10^3 #2.03 tons
oxydizer_mass_jarvis = 6*propellant_mass_jarvis

#HydroLox Spaceshuttle
propellant_mass_spaceshuttle = 19.1 * 10^3
structural_mass_spaceshuttle = 2.7 * 10^3

#MethaLox Starship
propellant_mass_starship = 18.61 *10^3 #18.61 tons
structural_mass_starship = 1.37 * 10^3 #1.37 tons


def tank_volume():
    ##
    return volume
# calculating tank thickness
def calculate_tank_thickness(tank_diameter, tank_length, propellant_pressure, allowable_stress):
    # Assuming a thin-walled cylinder for simplicity
    # Using the formula: t = (P * D) / (2 * σ)
    # where P is the internal pressure, D is the diameter, and σ is the allowable stress
    thickness = (propellant_pressure * tank_diameter) / (2 * allowable_stress)
    return thickness

def calculate_tank_mass(tank_diameter, tank_length, thickness, material_density):
    # Calculate the volume of the tank material (assuming a cylindrical shape)
    outer_radius = tank_diameter / 2 + thickness
    inner_radius = tank_diameter / 2
    volume = np.pi * (outer_radius**2 - inner_radius**2) * tank_length

    mass = volume * material_density
    return mass

