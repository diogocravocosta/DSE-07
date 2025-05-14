import numpy as np

#material density [kg/m3]
#material strength in cryogenic conditions [MPa]
materials_properties = {
    "Al-Li": {"density": 3000, "strength": 560*10**6, "young modulus": 69*10**9},
    "Annealed 304L Stainless Steel": {"density": 7800, "strength": 190*10**6, "young modulus": 200*10**9},
    "Annealed Ti-6Al-4V": {"density": 4400, "strength": 910*10**6, "young modulus": 110*10**9},
    "Inconel 718": {"density": 8300, "strength": 1190*10**6, "young modulus": 190*10**9}
}
#Constraints
buckling_stress_safety_factor = 1.5
tank_diameter = 7
LH2_LOX_O_F_ratio = 6/1
#LH2_volume_shrinkage_coefficient = 1.2/100 #%
#LOX_volume_shrinkage_coefficient = 1.43/100 #%
#CH4_volume_shrinkage_coefficient = 1.28/100 #%

LH2_pressure = 1035 #kPa (10 bar)
LOX_pressure = 2330.47 #kPa (23 bar)
Ch4_pressure = 400 #kPa

LH2_boiloff_margin = 1.012 * 1.02 #(3.2% extra)
LOX_boiloff_margin = 1.0143 * 1.02 #(3.46% extra)
CH4_boiloff_margin = 1.0128 * 1.02

LH2_density = 70.85 #kg/m3
LOX_density = 1141 #kg/m3
CH4_density = 422.62 #kg/m3


#HydroLOx Jarvis
structural_mass_jarvis = 2.95 * 10**3 #2.03 tons
propellant_mass_jarvis = 21.4 * 10**3 #14.7 tons
propellant_volume_jarvis = propellant_mass_jarvis/LH2_density*LH2_boiloff_margin
oxydizer_volume_jarvis = 6*propellant_mass_jarvis / LOX_density*LOX_boiloff_margin

#HydroLox Spaceshuttle
propellant_mass_spaceshuttle = 21.2 * 10**3
structural_mass_spaceshuttle = 3.03 * 10**3
propellant_volume_spaceshuttle = propellant_mass_spaceshuttle/LH2_density*LH2_boiloff_margin
oxydizer_volume_spaceshuttle = 6*propellant_mass_spaceshuttle/LOX_density*LOX_boiloff_margin

#MethaLox Starship
propellant_mass_starship = 27.16 * 10**3 #18.61 tons
structural_mass_starship = 2 * 10**3 #1.37 tons
propellant_volume_starship = propellant_mass_starship/CH4_density*CH4_boiloff_margin
oxydizer_volume_starship = 6*propellant_mass_starship / LOX_density*LOX_boiloff_margin

def calculate_tank_length(volume, diameter):
    # V_total = V_cylinder + V_hemispheres = π * r^2 * h + (4/3) * π * r^3
    radius = diameter / 2
    hemisphere_volume = (4/3) * np.pi * radius**3
    cylinder_volume = volume - hemisphere_volume
    length = cylinder_volume / (np.pi * radius**2)
    return length

# calculating tank thickness
def calculate_tank_thickness(tank_diameter, tank_length, young_modulus, propellant_pressure, allowable_stress, structural_mass):
    # Assuming a thin-walled cylinder for simplicity
    # Using the formula: t = (P * D) / (2 * σ)
    # where P is the internal pressure, D is the diameter, and σ is the allowable stress
    thickness_pressure = (propellant_pressure * tank_diameter) / (2 * allowable_stress)

    #Axial loading
    # Buckling stress with knockdown factor:
    # sigma_cr = γ * (0.605 * E * t) / R
    gamma_knockdown_factor = 0.9
    critical_buckling_stress = (structural_mass * 8.5)/(np.pi * (tank_diameter/2)**2)
    thickness_load = (critical_buckling_stress*buckling_stress_safety_factor*tank_diameter/2)/(0.605*gamma_knockdown_factor*young_modulus)
    if thickness_pressure > thickness_load:
        return thickness_pressure
    else:
        return thickness_load

def calculate_tank_mass(tank_diameter, tank_length, thickness, material_density):
    volume_shell = 2*np.pi*tank_diameter/2*thickness*(tank_length + 2*tank_diameter/2)
    mass = volume_shell * material_density
    return mass

########################################### Results ################################################################
# Choose a material
material = "Al-Li"

# Access density and strength
density = materials_properties[material]["density"]       # kg/m^3
strength = materials_properties[material]["strength"]     # Pa
young_modulus = materials_properties[material]["young modulus"]   # Pa

tank_length_jarvis = calculate_tank_length(propellant_volume_jarvis, tank_diameter)
tank_length_spaceshuttle = calculate_tank_length(propellant_volume_spaceshuttle, tank_diameter)
print(calculate_tank_length(propellant_volume_spaceshuttle, tank_diameter))
tank_length_starship = calculate_tank_length(propellant_volume_starship, tank_diameter)

#Jarvis tanks
print("Jarvis Tanks:")
print("Jarvis LH2 Tanks Volume:" + str(propellant_volume_jarvis))
print("Jarvis LOX Tanks Volume:" + str(oxydizer_volume_jarvis))
thickness_LH2_jarvis = calculate_tank_thickness(tank_diameter, tank_length_jarvis, young_modulus, LH2_pressure, strength, structural_mass_jarvis)
print("Thickness LH2 Tank Jarvis:" + str(thickness_LH2_jarvis))
thickness_LOX_jarvis = calculate_tank_thickness(tank_diameter, tank_length_jarvis, young_modulus, LOX_pressure, strength, structural_mass_jarvis)
print("Thickness LOX Tank Jarvis:" + str(thickness_LOX_jarvis))

