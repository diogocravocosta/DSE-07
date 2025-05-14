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
buckling_stress_safety_factor = 1.5
tank_diameter = 7
LH2_LOX_O_F_ratio = 6/1
#LH2_volume_shrinkage_coefficient = 1.2/100 #%
#LOX_volume_shrinkage_coefficient = 1.43/100 #%
#CH4_volume_shrinkage_coefficient = 1.28/100 #%

LH2_pressure = 1035 #kPa (10 bar)
LOX_pressure = 2330,47 #kPa (23 bar)
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
print(propellant_volume_spaceshuttle)
print(calculate_tank_length(propellant_volume_spaceshuttle, tank_diameter))

# calculating tank thickness
def calculate_tank_thickness(tank_diameter, tank_length, young_modulus, propellant_pressure, allowable_stress):
    # Assuming a thin-walled cylinder for simplicity
    # Using the formula: t = (P * D) / (2 * σ)
    # where P is the internal pressure, D is the diameter, and σ is the allowable stress
    thickness_pressure = (propellant_pressure * tank_diameter) / (2 * allowable_stress)
    gamma_knockdown_factor = 0.9
    thickness_load = (allowable_stress*buckling_stress_safety_factor*tank_diameter/2)/(0.605*0.9*young_modulus)
    if thickness_pressure > thickness_load:
        return thickness_pressure
    else:
        return thickness_load

def calculate_tank_mass(tank_diameter, tank_length, thickness, material_density):
    volume_shell = 2*np.pi *tank_diameter/2*thickness*(tank_length + 2* tank_diameter/2)
    mass = volume_shell * material_density
    return mass

