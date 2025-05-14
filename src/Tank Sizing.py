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
#LH2_volume_shrinkage_coefficient = 1.2/100 #%
#LOX_volume_shrinkage_coefficient = 1.43/100 #%
#CH4_volume_shrinkage_coefficient = 1.28/100 #%

LH2_pressure = 1035 #kPa (10 bar)
LOX_pressure = 2330.47 #kPa (23 bar)
Ch4_pressure = 400 #kPa (4bar)

LH2_boiloff_margin = 1.012 * 1.02 #(3.2% extra)
LOX_boiloff_margin = 1.0143 * 1.02 #(3.46% extra)
CH4_boiloff_margin = 1.0128 * 1.02

LH2_density = 70.85 #kg/m3
LOX_density = 1141 #kg/m3
CH4_density = 422.62 #kg/m3


#HydroLOx Jarvis
structural_mass_jarvis = 2.95 * 10**3 #2.03 tons
propellant_mass_jarvis = 21390.62298
wet_mass_jarvis = propellant_mass_jarvis + structural_mass_jarvis
LH2_mass_jarvis = 1/7.03*propellant_mass_jarvis
print("Mass LH2: " + str(LH2_mass_jarvis) + " kg")
LOX_mass_jarvis = 6.03/7.03*propellant_mass_jarvis
print("Mass LOX: " + str(LOX_mass_jarvis) + " kg")
LH2_volume_jarvis = LH2_mass_jarvis/LH2_density*LH2_boiloff_margin
print("Volume LH2: " + str(LH2_volume_jarvis) + " m3")
LOX_volume_jarvis = LOX_mass_jarvis/ LOX_density*LOX_boiloff_margin
print("Volume LOX: " + str(LOX_volume_jarvis) + " m3")

#HydroLox Spaceshuttle
propellant_mass_spaceshuttle = 21182.68314
structural_mass_spaceshuttle = 3.03 * 10**3
wet_mass_spaceshuttle = propellant_mass_spaceshuttle + structural_mass_spaceshuttle
LH2_mass_spaceshuttle = 1/7.03*propellant_mass_spaceshuttle
LOX_mass_spaceshuttle = 6.03/7.03*propellant_mass_spaceshuttle
LH2_volume_spaceshuttle = LH2_mass_spaceshuttle/LH2_density*LH2_boiloff_margin
LOX_volume_spaceshuttle = LOX_mass_spaceshuttle/ LOX_density*LOX_boiloff_margin

#MethaLox Starship
propellant_mass_starship = 27157.54779
structural_mass_starship = 2 * 10**3 #1.37 tons
wet_mass_starship = propellant_mass_starship + structural_mass_starship
CH4_mass_starship = 1/4.6*propellant_mass_starship
LOX_mass_starship = 3.6/4.6*propellant_mass_starship
CH4_volume_starship = CH4_mass_starship/CH4_density*CH4_boiloff_margin
LOX_volume_starship = LOX_mass_starship/LOX_density*LOX_boiloff_margin

def calculate_tank_length(volume, diameter):
    # V_total = V_cylinder + V_hemispheres = π * r^2 * h + (4/3) * π * r^3
    radius = diameter / 2
    hemisphere_volume = (4/3) * np.pi * radius**3
    cylinder_volume = volume - hemisphere_volume
    length = cylinder_volume / (np.pi * radius**2)
    return length

#print((4/3) * np.pi * (tank_diameter/2)**3)
print("Jarvis length tank: " + str(calculate_tank_length(LH2_volume_spaceshuttle, tank_diameter)))
print(LH2_volume_jarvis - (np.pi * (tank_diameter/2)**2 * 4.34 + 8/3* np.pi * (tank_diameter/2)**3))
# calculating tank thickness
def calculate_tank_thickness(tank_diameter, tank_length, young_modulus, propellant_pressure, allowable_stress, wet_mass):
    # Assuming a thin-walled cylinder for simplicity
    # Using the formula: t = (P * D) / (2 * σ)
    # where P is the internal pressure, D is the diameter, and σ is the allowable stress
    thickness_pressure = (propellant_pressure * tank_diameter) / (2 * allowable_stress)

    #Axial loading
    # Buckling stress with knockdown factor:
    # sigma_cr = γ * (0.605 * E * t) / R
    gamma_knockdown_factor = 0.9
    thickness_load = np.sqrt((wet_mass* 8.5 * 9.81 * buckling_stress_safety_factor)/(2*np.pi)*1/(0.605*gamma_knockdown_factor*young_modulus))
    if thickness_pressure > thickness_load:
        return thickness_pressure
    else:
        return thickness_load

def calculate_tank_mass(tank_diameter, tank_length, thickness, material_density):
    volume_shell = 2*np.pi*tank_diameter/2*thickness*(tank_length + 2*tank_diameter/2)
    mass = volume_shell * material_density
    return mass

def calculate_margin_of_safety(thickness, tank_diameter, wet_mass, E, gamma=0.9, g_load=8.5):
    g = 9.81  # gravity
    # Applied stress
    A_wall = 2 * np.pi * (tank_diameter/2) * thickness
    sigma_applied = (wet_mass * g * g_load) / A_wall

    # Critical buckling stress
    sigma_cr = gamma * 0.605 * E * thickness / (tank_diameter/2)

    # Margin of safety
    MoS = sigma_cr / sigma_applied - 1
    return MoS

########################################### Results ################################################################
# Choose a material
material = "Al-Li"

# Access density and strength
density = materials_properties[material]["density"]       # kg/m^3
strength = materials_properties[material]["strength"]     # Pa
young_modulus = materials_properties[material]["young modulus"]   # Pa

tank_length_jarvis_LH2 = calculate_tank_length(LH2_volume_jarvis, tank_diameter)
tank_length_jarvis_LOX = calculate_tank_length(LOX_volume_jarvis, tank_diameter)
tank_length_spaceshuttle_LH2 = calculate_tank_length(LH2_volume_spaceshuttle, tank_diameter)
tank_length_spaceshuttle_LOX = calculate_tank_length(LOX_volume_spaceshuttle, tank_diameter)
tank_length_starship_CH4 = calculate_tank_length(CH4_volume_starship, tank_diameter)
tank_length_starship_LOX = calculate_tank_length(LOX_volume_starship, tank_diameter)

#Jarvis tanks
print("Jarvis Tanks:")
print("Jarvis LH2 Tanks Volume: " + str(LH2_volume_jarvis))
print("Jarvis LOX Tanks Volume: " + str(LOX_volume_jarvis))
thickness_LH2_jarvis = calculate_tank_thickness(tank_diameter, tank_length_jarvis_LH2, young_modulus, LH2_pressure, strength, wet_mass_jarvis)
print("Thickness LH2 Tank Jarvis: " + str(thickness_LH2_jarvis))
thickness_LOX_jarvis = calculate_tank_thickness(tank_diameter, tank_length_jarvis_LOX, young_modulus, LOX_pressure, strength, wet_mass_jarvis)
print("Thickness LOX Tank Jarvis: " + str(thickness_LOX_jarvis))
mass_LH2_tank_jarvis = calculate_tank_mass(tank_diameter, tank_length_jarvis_LH2, thickness_LH2_jarvis, density)
print("Mass LH2 Tank Jarvis: " + str(mass_LH2_tank_jarvis))
mass_LOX_tank_jarvis = calculate_tank_mass(tank_diameter, tank_length_jarvis_LOX, thickness_LOX_jarvis, density)
print("Mass LOX Tank Jarvis: " + str(mass_LOX_tank_jarvis))

print("=============================================================================")
#Space Shuttle tanks
print("Space Shuttle Tanks:")
print("Space Shuttle LH2 Tanks Volume: " + str(LH2_volume_spaceshuttle))
print("Space Shuttle LOX Tanks Volume: " + str(LOX_volume_spaceshuttle))
thickness_LH2_spaceshuttle = calculate_tank_thickness(tank_diameter, tank_length_spaceshuttle_LH2, young_modulus, LH2_pressure, strength, wet_mass_spaceshuttle)
print("Thickness LH2 Tank Space Shuttle: " + str(thickness_LH2_spaceshuttle))
thickness_LOX_spaceshuttle = calculate_tank_thickness(tank_diameter, tank_length_spaceshuttle_LOX, young_modulus, LOX_pressure, strength, wet_mass_spaceshuttle)
print("Thickness LOX Tank Space Shuttle: " + str(thickness_LOX_spaceshuttle))
mass_LH2_tank_spaceshuttle = calculate_tank_mass(tank_diameter, tank_length_spaceshuttle_LH2, thickness_LH2_spaceshuttle, density)
print("Mass LH2 Tank Space Shuttle: " + str(mass_LH2_tank_spaceshuttle))
mass_LOX_tank_spaceshuttle = calculate_tank_mass(tank_diameter, tank_length_spaceshuttle_LOX, thickness_LOX_spaceshuttle, density)
print("Mass LOX Tank Space Shuttle: " + str(mass_LOX_tank_spaceshuttle))

print("=============================================================================")
#Starship tanks
print("Starship Tanks:")
print("Starship CH4 Tanks Volume: " + str(CH4_volume_starship))
print("Starship LOX Tanks Volume: " + str(LOX_volume_starship))
thickness_CH4_starship = calculate_tank_thickness(tank_diameter, tank_length_starship_CH4, young_modulus, Ch4_pressure, strength, wet_mass_starship)
print("Thickness LH2 Tank Starship: " + str(thickness_CH4_starship))
thickness_LOX_starship = calculate_tank_thickness(tank_diameter, tank_length_starship_LOX, young_modulus, LOX_pressure, strength, wet_mass_starship)
print("Thickness LOX Tank Starship: " + str(thickness_LOX_starship))
mass_CH4_tank_starship = calculate_tank_mass(tank_diameter, tank_length_starship_CH4, thickness_CH4_starship, density)
print("Mass CH4 Tank Starship: " + str(mass_CH4_tank_starship))
mass_LOX_tank_starship = calculate_tank_mass(tank_diameter, tank_length_starship_LOX, thickness_LOX_starship, density)
print("Mass LOX Tank Starship: " + str(mass_LOX_tank_starship))

print("=============================================================================")
SM_jarvis_LH2 = calculate_margin_of_safety(thickness_LH2_jarvis, tank_diameter, wet_mass_jarvis, young_modulus, gamma=0.9, g_load=8.5)
print("SM LH2 Jarvis: " + str(calculate_margin_of_safety(thickness_LH2_jarvis, tank_diameter, wet_mass_jarvis, young_modulus, gamma=0.9, g_load=8.5)))
SM_jarvis_LOX = calculate_margin_of_safety(thickness_LH2_jarvis, tank_diameter, wet_mass_jarvis, young_modulus, gamma=0.9, g_load=8.5)
print("SM LOX Jarvis: " +str(calculate_margin_of_safety(thickness_LH2_jarvis, tank_diameter, wet_mass_jarvis, young_modulus, gamma=0.9, g_load=8.5)))
SM_spaceshuttle_LH2 = calculate_margin_of_safety(thickness_LH2_spaceshuttle, tank_diameter, wet_mass_spaceshuttle, young_modulus, gamma=0.9, g_load=8.5)
print("SM LH2 Space Shuttle: " + str(calculate_margin_of_safety(thickness_LH2_spaceshuttle, tank_diameter, wet_mass_spaceshuttle, young_modulus, gamma=0.9, g_load=8.5)))
SM_spaceshuttle_LOX = calculate_margin_of_safety(thickness_LOX_spaceshuttle, tank_diameter, wet_mass_spaceshuttle, young_modulus, gamma=0.9, g_load=8.5)
print("SM LOX Space Shuttle: " +str(calculate_margin_of_safety(thickness_LOX_spaceshuttle, tank_diameter, wet_mass_spaceshuttle, young_modulus, gamma=0.9, g_load=8.5)))
SM_starship_CH4 = calculate_margin_of_safety(thickness_CH4_starship, tank_diameter, wet_mass_starship, young_modulus, gamma=0.9, g_load=8.5)
print("SM CH4 Starship: " + str(calculate_margin_of_safety(thickness_CH4_starship, tank_diameter, wet_mass_starship, young_modulus, gamma=0.9, g_load=8.5)))
SM_starship_LOX = calculate_margin_of_safety(thickness_LOX_starship, tank_diameter, wet_mass_starship, young_modulus, gamma=0.9, g_load=8.5)
print("SM LOX Space Shuttle: " +str(calculate_margin_of_safety(thickness_LOX_starship, tank_diameter, wet_mass_starship, young_modulus, gamma=0.9, g_load=8.5)))
print(0.603*(0.9*young_modulus*thickness_LH2_jarvis)/(tank_diameter))



def size_tank_thickness(model, wet_mass, propellant_pressure, fuel_mass, oxidizer_mass, tank_length_LH2, tank_length_CH4, tank_length_LOX, tank_diameter, E, strength, gamma=0.9):
    t = 0.001  # Start with 1 mm
    while True:
        sigma_axial = wet_mass* 8.5 * 9.81 * buckling_stress_safety_factor / (2 * np.pi * tank_diameter/2 * t)
        #cg location
        if model == "LH2":
            cg_location = (fuel_mass * (tank_length_LH2 + tank_length_LOX) + oxidizer_mass * tank_length_LOX)/(fuel_mass + oxidizer_mass)
            sigma_bend = (3*9.81* (fuel_mass + oxidizer_mass) * cg_location) / (2 * np.pi * (tank_diameter/2)**2 * t)
        elif model == "CH4":
            cg_location = (fuel_mass * (tank_length_CH4) + 12500 * (tank_length_CH4 + tank_length_LOX + tank_length_LH2) + oxidizer_mass * (tank_length_CH4 + tank_length_LOX)) / (fuel_mass + oxidizer_mass + 12500)
            sigma_bend = (3 * 9.81 * (fuel_mass + oxidizer_mass) * cg_location) / (2 * np.pi * (tank_diameter / 2) ** 2 * t)

        sigma_cr = gamma * 0.605 * E * t / (tank_diameter/2)
        interaction = (sigma_axial + sigma_bend) / sigma_cr
        if interaction <= 1.0:
            break
        t += 0.0001  # Increment 0.1 mm
    #check hoop stress
    thickness_pressure = (propellant_pressure * tank_diameter) / (2 * strength)
    return t

print(size_tank_thickness("LH2", wet_mass_jarvis, propellant_pressure, fuel_mass, oxidizer_mass, tank_length_LH2, tank_length_CH4, tank_length_LOX, tank_diameter, E, strength, gamma=0.9))