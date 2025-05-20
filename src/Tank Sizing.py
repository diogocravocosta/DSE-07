import numpy as np

#material density [kg/m3]
#material strength in cryogenic conditions [MPa]
materials_properties = {
    "Al-Li": {"density": 3000, "strength": 730*10**6, "young modulus": 75*10**9},
    "Annealed 304L Stainless Steel": {"density": 7800, "strength": 1060*10**6, "young modulus": 200*10**9},
    "Annealed Ti-6Al-4V": {"density": 4400, "strength": 1100*10**6, "young modulus": 115*10**9},
    "Inconel 718": {"density": 8300, "strength": 1190*10**6, "young modulus": 190*10**9},
    "CFRP": {"density": 1500, "strength": 1500*10**6, "young modulus": 140*10**9}
}
#Constraints
safety_factor = 1.5
tank_diameter = 7
#LH2_volume_shrinkage_coefficient = 1.2/100 #%
#LOX_volume_shrinkage_coefficient = 1.43/100 #%
#CH4_volume_shrinkage_coefficient = 1.28/100 #%

LH2_pressure = 500 #kPa (5 bar)
LOX_pressure = 270 #kPa (2.7 bar)
Ch4_pressure = 400 #kPa (4bar)

LH2_boiloff_margin = 1.1 #1.012 * 1.02 #(3.2% extra)
LOX_boiloff_margin = 1.0143 * 1.02 #(3.46% extra)
CH4_boiloff_margin = 1.0128 * 1.02

LH2_density = 70.85 #kg/m3
LOX_density = 1141 #kg/m3
CH4_density = 422.62 #kg/m3

print("HydroLox Jarvis")
#HydroLOx Jarvis
structural_mass_jarvis_LH2 = 33648.68751
propellant_mass_jarvis_LH2 = 244372.1506
wet_mass_jarvis_LH2 = propellant_mass_jarvis_LH2 + structural_mass_jarvis_LH2
LH2_mass_jarvis_LH2 = 1/7.03*propellant_mass_jarvis_LH2 +18000
print("Mass LH2: " + str(LH2_mass_jarvis_LH2) + " kg")
LOX_mass_jarvis_LH2 = 6.03/7.03*propellant_mass_jarvis_LH2
print("Mass LOX: " + str(LOX_mass_jarvis_LH2) + " kg")
LH2_volume_jarvis_LH2 = LH2_mass_jarvis_LH2/LH2_density*LH2_boiloff_margin
print("Volume LH2: " + str(LH2_volume_jarvis_LH2) + " m3")
LOX_volume_jarvis_LH2 = LOX_mass_jarvis_LH2/ LOX_density*LOX_boiloff_margin
print("Volume LOX: " + str(LOX_volume_jarvis_LH2) + " m3")
print("-----------------------------------------------------------------------------------------------------------------")

print("MethaLOX Jarvis")
#MethaLOx Jarvis
structural_mass_jarvis_CH4 = 57257.86121
propellant_mass_jarvis_CH4 = 528534.1034
payload_mass_jarvis_CH4 = 18000
wet_mass_jarvis_CH4 = propellant_mass_jarvis_CH4 + structural_mass_jarvis_CH4 + payload_mass_jarvis_CH4
CH4_mass_jarvis = 1/4.6*propellant_mass_jarvis_CH4
print("Mass CH4: " + str(CH4_mass_jarvis) + " kg")
LOX_mass_jarvis_CH4 = 3.6/4.6*propellant_mass_jarvis_CH4
print("Mass LOX: " + str(LOX_mass_jarvis_CH4) + " kg")
CH4_volume_jarvis = CH4_mass_jarvis/CH4_density*CH4_boiloff_margin
print("Volume CH4: " + str(CH4_volume_jarvis) + " m3")
LOX_volume_jarvis_CH4 = LOX_mass_jarvis_CH4/LOX_density*LOX_boiloff_margin
print("Volume LOX: " + str(LOX_volume_jarvis_CH4) + " m3")
LH2_volume_jarvis_CH4 = payload_mass_jarvis_CH4/LH2_density*LH2_boiloff_margin
print("Volume LH2: " + str(LH2_volume_jarvis_CH4) + " m3")
print("-----------------------------------------------------------------------------------------------------------------")

print("MethaLOX Starship")
#MethaLox Starship
propellant_mass_starship_CH4 = 337413.4919
payload_mass_starship_CH4 = 18000
structural_mass_starship_CH4 =  36553.12829
wet_mass_starship_CH4 = propellant_mass_starship_CH4 + structural_mass_starship_CH4 + payload_mass_starship_CH4
CH4_mass_starship = 1/4.6*propellant_mass_starship_CH4
print("Mass CH4: " + str(CH4_mass_starship) + " kg")
LOX_mass_starship_CH4 = 3.6/4.6*propellant_mass_starship_CH4
print("Mass LOX: " + str(LOX_mass_starship_CH4) + " kg")
CH4_volume_starship = CH4_mass_starship/CH4_density*CH4_boiloff_margin
print("Volume CH4: " + str(CH4_volume_starship) + " m3")
LOX_volume_starship_CH4 = LOX_mass_starship_CH4/LOX_density*LOX_boiloff_margin
print("Volume LOX: " + str(LOX_volume_starship_CH4) + " m3")
LH2_volume_starship_CH4 = payload_mass_starship_CH4/LH2_density*LH2_boiloff_margin
print("Volume LH2: " + str(LH2_volume_starship_CH4) + " m3")
print("-----------------------------------------------------------------------------------------------------------------")

print("HydroLOX Starship")
#HydroLOX Starship
propellant_mass_starship_LH2 = 179895.8284
structural_mass_starship_LH2 = 24770.65615
wet_mass_starship_LH2 = propellant_mass_starship_LH2 + structural_mass_starship_LH2
LH2_mass_starship_LH2 = 1/7.03*propellant_mass_starship_LH2 + 18000
print("Mass LH2: " + str(LH2_mass_starship_LH2) + " kg")
LOX_mass_starship_LH2 = 6.03/7.03*propellant_mass_starship_LH2
print("Mass LOX: " + str(LOX_mass_starship_LH2) + " kg")
LH2_volume_starship_LH2 = LH2_mass_starship_LH2/LH2_density*LH2_boiloff_margin
print("Volume Liquid LH2: " + str(LH2_mass_starship_LH2/LH2_density))
print("Volume LH2: " + str(LH2_volume_starship_LH2) + " m3")
LOX_volume_starship_LH2 = LOX_mass_starship_LH2/LOX_density*LOX_boiloff_margin
print("Volume LOX: " + str(LOX_volume_starship_LH2) + " m3")
print("-----------------------------------------------------------------------------------------------------------------")
print("-----------------------------------------------------------------------------------------------------------------")
print("-----------------------------------------------------------------------------------------------------------------")
def calculate_tank_length(model, tank_model, volume, tank_diameter):
    # V_total = V_cylinder + V_ellipsoidal_caps = π * r^2 * h + 2 * (π/24) * D^3
    if model == "LH2_LOX":
        if tank_model == "LOX":
            radius = tank_diameter / 2
            cap_volume = (np.pi / 24) * tank_diameter ** 3
            ellipsoidal_caps_volume = 2 * cap_volume
            cylinder_volume = volume - ellipsoidal_caps_volume
            length = cylinder_volume / (np.pi * radius**2)
        elif tank_model == "LH2":
            radius = tank_diameter / 2
            cap_volume = (np.pi / 24) * tank_diameter ** 3
            cylinder_volume = volume - cap_volume + cap_volume
            length = cylinder_volume / (np.pi * radius ** 2)
    elif model == "CH4_LOX":
        if tank_model == "CH4":
            radius = tank_diameter / 2
            cap_volume = (np.pi / 24) * tank_diameter ** 3
            ellipsoidal_caps_volume = 2 * cap_volume
            cylinder_volume = volume - ellipsoidal_caps_volume
            length = cylinder_volume / (np.pi * radius**2)
        elif tank_model == "LOX":
            radius = tank_diameter / 2
            cap_volume = (np.pi / 24) * tank_diameter ** 3
            cylinder_volume = volume - cap_volume + cap_volume
            length = cylinder_volume / (np.pi * radius ** 2)
        elif tank_model == "LH2":
            radius = tank_diameter / 2
            cap_volume = (np.pi / 24) * tank_diameter ** 3
            cylinder_volume = volume - cap_volume + cap_volume
            length = cylinder_volume / (np.pi * radius ** 2)
    return length

print("Jarvis length tank: " + str(calculate_tank_length("LH2_LOX", "LH2", LH2_volume_jarvis_LH2, tank_diameter)))
print("Jarvis volume tank: " + str(np.pi * (tank_diameter/2)**2 * calculate_tank_length("LH2_LOX", "LH2", LH2_volume_jarvis_LH2, tank_diameter)))
print("Jarvis volume LH2: " + str(LH2_volume_jarvis_LH2))
def calculate_tank_thickness(wet_mass, propellant_pressure, fuel_mass, tank_length, tank_diameter, E, strength, gamma=0.65):
    t = 0.001  # Start with 1 mm
    while True:
        #axial stress due to the launch acceleration
        sigma_axial = wet_mass* 8.5 * 9.81 / (2 * np.pi * tank_diameter/2 * t)
        #bending stress due to lateral loads
        sigma_bend = (3*9.81* fuel_mass * tank_length/2) / (2 * np.pi * (tank_diameter/2)**2 * t)
        sigma_cr = gamma * 0.605 * E * t / (tank_diameter/2)

        interaction = safety_factor * (sigma_axial + sigma_bend) / sigma_cr

        if interaction <= 1.0:
            break
        t += 0.0001  # Increment 0.1 mm
    print("the axial stress is: " + str(sigma_axial))
    print("the bending stress is: " + str(sigma_bend))
    #check hoop stress
    thickness_pressure = (propellant_pressure * tank_diameter) / (2 * strength)
    if thickness_pressure>t:
        t = thickness_pressure
        print("Hoop stress leading")
    return t

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
material = "Annealed 304L Stainless Steel"

# Access density and strength
density = materials_properties[material]["density"]       # kg/m^3
strength = materials_properties[material]["strength"]     # Pa
young_modulus = materials_properties[material]["young modulus"]   # Pa

tank_length_jarvis_LH2 = calculate_tank_length("LH2_LOX", "LH2", LH2_volume_jarvis_LH2, tank_diameter)
tank_length_jarvis_LH2_LOX = calculate_tank_length("LH2_LOX", "LOX", LOX_volume_jarvis_LH2, tank_diameter)

tank_length_jarvis_CH4 = calculate_tank_length("CH4_LOX", "CH4", CH4_volume_jarvis, tank_diameter)
tank_length_jarvis_CH4_LOX = calculate_tank_length("CH4_LOX", "LOX", LOX_volume_jarvis_CH4, tank_diameter)
tank_length_jarvis_CH4_LH2 = calculate_tank_length("CH4_LOX", "LH2", LH2_volume_jarvis_CH4, tank_diameter)

tank_length_starship_CH4 = calculate_tank_length("CH4_LOX", "CH4", CH4_volume_starship, tank_diameter)
tank_length_starship_CH4_LOX = calculate_tank_length("CH4_LOX", "LOX", LOX_volume_starship_CH4, tank_diameter)
tank_length_starship_CH4_LH2 = calculate_tank_length("CH4_LOX", "LH2", LH2_volume_starship_CH4, tank_diameter)

tank_length_starship_LH2 = calculate_tank_length("LH2_LOX", "LH2", LH2_volume_starship_LH2, tank_diameter)
print("############### Tank LH2 length: " + str(tank_length_starship_LH2))
print("Volume HydroLox Tank: " + str(LOX_volume_starship_LH2))
tank_length_starship_LH2_LOX = calculate_tank_length("LH2_LOX", "LOX", LOX_volume_starship_LH2, tank_diameter)
#print("TANK VOLUME: " + str(np.pi * (tank_diameter/2)^2 * tank_length_jarvis_LH2_LOX + 2 * (np.pi/24) * (tank_diameter)^3))
print("############### Tank LOX length: " + str(tank_length_starship_LH2_LOX))
#Jarvis HydroLOX tanks
print("Jarvis HydroLox Tanks:")
print("Jarvis LH2 Tanks Volume: " + str(LH2_volume_jarvis_LH2))
print("Jarvis LOX Tanks Volume: " + str(LOX_volume_jarvis_LH2))
thickness_LH2_jarvis_LH2 = calculate_tank_thickness(wet_mass_jarvis_LH2, LH2_pressure, LH2_mass_jarvis_LH2, tank_length_jarvis_LH2, tank_diameter, young_modulus, strength, gamma=0.65)
print("Thickness LH2 Tank Jarvis: " + str(thickness_LH2_jarvis_LH2))
thickness_LOX_jarvis_LH2 = calculate_tank_thickness(wet_mass_jarvis_LH2, LOX_pressure, LOX_mass_jarvis_LH2, tank_length_jarvis_LH2_LOX, tank_diameter, young_modulus, strength, gamma=0.65)
print("Thickness LOX Tank Jarvis: " + str(thickness_LOX_jarvis_LH2))
mass_LH2_tank_jarvis_LH2 = calculate_tank_mass(tank_diameter, tank_length_jarvis_LH2, thickness_LH2_jarvis_LH2, density)
print("Mass LH2 Tank Jarvis: " + str(mass_LH2_tank_jarvis_LH2))
mass_LOX_tank_jarvis_LH2 = calculate_tank_mass(tank_diameter, tank_length_jarvis_LH2_LOX, thickness_LOX_jarvis_LH2, density)
print("Mass LOX Tank Jarvis: " + str(mass_LOX_tank_jarvis_LH2))

print("=============================================================================")

#Jarvis MethaLOX tanks
print("Jarvis MethaLOX Tanks:")
print("Jarvis CH4 Tanks Volume: " + str(CH4_volume_jarvis))
print("Jarvis LH2 Tanks Volume: " + str(LH2_volume_jarvis_CH4))
print("Jarvis LOX Tanks Volume: " + str(LOX_volume_jarvis_CH4))
thickness_LH2_jarvis_CH4 = calculate_tank_thickness(wet_mass_jarvis_CH4, LH2_pressure, payload_mass_jarvis_CH4, tank_length_jarvis_CH4_LH2, tank_diameter, young_modulus, strength, gamma=0.65)
print("Thickness LH2 Tank Jarvis: " + str(thickness_LH2_jarvis_CH4))
thickness_CH4_jarvis_CH4 = calculate_tank_thickness(wet_mass_jarvis_CH4, LH2_pressure, CH4_mass_jarvis, tank_length_jarvis_CH4, tank_diameter, young_modulus, strength, gamma=0.65)
print("Thickness CH4 Tank Jarvis: " + str(thickness_CH4_jarvis_CH4))
thickness_LOX_jarvis_CH4 = calculate_tank_thickness(wet_mass_jarvis_CH4, LOX_pressure, LOX_mass_jarvis_CH4, tank_length_jarvis_CH4_LOX, tank_diameter, young_modulus, strength, gamma=0.65)
print("Thickness LOX Tank Jarvis: " + str(thickness_LOX_jarvis_CH4))
mass_LH2_tank_jarvis_CH4 = calculate_tank_mass(tank_diameter, tank_length_jarvis_CH4_LH2, thickness_LH2_jarvis_CH4, density)
print("Mass LH2 Tank Jarvis: " + str(mass_LH2_tank_jarvis_CH4))
mass_CH4_tank_jarvis_CH4 = calculate_tank_mass(tank_diameter, tank_length_jarvis_CH4, thickness_CH4_jarvis_CH4, density)
print("Mass CH4 Tank Jarvis: " + str(mass_CH4_tank_jarvis_CH4))
mass_LOX_tank_jarvis_CH4 = calculate_tank_mass(tank_diameter, tank_length_jarvis_CH4_LOX, thickness_LOX_jarvis_CH4, density)
print("Mass LOX Tank Jarvis: " + str(mass_LOX_tank_jarvis_CH4))

print("=============================================================================")

#Starship MethaLOX tanks
print("Starship MethaLOX Tanks:")
print("Starship CH4 Tanks Volume: " + str(CH4_volume_starship))
print("Starship LOX Tanks Volume: " + str(LOX_volume_starship_CH4))
print("Starship LH2 Tanks Volume: " + str(LH2_volume_starship_CH4))
thickness_CH4_starship_CH4 = calculate_tank_thickness(wet_mass_starship_CH4, Ch4_pressure, CH4_mass_starship, tank_length_starship_CH4, tank_diameter, young_modulus, strength, gamma=0.65)
print("Thickness CH4 Tank Starship: " + str(thickness_CH4_starship_CH4))
thickness_LOX_starship_CH4 = calculate_tank_thickness(wet_mass_starship_CH4, LOX_pressure, LOX_mass_starship_CH4, tank_length_starship_CH4_LOX, tank_diameter, young_modulus, strength, gamma=0.65)
print("Thickness LOX Tank Starship: " + str(thickness_LOX_starship_CH4))
thickness_LH2_starship_CH4 = calculate_tank_thickness(wet_mass_starship_CH4, LH2_pressure, payload_mass_starship_CH4, tank_length_starship_CH4_LH2, tank_diameter, young_modulus, strength, gamma=0.65)
print("Thickness LH2 Tank Starship: " + str(thickness_LH2_starship_CH4))
mass_CH4_tank_starship_CH4 = calculate_tank_mass(tank_diameter, tank_length_starship_CH4, thickness_CH4_starship_CH4, density)
print("Mass CH4 Tank Starship: " + str(mass_CH4_tank_starship_CH4))
mass_LOX_tank_starship_CH4 = calculate_tank_mass(tank_diameter, tank_length_starship_CH4_LOX, thickness_LOX_starship_CH4, density)
print("Mass LOX Tank Starship: " + str(mass_LOX_tank_starship_CH4))
mass_LH2_tank_starship_CH4 = calculate_tank_mass(tank_diameter, tank_length_starship_CH4_LH2, thickness_LH2_starship_CH4, density)
print("Mass LH2 Tank Starship: " + str(mass_LH2_tank_starship_CH4))

print("=============================================================================")

#Starship HydroLOX tanks
print("Starship HydroLOX Tanks:")
print("Starship LH2 Tanks Mass: " + str(LH2_mass_starship_LH2))
print("Starship LH2 Tanks Length: " + str(tank_length_starship_LH2))
print("Starship LOX Tanks Volume: " + str(LOX_volume_starship_LH2))
print("Starship LH2 Tanks Volume: " + str(LH2_volume_starship_LH2))
thickness_LOX_starship_LH2 = calculate_tank_thickness(wet_mass_starship_LH2, LOX_pressure, LOX_mass_starship_LH2, tank_length_starship_LH2_LOX, tank_diameter, young_modulus, strength, gamma=0.65)
print("Thickness LOX Tank Starship: " + str(thickness_LOX_starship_LH2))
thickness_LH2_starship_LH2 = calculate_tank_thickness(wet_mass_starship_LH2, LH2_pressure, LH2_mass_starship_LH2, tank_length_starship_LH2, tank_diameter, young_modulus, strength, gamma=0.65)
print("Thickness LH2 Tank Starship: " + str(thickness_LH2_starship_LH2))
mass_LOX_tank_starship_LH2 = calculate_tank_mass(tank_diameter, tank_length_starship_LH2_LOX, thickness_LOX_starship_LH2, density)
print("Mass LOX Tank Starship: " + str(mass_LOX_tank_starship_LH2))
mass_LH2_tank_starship_LH2= calculate_tank_mass(tank_diameter, tank_length_starship_LH2, thickness_LH2_starship_LH2, density)
print("Mass LH2 Tank Starship: " + str(mass_LH2_tank_starship_LH2))

print("=============================================================================")
# Jarvis HydroLOX
SM_jarvis_LH2 = calculate_margin_of_safety(thickness_LH2_jarvis_LH2, tank_diameter, wet_mass_jarvis_LH2, young_modulus, gamma=0.9)
print("SM LH2 Jarvis: " + str(SM_jarvis_LH2))

SM_jarvis_LOX = calculate_margin_of_safety(thickness_LOX_jarvis_LH2, tank_diameter, wet_mass_jarvis_LH2, young_modulus, gamma=0.9)
print("SM LOX Jarvis: " + str(SM_jarvis_LOX))

# Starship CH4
SM_starship_CH4 = calculate_margin_of_safety(thickness_CH4_starship_CH4, tank_diameter, wet_mass_starship_CH4, young_modulus, gamma=0.9)
print("SM CH4 Starship: " + str(SM_starship_CH4))

SM_starship_LOX = calculate_margin_of_safety(thickness_LOX_starship_CH4, tank_diameter, wet_mass_starship_CH4, young_modulus, gamma=0.9)
print("SM LOX Starship: " + str(SM_starship_LOX))

# Starship HydroLOX
SM_starship_LH2 = calculate_margin_of_safety(thickness_LH2_starship_LH2, tank_diameter, wet_mass_starship_LH2, young_modulus, gamma=0.9)
print("SM LH2 Starship: " + str(SM_starship_LH2))

SM_starship_LOX_LH2 = calculate_margin_of_safety(thickness_LOX_starship_LH2, tank_diameter, wet_mass_starship_LH2, young_modulus, gamma=0.9)
print("SM LOX Starship (HydroLOX): " + str(SM_starship_LOX_LH2))


def calculate_outer_area(tank_length, tank_diameter):
    radius = tank_diameter / 2
    area_cylinder = 2 * np.pi * radius * tank_length
    area_caps = 2 * np.pi * radius**2  # two ellipsoidal caps
    total_area = area_cylinder + area_caps
    return total_area

def calculate_double_ellipse_area(tank_length, tank_diameter):
    radius = tank_diameter / 2
    area_caps = 4 * np.pi * radius ** 2  # two ellipsoidal caps
    total_area = area_caps
    return total_area
print("Outer Area Starship LH2" + str(calculate_outer_area(tank_length_starship_LH2, tank_diameter) + calculate_outer_area(tank_length_starship_LH2_LOX, tank_diameter)- 0.5*calculate_double_ellipse_area(tank_length_starship_LH2_LOX, tank_diameter)))
print("Outer Area Starship CH4" + str(calculate_outer_area(tank_length_starship_CH4_LH2, tank_diameter) + calculate_outer_area(tank_length_starship_CH4, tank_diameter) + calculate_outer_area(tank_length_starship_CH4_LOX, tank_diameter) - 2*calculate_double_ellipse_area(tank_length_starship_CH4_LOX, tank_diameter)))

print("Outer Area Jarvis LH2" + str(calculate_outer_area(tank_length_jarvis_LH2, tank_diameter) + calculate_outer_area(tank_length_jarvis_LH2_LOX, tank_diameter)- 0.5*calculate_double_ellipse_area(tank_length_jarvis_LH2_LOX, tank_diameter)))
print("Outer Area Jarvis CH4" + str(calculate_outer_area(tank_length_jarvis_CH4_LH2, tank_diameter) + calculate_outer_area(tank_length_jarvis_CH4, tank_diameter) + calculate_outer_area(tank_length_jarvis_CH4_LOX, tank_diameter) - 2*calculate_double_ellipse_area(tank_length_jarvis_CH4_LOX, tank_diameter)))
######################################################### OLD CODE ########################################################################

#def calculate_tank_thickness(tank_diameter, tank_length, young_modulus, propellant_pressure, allowable_stress, wet_mass):
    # Assuming a thin-walled cylinder for simplicity
    # Using the formula: t = (P * D) / (2 * σ)
    # where P is the internal pressure, D is the diameter, and σ is the allowable stress
    #thickness_pressure = (propellant_pressure * tank_diameter) / (2 * allowable_stress)

    #Axial loading
    # Buckling stress with knockdown factor:
    # sigma_cr = γ * (0.605 * E * t) / R
    #gamma_knockdown_factor = 0.9
    #thickness_load = np.sqrt((wet_mass* 8.5 * 9.81 * buckling_stress_safety_factor)/(2*np.pi)*1/(0.605*gamma_knockdown_factor*young_modulus))
    #if thickness_pressure > thickness_load:
        #return thickness_pressure
    #else:
        #return thickness_load

