import numpy as np

materials_properties = {
    "Annealed 304L Stainless Steel": {"density": 7800, "strength": 1060*10**6, "young modulus": 200*10**9},
}
#Geometry
radius_ratio = 0.2
phi = 10

#Constraints
safety_factor = 2.0
safety_factor_pressure = 1
tank_diameter = 7
payload_mass = 15000
thrust_engines = 2129573.909

LH2_pressure = 1000* safety_factor_pressure *1000#kPa (10bar)
LOX_pressure = 270*2 * 1000#kPa (2.7 bar)

LH2_boiloff_margin = 1.1 #1.012 * 1.02 #(3.2% extra)
LOX_boiloff_margin = 1.0143 * 1.02 #(3.46% extra)

LH2_density = 77 #kg/m3
LOX_density = 1340 #kg/m3

structural_mass = 20642.21346
propellant_mass = 149913.1903
wet_mass = propellant_mass + structural_mass
LH2_mass = 1/7.0*propellant_mass + payload_mass
print("Mass LH2: " + str(LH2_mass) + " kg")
LOX_mass = 6.0/7.0*propellant_mass
print("Mass LOX: " + str(LOX_mass) + " kg")
LH2_volume = LH2_mass/LH2_density*LH2_boiloff_margin
print("Volume LH2: " + str(LH2_volume) + " m3")
LOX_volume = LOX_mass/ LOX_density*LOX_boiloff_margin
print("Volume LOX: " + str(LOX_volume) + " m3")
print("-----------------------------------------------------------------------------------------------------------------")

def calculate_tank_length_cylinder(tank_model, volume, tank_diameter):
    # V_total = V_cylinder + V_ellipsoidal_caps = π * r^2 * h + 2 * (π/24) * D^3
    if tank_model == "LOX":
        radius = tank_diameter / 2
        cap_volume = (np.pi / 24) * tank_diameter ** 3
        ellipsoidal_caps_volume = 2 * cap_volume
        cylinder_volume = volume - 2 * ellipsoidal_caps_volume
        length = cylinder_volume / (np.pi * radius ** 2)
    elif tank_model == "LH2":
        radius = tank_diameter / 2
        ellipsoidal_cap_volume = (np.pi / 24) * tank_diameter ** 3
        cylinder_volume = volume
        length = cylinder_volume / (np.pi * radius ** 2)
    return length

def calculate_tank_length(tank_model, volume, radius_ratio, phi):
    # R = ((V*3*tan(phi))/(pi * (1-radius_ratio) * (1+radius_ratio+radius_ratio^2)))^(1/3)
    # V_total = pi/3 * h (r^2 + r*R + R^2)
    phi = phi/180*np.pi
    if tank_model == "LOX":
        R = ((volume*3*np.tan(phi))/(np.pi * (1-radius_ratio) * (1+radius_ratio+radius_ratio**2)))**(1/3)
        r = radius_ratio * R
        length = ((1-radius_ratio) * R)/np.tan(phi)
    elif tank_model == "LH2":
        R = ((volume*3*np.tan(phi))/(np.pi * (1-radius_ratio) * (1+radius_ratio+radius_ratio**2)))**(1/3)
        r = radius_ratio * R
        length = ((1-radius_ratio) * R)/np.tan(phi)
    return length, R, r

def calculate_tank_thickness(wet_mass, propellant_pressure, fuel_mass, tank_length, R, r, phi, E, strength, thrust, gamma=0.65):
    phi = phi/180*np.pi
    max_thickness = 0.1  # 10 cm limit
    t = 0.001  # Start with 1 mm
    t_press = (propellant_pressure * R)/(np.cos(phi) * (strength * 0.85 - propellant_pressure/6894.76))
    while t < max_thickness:
        #axial stress due to the launch acceleration
        V = np.pi*r*tank_length*t
        sigma_axial = (thrust)/(2*np.pi*R*t*np.cos(phi))

        #bending stress due to lateral loads
        sigma_bend = (7800*V*2.2*9.81*(t**2/2))/(t**3/12)

        sigma_cr_axial = gamma * 0.605 * E * t / (R)
        p_cr = E/(1-0.3**2)*(t/R)**2
        interaction = safety_factor * (sigma_axial/sigma_cr_axial + sigma_bend/sigma_cr_axial + propellant_pressure/p_cr)

        if interaction <= 1.0:
            break
        t += 0.0001  # Increment 0.1 mm
    else:
        raise ValueError("Could not find a suitable tank thickness. Check your assumptions or inputs.")
    print("the axial stress is: " + str(sigma_axial))
    print("the bending stress is: " + str(sigma_bend))
    #check hoop stress
    if t_press>t:
        t = t_press
        print("Hoop stress leading")
    return t

def calculate_tank_mass(R, r, tank_length, thickness, material_density):
    volume_shell_wall = np.pi * thickness * (R + r) * ((R+r)**2 + tank_length**2)**(1/2)
    mass = volume_shell_wall * material_density
    return mass


######################################################################################################################

print("######################################################################################################")
print("Results: ")
# Choose a material
material = "Annealed 304L Stainless Steel"

# Access density and strength
density = materials_properties[material]["density"]       # kg/m^3
strength = materials_properties[material]["strength"]     # Pa
young_modulus = materials_properties[material]["young modulus"]   # Pa

tank_length_LH2, R_LH2, r_LH2 = calculate_tank_length("LH2", LH2_volume, radius_ratio, phi)
print("LH2 Tank Length: " + str(tank_length_LH2))
print("LH2 Tank Bottom Diameter: " + str(R_LH2*2))
print("LH2 Tank Top Diameter: " + str(r_LH2*2))
tank_length_LOX, R_LOX, r_LOX = calculate_tank_length("LOX", LOX_volume, radius_ratio, phi)
print("LOX Tank Length: " + str(tank_length_LOX))
print("LOX Tank Bottom Diameter: " + str(R_LOX * 2))
print("LOX Tank Top Diameter: " + str(r_LOX * 2))
print("Volume LOX: " + str(LOX_volume))

print("LH2 Tanks Volume: " + str(LH2_volume))
print("LOX Tanks Volume: " + str(LOX_volume))
thickness_LH2 = calculate_tank_thickness(wet_mass, LH2_pressure, LH2_mass, tank_length_LH2, R_LH2, r_LH2, phi, young_modulus, strength, thrust_engines, gamma=0.65)
print("Thickness LH2 Tank: " + str(thickness_LH2))
thickness_LOX = calculate_tank_thickness(wet_mass, LOX_pressure, LOX_mass, tank_length_LOX, R_LOX, r_LOX, phi, young_modulus, strength, thrust_engines, gamma=0.65)
print("Thickness LOX Tank: " + str(thickness_LOX))
mass_LH2_tank = calculate_tank_mass(R_LH2, r_LH2, tank_length_LH2, thickness_LH2, density)
print("Mass LH2 Tank: " + str(mass_LH2_tank))
mass_LOX_tank = calculate_tank_mass(R_LOX, r_LOX, tank_length_LOX, thickness_LOX, density)
print("Mass LOX Tank: " + str(mass_LOX_tank))