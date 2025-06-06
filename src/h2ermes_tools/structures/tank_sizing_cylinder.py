import numpy as np

materials_properties = {
    "Annealed 304L Stainless Steel": {
        "density": 7800,
        "strength": 1060 * 10**6,
        "young modulus": 200 * 10**9,
    },
}

# Constraints
safety_factor = 1.5
safety_factor_pressure = 1.5
tank_diameter = 7
payload_mass = 15000
thrust_engines = 2129573.909
LH2_pressure = 1000 * safety_factor_pressure  # kPa (10bar)
LOX_pressure = 270 * 2  # kPa (2.7 bar)
LH2_boiloff_margin = 1.1  # 1.012 * 1.02 #(3.2% extra)
LOX_boiloff_margin = 1.0143 * 1.02  # (3.46% extra)
LH2_density = 77  # kg/m3
LOX_density = 1141  # kg/m3
structural_mass = 28040.57293
propellant_mass = 203643.4588
wet_mass = propellant_mass + structural_mass
LH2_mass = 1 / 7.03 * propellant_mass + payload_mass
print("Mass LH2: " + str(LH2_mass) + " kg")
LOX_mass = 6.03 / 7.03 * propellant_mass
print("Mass LOX: " + str(LOX_mass) + " kg")
LH2_volume = LH2_mass / LH2_density * LH2_boiloff_margin
print("Volume LH2: " + str(LH2_volume) + " m3")
LOX_volume = LOX_mass / LOX_density * LOX_boiloff_margin
print("Volume LOX: " + str(LOX_volume) + " m3")
print(
    "-----------------------------------------------------------------------------------------------------------------"
)


def calculate_tank_length(tank_model, volume, tank_diameter):
    # V_total = V_cylinder + V_ellipsoidal_caps = π * r^2 * h + 2 * (π/24) * D^3
    if tank_model == "LOX":
        radius = tank_diameter / 2
        cap_volume = (np.pi / 24) * tank_diameter**3
        ellipsoidal_caps_volume = 2 * cap_volume
        cylinder_volume = volume - 2 * ellipsoidal_caps_volume
        length = cylinder_volume / (np.pi * radius**2)
    elif tank_model == "LH2":
        radius = tank_diameter / 2
        ellipsoidal_cap_volume = (np.pi / 24) * tank_diameter**3
        cylinder_volume = volume
        length = cylinder_volume / (np.pi * radius**2)
    return length


def calculate_tank_thickness(
    wet_mass,
    propellant_pressure,
    fuel_mass,
    tank_length,
    tank_diameter,
    E,
    strength,
    thrust,
    gamma=0.65,
):
    max_thickness = 0.1  # 10 cm limit
    t = 0.001  # Start with 1 mm
    while t < max_thickness:
        # pressure
        p_press = (
            2 * np.pi * E * t**2 * (0.605 * gamma)
            + propellant_pressure * np.pi * (tank_diameter / 2) ** 2
        )
        p_cr = (
            0.926
            * E
            * np.sqrt(gamma)
            / (
                ((tank_diameter / 2) / t) ** (5 / 2)
                * (tank_length / (tank_diameter / 2))
            )
        )
        # axial stress due to the launch acceleration
        sigma_axial = thrust_engines / (2 * np.pi * tank_diameter / 2 * t)
        # bending stress due to lateral loads
        sigma_bend = (
            np.pi * (tank_diameter / 2) * E * (t) ** 2 * (0.605 * gamma)
            + 0.8 * propellant_pressure * np.pi * (tank_diameter / 2) ** 3
        )
        sigma_cr_axial = gamma * 0.605 * E * t / (tank_diameter / 2)
        D = E * t**3 / (12 * (1 - 0.3**2))
        z = 0.95 * tank_length**2 / (tank_length / 2 * t)
        k_x = gamma * z * 4 * 1.73 / np.pi**2
        sigma_cr_bend = (
            k_x * (np.pi) ** 3 * D * (tank_diameter / 2) ** 2 / (tank_length) ** 2
        )
        interaction = safety_factor * (
            sigma_axial / sigma_cr_axial
            + sigma_bend / sigma_cr_bend
            + propellant_pressure / p_cr
        )
        if interaction <= 1.0:
            break
        t += 0.0001  # Increment 0.1 mm
    else:
        raise ValueError(
            "Could not find a suitable tank thickness. Check your assumptions or inputs."
        )
    print("the axial stress is: " + str(sigma_axial))
    print("the bending stress is: " + str(sigma_bend))
    # check hoop stress
    thickness_pressure = (propellant_pressure * 1000 * tank_diameter) / (2 * strength)
    if thickness_pressure > t:
        t = thickness_pressure
        print("Hoop stress leading")
    return t


def calculate_tank_mass(tank_diameter, tank_length, thickness, material_density):
    volume_shell = (
        2
        * np.pi
        * tank_diameter
        / 2
        * thickness
        * (tank_length + 2 * tank_diameter / 2)
    )
    mass = volume_shell * material_density
    return mass


######################################################################################################################
# Choose a material
material = "Annealed 304L Stainless Steel"
# Access density and strength
density = materials_properties[material]["density"]  # kg/m^3
strength = materials_properties[material]["strength"]  # Pa
young_modulus = materials_properties[material]["young modulus"]  # Pa
tank_length_LH2 = calculate_tank_length("LH2", LH2_volume, tank_diameter)
print("LH2 Tank Length: " + str(tank_length_LH2))
tank_length_LOX = calculate_tank_length("LOX", LOX_volume, tank_diameter)
print("LOX Tank Length: " + str(tank_length_LOX))
print("Volume LOX: " + str(LOX_volume))
print("LH2 Tanks Volume: " + str(LH2_volume))
print("LOX Tanks Volume: " + str(LOX_volume))
thickness_LH2 = calculate_tank_thickness(
    wet_mass,
    LH2_pressure,
    LH2_mass,
    tank_length_LH2,
    tank_diameter,
    young_modulus,
    strength,
    thrust_engines,
    gamma=0.65,
)
print("Thickness LH2 Tank: " + str(thickness_LH2))
thickness_LOX = calculate_tank_thickness(
    wet_mass,
    LOX_pressure,
    LOX_mass,
    tank_length_LOX,
    tank_diameter,
    young_modulus,
    strength,
    thrust_engines,
    gamma=0.65,
)
print("Thickness LOX Tank: " + str(thickness_LOX))
mass_LH2_tank = calculate_tank_mass(
    tank_diameter, tank_length_LH2, thickness_LH2, density
)
print("Mass LH2 Tank: " + str(mass_LH2_tank))
mass_LOX_tank = calculate_tank_mass(
    tank_diameter, tank_length_LOX, thickness_LOX, density
)
print("Mass LOX Tank: " + str(mass_LOX_tank))
