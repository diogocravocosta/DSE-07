import numpy as np

materials_properties = {
    "Annealed 304L Stainless Steel": {
        "density": 7800,
        "strength": 1060 * 10**6,
        "young modulus": 200 * 10**9,
    },
}
# Geometry
R = 5  # m
radius_ratio = 0.5
phi = np.deg2rad(10)  # phi in radians

# Constraints
safety_factor = 1.5
safety_factor_pressure = 2.0
tank_diameter = 7
payload_mass = 15000
thrust_engines = 2129573.909

LH2_pressure = 1000 * safety_factor_pressure * 1000  # Pa (10bar)
LOX_pressure = 270 * 2 * 1000  # Pa (2.7 bar)

LH2_boiloff_margin = 1.1  # 10% ullage
LOX_boiloff_margin = 1.0143 * 1.02  # 3% ullage

LH2_density = 77  # kg/m3
LOX_density = 1340  # kg/m3

structural_mass = 20642.21346
propellant_mass = 149913.1903
wet_mass = propellant_mass + structural_mass
LH2_mass = 1 / 7.0 * propellant_mass + payload_mass
print("Mass LH2: " + str(LH2_mass) + " kg")
LOX_mass = 6.0 / 7.0 * propellant_mass
print("Mass LOX: " + str(LOX_mass) + " kg")
LH2_volume = LH2_mass / LH2_density * LH2_boiloff_margin
print("Volume LH2: " + str(LH2_volume) + " m3")
LOX_volume = LOX_mass / LOX_density * LOX_boiloff_margin
print("Volume LOX: " + str(LOX_volume) + " m3")
print(
    "-----------------------------------------------------------------------------------------------------------------"
)


def tank_overall_dimensions():
    ratio_radius = 0.5
    phi = np.arange(0, 91, 1)
    R = (
        (520 + 99)
        * 3
        * np.tan(np.radians(phi))
        / (np.pi * (1 - ratio_radius) * (1 + ratio_radius + ratio_radius**2))
    ) * (1 / 3)
    r = R * ratio_radius
    h = (R - r) / np.tan(np.radians(phi))
    sols = 0
    for i, val in enumerate(R):
        if 2 * val >= 10:
            continue
        if 2 * val < 7:
            continue
        sols += 1
        print(f"For {ratio_radius=}, {phi[i]=}:\n{R[i]=}\n{r[i]=}\n{h[i]=}\n")
    print(sols)


def calculate_tank_length_LOX(tank_model, volume, radius_ratio, phi, R):
    r = radius_ratio * R
    h = (3 * volume - 2 * np.pi * R + 2 * np.pi * r) / (np.pi * (R**2 + r * R + r**2))
    # total_length = h + R
    return h, R, r


def calculate_tank_length_LH2(tank_model, volume, radius_ratio, phi, LOX_radius):
    R = LOX_radius
    r = radius_ratio * R
    h = (3 * volume - 2 * np.pi * r**3 - 2 * np.pi * R**3) / (
        np.pi * (R**2 + r * R + r**2)
    )
    # total_length = h + r + R
    return h, R, r


def calculate_tank_thickness(
    wet_mass,
    propellant_pressure,
    fuel_mass,
    tank_length,
    R,
    r,
    phi,
    E,
    strength,
    thrust,
    gamma=0.65,
):
    average_radius = (R + r) / (2 * np.cos(phi))
    max_thickness = 0.1  # 10 cm limit
    t = 0.001  # Start with 1 mm

    t_press = (propellant_pressure / 6894.76 * R) / (
        np.cos(phi) * (strength * 0.85 - propellant_pressure / 6894.76)
    )  # Pressure vessel design

    while t < max_thickness:
        # axial stress due to the launch acceleration
        V = np.pi * r * tank_length * t
        # sigma_axial = (thrust) / (2 * np.pi * R * t * np.cos(phi))  # Roak textbook

        # bending stress due to lateral loads
        # sigma_bend = (7800 * V * 2.2 * 9.81 * (t**2 / 2)) / (t**3 / 12)

        # Combined Loading
        N_axial = np.cos(phi) * thrust - np.sin(phi) * 7800 * V * 2 * 9.81
        M_bend = (
            7800 * V * (np.sin(phi) * 6 + np.cos(phi) * 2) * 9.81 * tank_length / 2
        )  # (t**2 / 2)

        M_cr_bend = (0.33 / (3 * (1 - 0.33**2)) + 0.1) * (
            2 * np.pi * E * r * t**2 * np.cos(phi) ** 2
        ) + np.pi * r**3 * propellant_pressure / 2

        P_cr_axial = (0.33 / (3 * (1 - 0.33**2)) + 0.1) * (
            2 * np.pi * E * t**2 * np.cos(phi) ** 2
        ) + np.pi * r**2 * propellant_pressure

        p_cr = (0.92 * E * 0.75) / (
            (tank_length / average_radius) * (average_radius / t) ** (5 / 2)
        )  # in N
        interaction = safety_factor * (N_axial / P_cr_axial + M_bend / M_cr_bend)

        if interaction <= 1.0:
            break
        t += 0.0001  # Increment 0.1 mm
    else:
        raise ValueError(
            "Could not find a suitable tank thickness. Check your assumptions or inputs."
        )
    print("the axial load is: " + str(N_axial))
    print("the bending load is: " + str(M_bend))
    print("The PRESSURE IS: " + str(p_cr))
    # check hoop stress
    if t_press > t:
        t = t_press
        print("Hoop stress leading")
    return t


def calculate_tank_mass(R, r, tank_length, thickness, material_density):
    volume_shell_wall = (
        np.pi * thickness * (R + r) * ((R + r) ** 2 + tank_length**2) ** (1 / 2)
    )
    mass = volume_shell_wall * material_density
    return mass


def check_vibrations(tank_mass, thickness, E, tank_length):
    omega = np.sqrt(
        (3 * E * (np.pi * (2 * thickness) / 64) / (tank_mass * tank_length**3))
    )  # with k = 3EI/L^3 modelled as cantiliver beam
    f_natural = 1 / (2 * np.pi) * omega
    print(f_natural)
    return f_natural


######################################################################################################################

print(
    "######################################################################################################"
)
print("Results: ")
# Choose a material
material = "Annealed 304L Stainless Steel"

# Access density and strength
density = materials_properties[material]["density"]  # kg/m^3
strength = materials_properties[material]["strength"]  # Pa
young_modulus = materials_properties[material]["young modulus"]  # Pa

tank_length_LOX, R_LOX, r_LOX = calculate_tank_length_LOX(
    "LOX", LOX_volume, 0.9, phi, R
)
print("LOX Tank Length: " + str(tank_length_LOX))
print("LOX Tank Bottom Diameter: " + str(R_LOX * 2))
print("LOX Tank Top Diameter: " + str(r_LOX * 2))
print(
    "----------------------------------------------------------------------------------------------"
)
tank_length_LH2, R_LH2, r_LH2 = calculate_tank_length_LH2(
    "LH2", LH2_volume, radius_ratio, phi, r_LOX
)
print("LH2 Tank Length: " + str(tank_length_LH2))
print("LH2 Tank Bottom Diameter: " + str(R_LH2 * 2))
print("LH2 Tank Top Diameter: " + str(r_LH2 * 2))
print(
    "----------------------------------------------------------------------------------------------"
)
print("LH2 Tanks Volume: " + str(LH2_volume))
print("LOX Tanks Volume: " + str(LOX_volume))
thickness_LH2 = calculate_tank_thickness(
    wet_mass,
    LH2_pressure,
    LH2_mass,
    tank_length_LH2,
    R_LH2,
    r_LH2,
    phi,
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
    R_LOX,
    r_LOX,
    phi,
    young_modulus,
    strength,
    thrust_engines,
    gamma=0.65,
)
print("Thickness LOX Tank: " + str(thickness_LOX))
mass_LH2_tank = calculate_tank_mass(
    R_LH2, r_LH2, tank_length_LH2, thickness_LH2, density
)
print("Mass LH2 Tank: " + str(mass_LH2_tank))
mass_LOX_tank = calculate_tank_mass(
    R_LOX, r_LOX, tank_length_LOX, thickness_LOX, density
)
print("Mass LOX Tank: " + str(mass_LOX_tank))

# tank_overall_dimensions()
natural_frequency = check_vibrations(
    mass_LH2_tank, thickness_LH2, young_modulus, tank_length_LH2
)
print("Natural Frequency: " + str(natural_frequency))
