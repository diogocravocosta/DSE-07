import numpy as np
import matplotlib.pyplot as plt

def calculate_frustum_tank_length(volume: float,
                                  top_radius: float,
                                  bottom_radius: float,
                                  caps_height_radius_ratio: float,
                                  subtract_bottom_cap: bool,
                                  subtract_top_cap: bool
                                  ) -> float:
    """
    Calculate the length of a frustum tank given its volume, top and bottom radius,
    and the ratio of the height of the caps to the radius.
    Args:
        volume:
        top_radius:
        bottom_radius:
        caps_height_radius_ratio:

    Returns:
        The length of the frustum tank.
    """
    bottom_cap_volume = 2/3 * np.pi * bottom_radius**3 * caps_height_radius_ratio
    top_cap_volume = 2/3 * np.pi * top_radius**3 * caps_height_radius_ratio

    frustum_volume = volume

    if subtract_bottom_cap:
        frustum_volume -= bottom_cap_volume
    else:
        frustum_volume += bottom_cap_volume

    if subtract_top_cap:
        frustum_volume -= top_cap_volume
    else:
        frustum_volume += top_cap_volume

    if frustum_volume < 0:
        raise ValueError("Volume is too small for the given caps height ratio.")

    frustum_height = 3 * frustum_volume / (np.pi * (top_radius**2 + bottom_radius**2 + top_radius * bottom_radius))
    return frustum_height


def calculate_tank_thickness(
    wet_mass,
    propellant_pressure,
    max_propellant_pressure,
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
    """
    Compute required wall thickness of a conical tank using given rocket and tank parameters.
    """

    # Average radius (approximate effective radius for stress)
    radius = (R + r) / 2
    # Allowable stress
    sigma_allow = strength / gamma

    def von_mises(t):
        # Hoop stress due to internal pressure
        sigma_h = (propellant_pressure * radius / (t)) * ((1 - np.sin(phi)**2) / np.sin(phi))
        # Meridional (axial) stress due to internal pressure
        sigma_m_p = (propellant_pressure * radius) / (2*t * np.sin(phi))
        # Meridional (axial) stress due to thrust
        A = 2 * np.pi * radius * t  # conical shell approx area
        sigma_m_l = safety_factor * thrust / A
        # Meridional (bending) stress 
        #sigma_b = M_bend * y / I
        sigma_m = sigma_m_p + sigma_m_l
        # Von Mises stress
        return np.sqrt(sigma_h**2 + sigma_m**2 - sigma_h * sigma_m)

    # Iterate to find minimum thickness
    t_guess = 0.001  # start at 1 mm
    while True:
        sigma_vm = von_mises(t_guess)
        if sigma_vm <= sigma_allow:
            break
        t_guess += 0.0001  # increment in 0.1 mm steps

    #Check for Axial Buckling 
    N_cr = 0.33 * (2*np.pi*E*t_guess**2*np.cos(phi))/np.sqrt(3*(1-0.33**3))
    if safety_factor * thrust / (2 * np.pi * radius * t_guess) > N_cr:
        print("THE STRUCTURE BUCKLES and sigma critical is: " + str(N_cr))
    
    #Check for hoop stress at maximum boil off pressure:
    t_press = (max_propellant_pressure * radius)/(strength) * (1-np.sin(phi)**2)/(np.sin(phi))
    # check hoop stress
    if t_press > t_guess:
        t_guess = t_press
        print("Hoop stress leading and hoop stress is: " + str((max_propellant_pressure * radius)/(t_press) * (1-np.sin(phi)**2)/(np.sin(phi))))
    return t_guess


def calculate_tank_mass(bottom_radius: float,
                        top_radius: float,
                        tank_length: float,
                        thickness: float,
                        material_density: float,
                        cap_height_radius_ratio: float) -> float:

    volume_shell_wall = (
        (np.pi * tank_length / 3)
        * (bottom_radius**2 + bottom_radius*top_radius+top_radius**2
           - (bottom_radius-thickness)**2 - (bottom_radius-thickness)*(top_radius-thickness) - (top_radius-thickness)**2)
        + 2*np.pi/3*((bottom_radius+thickness)**3 - bottom_radius**3) * cap_height_radius_ratio
        + 2*np.pi/3*((top_radius+thickness)**3 -top_radius**3) * cap_height_radius_ratio
    )
    mass = volume_shell_wall * material_density
    return mass


def check_vibrations(tank_mass, thickness, E, tank_length):
    omega = np.sqrt(
        (3 * E * (np.pi * (2 * thickness) / 64) / (tank_mass * tank_length**3))
    )  # with k = 3EI/L^3 modelled as cantilever beam
    f_natural = 1 / (2 * np.pi) * omega
    return f_natural


def plot_stress_vs_thickness(P, R, r, phi, thrust, strength, gamma):
    thicknesses = np.linspace(0.001, 0.1, 500)
    stresses = []

    for t in thicknesses:
        radius = (R + r) / 2
        sigma_h = (P * radius / (t)) * ((2 - np.sin(phi)**2) / np.sin(phi))
        sigma_m_p = (P * radius) / (t * np.sin(phi))
        A = 2 * np.pi * radius * t
        sigma_m_l = thrust / A
        sigma_m = sigma_m_p + sigma_m_l
        sigma_vm = np.sqrt(sigma_h**2 + sigma_m**2 - sigma_h * sigma_m)
        stresses.append(sigma_vm / 1e6)  # MPa

    plt.plot(thicknesses * 1000, stresses)
    plt.axhline(y=strength/gamma / 1e6, color='r', linestyle='--', label='Allowable')
    plt.xlabel("Wall Thickness (mm)")
    plt.ylabel("Von Mises Stress (MPa)")
    plt.title("Stress vs Thickness")
    plt.grid(True)
    plt.legend()
    plt.show()

# def tank_overall_dimensions():
#     ratio_radius = 0.5
#     phi = np.arange(0, 91, 1)
#     R = (
#         (520 + 99)
#         * 3
#         * np.tan(np.radians(phi))
#         / (np.pi * (1 - ratio_radius) * (1 + ratio_radius + ratio_radius**2))
#     ) * (1 / 3)
#     r = R * ratio_radius
#     h = (R - r) / np.tan(np.radians(phi))
#     sols = 0
#     for i, val in enumerate(R):
#         if 2 * val >= 10:
#             continue
#         if 2 * val < 7:
#             continue
#         sols += 1
#         print(f"For {ratio_radius=}, {phi[i]=}:\n{R[i]=}\n{r[i]=}\n{h[i]=}\n")
#     print(sols)

######################################################################################################################
if __name__ == "__main__":
    materials_properties = {
        "Annealed 304L Stainless Steel": {
            "density"      : 7800,
            "strength"     : 1060 * 10 ** 6,
            "young modulus": 200 * 10 ** 9,
        },
    }
    # Geometry
    radius_ratio = 0.5
    top_radius_ratio = 0.52
    bottom_radius_ratio = radius_ratio / top_radius_ratio

    bottom_radius = 5  # m
    middle_radius = bottom_radius * bottom_radius_ratio
    top_radius = middle_radius * top_radius_ratio
    cap_height_radius_ratio = 0.5  # ratio of the height of the caps to the radius
    phi = np.deg2rad(10)  # phi in radians

    # Constraints
    safety_factor = 1.5
    safety_factor_pressure = 2.0
    # tank_diameter = 7
    payload_mass = 15000
    LH2_pressure = 200 * safety_factor_pressure * 1000  # Pa (10bar)
    max_pressure_boiloff = 200 * 1000 * safety_factor_pressure
    LOX_pressure = 280 * safety_factor_pressure * 1000  # Pa (2.7 bar)

    LH2_ullage_margin = 1.1  # 10% ullage
    LOX_ullage_margin = 1.1  # 10% ullage

    LH2_density = 77  # kg/m3
    LOX_density = 1340  # kg/m3

    structural_mass = 20642.21346
    propellant_mass = 149913.1903
    wet_mass = propellant_mass + structural_mass
    thrust_engines = 1.2 * wet_mass
    LH2_mass = 1 / 7.0 * propellant_mass + payload_mass
    print("Mass LH2: " + str(LH2_mass) + " kg")
    LOX_mass = 6.0 / 7.0 * propellant_mass
    print("Mass LOX: " + str(LOX_mass) + " kg")
    LH2_volume = LH2_mass / LH2_density * LH2_ullage_margin
    print("Volume LH2: " + str(LH2_volume) + " m3")
    LOX_volume = LOX_mass / LOX_density  * LOX_ullage_margin
    print("Volume LOX: " + str(LOX_volume) + " m3")
    print(
        "-----------------------------------------------------------------------------------------------------------------"
    )

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

    tank_length_LH2 = calculate_frustum_tank_length(LH2_volume,
                                                    top_radius,
                                                    middle_radius,
                                                    cap_height_radius_ratio,
                                                    subtract_bottom_cap=False,
                                                    subtract_top_cap=False)
    print("LH2 Tank Length: " + str(tank_length_LH2) + " m")
    print("LH2 Tank Bottom Diameter: " + str(middle_radius * 2) + " m")
    print("LH2 Tank Top Diameter: " + str(top_radius * 2) + " m")

    print(
        "----------------------------------------------------------------------------------------------"
    )
    tank_length_LOX = calculate_frustum_tank_length(LOX_volume,
                                                    middle_radius,
                                                    bottom_radius,
                                                    cap_height_radius_ratio,
                                                    subtract_bottom_cap=False,
                                                    subtract_top_cap=True)
    print("LOX Tank Length: " + str(tank_length_LOX) + " m")
    print("LOX Tank Bottom Diameter: " + str(bottom_radius * 2) + " m")
    print("LOX Tank Top Diameter: " + str(middle_radius * 2) + " m")
    print(
        "----------------------------------------------------------------------------------------------"
    )
    print("LH2 Tanks Volume: " + str(LH2_volume) + " m^3")
    print("LOX Tanks Volume: " + str(LOX_volume) + " m^3")
    thickness_LH2 = calculate_tank_thickness(
        wet_mass,
        LH2_pressure,
        max_pressure_boiloff,
        LH2_mass,
        tank_length_LH2,
        bottom_radius,
        middle_radius,
        phi,
        young_modulus,
        strength,
        thrust_engines,
        gamma=0.65,
    )
    print("Thickness LH2 Tank: " + str(thickness_LH2) + " m")
    thickness_LOX = calculate_tank_thickness(
        wet_mass,
        LOX_pressure,
        LOX_pressure,
        LOX_mass,
        tank_length_LOX,
        middle_radius,
        top_radius,
        phi,
        young_modulus,
        strength,
        thrust_engines,
        gamma=0.65,
    )
    print("Thickness LOX Tank: " + str(thickness_LOX) + " m")
    mass_LH2_tank = calculate_tank_mass(
        bottom_radius, middle_radius, tank_length_LH2, thickness_LH2, density, cap_height_radius_ratio
    )
    print("Mass LH2 Tank: " + str(mass_LH2_tank) + " kg")
    mass_LOX_tank = calculate_tank_mass(
        middle_radius, top_radius, tank_length_LOX, thickness_LOX, density, cap_height_radius_ratio
    )
    print("Mass LOX Tank: " + str(mass_LOX_tank) + " kg")

    # tank_overall_dimensions()
    natural_frequency = check_vibrations(
        mass_LH2_tank, thickness_LH2, young_modulus, tank_length_LH2
    )
    print("Natural Frequency: " + str(natural_frequency) + " Hz")
    plot_stress_vs_thickness(LOX_pressure, bottom_radius, middle_radius, phi, thrust_engines, strength, gamma=0.65)
    plot_stress_vs_thickness(LH2_pressure, middle_radius, top_radius, phi, thrust_engines, strength, gamma=0.65)