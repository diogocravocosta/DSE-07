import numpy as np
import matplotlib.pyplot as plt
from data import constants as cs
from data import material as mat
import rainflow
from scipy.stats import linregress


#--------------------------------------------------------
#Functions for Stress and Crack Growth Calculations
#--------------------------------------------------------

def mechanical_stress(force, radius, thickness, phi):
    """
    Calculates the mechanical stress in the wall of the conical tank for an axial load.
    Takes into account force, radius, thickness and degree of curvature.
    """
    sigma_conical = force / (2 * np.pi * radius * thickness * np.cos(np.radians(phi)))
    return sigma_conical


def thermal_stress(delta_T, young_mod, thermal_expansion_coeff, phi):
    """
    Calculates the thermal stress in the tank by assuming that the tank is clamped at both ends
    and then uses mechanical stress formula to see effect of thermal induced stress.
    """
    stress = delta_T * young_mod * thermal_expansion_coeff
    thermal_stress = stress / np.cos(np.radians(phi))
    return thermal_stress

    #def pressure_stress(Pressure_vapor,radius_tank,thickness_tank, phi):
    """
    Calculates the stress caused by the pressure inside a conical tank. 
    """
    # Conical stress formula
    conical_stress = Pressure_vapor / 0.85 * (
            radius_tank / thickness_tank / np.cos(np.radians(phi)) / 0.85 + 1 / 6894.76)
    return conical_stress


def nasa_pressure_combined_load_stress(material, thickness_tank, phi, radius_small, pressure_tank):
    """
    Formula to calculate the critical load in conical tank formulated by NASA and then getting stress in walls by using mechanical stress formula
    """
    gamma = 0.33
    delta_gamma = 0.12
    mu = 0.33  # GET EXACT VALUE
    term1 = gamma / np.sqrt(3 * (1 - mu ** 2))
    stiffness_term = (term1 + delta_gamma) * (
            2 * np.pi * material.E * thickness_tank ** 2 * np.cos(np.radians(phi)) ** 2)
    pressure_term = np.pi * radius_small ** 2 * pressure_tank
    P_cr = stiffness_term + pressure_term
    stress = mechanical_stress(P_cr, radius_small, thickness_tank, phi)
    return stress


def crack_growth(paris_coeff_C, paris_exp_m, sigma_max, sigma_min, Y_geometry_factor, a_crack_depth):
    """
    Uses the paris equation model crack growth for a sigma max and sigma min for a specified crack depth a. 
    """
    delta_K = Y_geometry_factor * (sigma_max - sigma_min) / 10 ** 6 * (np.pi * a_crack_depth) ** 0.5

    da_dn = paris_coeff_C * (delta_K ** paris_exp_m)
    return da_dn


def R_calculation(loading):
    """
    Calculates the load ratio R for a complex loading series.
    """
    sigma_max = max(loading)
    s = []
    for i in range(len(loading)):
        if loading[i] != 0:
            s.append(loading[i])
    sigma_min = min(s)
    R_load_ratio = sigma_min / sigma_max
    return R_load_ratio


def t_critical(yield_stress, pressure_vapor, phi, radius_tank):
    '''
    Calculates minimum conical tank wall thickness to not fail for a specified pressure load.
    '''
    thickness_crit = pressure_vapor * radius_tank / (
            np.cos(np.radians(phi)) * (yield_stress * 0.85 - pressure_vapor / 6894.76))
    return thickness_crit


def y_calculation(a_crack_depth, thickness):
    """
    Calculates the Y factor that goes into paris equation. Essentially is a geometry factor.
    """
    return 1 + 1.5 * (a_crack_depth / thickness) ** 2


def miners(stress_range, miner_coeff_C, miner_coeff_m, stress_cycle, safety_factor, launches):
    """
    Calculates if the material is safe from fracture using the Miner rule for fracture. Takes into account stress range and cycle calculated with rainflow.
    """
    sigma = cycle_launch(stress_range, safety_factor)
    cycle = cycle_launch(stress_cycle, launches)

    Di = []
    for i in range(len(sigma)):
        damage_calc = cycle[i] * sigma[i] ** miner_coeff_m / miner_coeff_C
        Di.append(damage_calc)

    damage = sum(Di)

    return Di, cycle, damage


def cycle_launch(stress_cycle, launches):
    """
    Multiplies a given vector (stress_cycle) by a constant (launches)
    """
    cycle = np.zeros(len(stress_cycle))
    for i in range(len(stress_cycle)):
        cycle[i] = stress_cycle[i] * launches
    return cycle


def critical_crack_depth_calc(material, Y_geometry_factor, sigma_max):
    """
    Calculates the thickness at which a crack will propagate by taking into account max applied sigma and material toughness.
    """
    critical_crack_depth = (material.Kic / (Y_geometry_factor * sigma_max * 1e-6)) ** 2 / np.pi
    return critical_crack_depth


def fatigue_paris_estimation(
        a_crack_depth,
        thickness,
        count,
        sigma_global_loading,
        stress_range,
        stress_cycle,
        launches,
        paris_coeff_C,
        paris_exp_m,
        plot,
        material,
        min_launches):
    """
    Calculates the fatigue life due to paris equation.
    Also takes into account the thickness at which failure happens due to 3 separate loading conditions.
    """
    Y_geometry_factor = y_calculation(a_crack_depth, thickness)

    critical_crack_depth = critical_crack_depth_calc(material, Y_geometry_factor, max(sigma_global_loading))
    # print('Critical crack depth is:', critical_crack_depth, "m")
    stress_cycle = cycle_launch(stress_cycle, launches)
    stress_range = cycle_launch(stress_range, 1e6)  # convert to pascals
    a_crack = []
    for i in range(len(stress_range)):
        for j in range(int(stress_cycle[i])):
            Y_geometry_factor = y_calculation(a_crack_depth,
                                              thickness)  # Calculate Y factor based on current crack size
            da_dn = crack_growth(paris_coeff_C,
                                 paris_exp_m,
                                 stress_range[i],
                                 0,
                                 Y_geometry_factor,
                                 a_crack_depth)

            a_crack_depth = a_crack_depth + da_dn
            a_crack.append(a_crack_depth)
            count += 1
            if a_crack_depth > critical_crack_depth:
                # print('Expected lifecycles: ', count)
                if count < min_launches:
                    thickness = thickness + 0.001
                break

    if plot == True:
        cycle = np.arange(0, count, 1)
        plt.figure(figsize=(8, 5))
        plt.plot(cycle, a_crack, marker='o')
        plt.xlabel('Cycle')
        plt.ylabel('Crack Size (mm)')
        plt.title('Crack Growth vs Cycle')
        plt.grid(True)
        plt.tight_layout()
        plt.show()
    return count


def fatigue_miner_estimation(stress_range, miner_c, miner_m, stress_cycle, safety_factor, launches, plot):
    """
    Calculates the fatigue life due to Miners cycle.
    Takes into account a safety factor and a conservative number of launches.
    """
    Di, cycle, damage = miners(stress_range, miner_c, miner_m, stress_cycle, safety_factor, launches)
    if damage > 1:
        pass
        # print('Failure was predicted due to Miners rule.')
    while damage > 1:
        launches = launches - 1
        damage = miners(stress_range, miner_c, miner_m, stress_cycle, safety_factor, launches)[2]
    # print('Damage count is:', damage, ". Expected number of launches:", launches)
    if plot == True:
        plt.figure(figsize=(8, 5))
        plt.plot(cycle, Di, marker='o')
        plt.xlabel('Cycle')
        plt.ylabel('Damage per Cycle (Di)')
        plt.title('Miner\'s Rule: Damage per Cycle vs Cycle')
        plt.grid(True)
        plt.tight_layout()
        plt.show()
    return damage, launches


def loading_phases(delta_T_earth,
                   delta_T_tank,
                   delta_T_tank_extreme,
                   material,
                   phi,
                   thickness_tank,
                   tank_radius,
                   pressure_launch_tank,
                   pressure_after_dock,
                   launch_mass,
                   cs,
                   g_launch_force_ratio,
                   thrust_engines,
                   pressure_vent,
                   pressure_dock_vent,
                   g_reentry_force_ratio,
                   dry_mass,
                   time_mission,
                   plot):
    """
    Calculates the loading conditions for different phases in the mission.
    Also calculates the individual load ratios for the different conditions (thermal, pressure, mechanical).
    """
    #------------------------------------------------------------
    # Different phases of mission and their loading conditions
    #------------------------------------------------------------
    # Before launch
    sigma_thermal_start = thermal_stress(delta_T_earth, material.E, material.cte, phi)  # is valid for this.
    sigma_pressure_start = nasa_pressure_combined_load_stress(material, thickness_tank, phi, tank_radius,
                                                              pressure_launch_tank)
    sigma_axial_start = mechanical_stress(launch_mass * cs.g_0, tank_radius, thickness_tank, phi)
    stress_comb_start = sigma_thermal_start + sigma_pressure_start + sigma_axial_start

    # 1st stage max q
    sigma_thermal_maxq = 0  #thermal_stress(delta_T_earth, material.E, material.cte,tank_radius, t, phi) # is valid for this.
    sigma_pressure_maxq = nasa_pressure_combined_load_stress(material, thickness_tank, phi, tank_radius,
                                                             pressure_launch_tank)
    sigma_axial_maxq = mechanical_stress(launch_mass * g_launch_force_ratio * cs.g_0, tank_radius, thickness_tank, phi)
    stress_comb_maxq = sigma_thermal_maxq + sigma_pressure_maxq + sigma_axial_maxq

    # During second stage fire
    sigma_thermal_launch = 0  #thermal_stress(delta_T_earth, material.E, material.cte,tank_radius, t, phi)
    sigma_pressure_launch = nasa_pressure_combined_load_stress(material, thickness_tank, phi, tank_radius,
                                                               pressure_launch_tank)
    sigma_axial_launch = mechanical_stress(thrust_engines, tank_radius, thickness_tank, phi)
    stress_comb_launch = sigma_thermal_launch + sigma_pressure_launch + sigma_axial_launch

    # Orbit - Right before docking
    sigma_thermal_orbit = thermal_stress(delta_T_tank, material.E, material.cte, phi)
    sigma_pressure_orbit = nasa_pressure_combined_load_stress(material, thickness_tank, phi, tank_radius, pressure_vent)
    stress_comb_orbit = sigma_thermal_orbit + sigma_pressure_orbit

    # Orbit - After docking
    sigma_thermal_dock = thermal_stress(delta_T_tank, material.E, material.cte, phi)
    sigma_pressure_dock = nasa_pressure_combined_load_stress(material, thickness_tank, phi, tank_radius,
                                                             pressure_after_dock)
    stress_comb_dock = sigma_thermal_orbit + sigma_pressure_orbit

    # Worst point during re-entry
    sigma_thermal_reentry = thermal_stress(delta_T_tank_extreme, material.E, material.cte, phi)
    sigma_axial_reentry = mechanical_stress(g_reentry_force_ratio * cs.g_0 * dry_mass, tank_radius, thickness_tank, phi)
    sigma_pressure_reentry = nasa_pressure_combined_load_stress(material, thickness_tank, phi, tank_radius,
                                                                pressure_vent)
    stress_comb_reentry = sigma_thermal_reentry + sigma_pressure_reentry + sigma_axial_reentry

    # On launch pad
    sigma_thermal_launchpad = thermal_stress(delta_T_earth, material.E, material.cte, phi)
    sigma_axial_launchpad = mechanical_stress(dry_mass * cs.g_0, tank_radius, thickness_tank, phi)
    stress_comb_launchpad = sigma_thermal_launchpad + sigma_axial_launchpad

    #------------------------------------------------------------
    # Calculate load ratio's R
    #------------------------------------------------------------

    sigma_global_loading = [stress_comb_start, stress_comb_maxq, stress_comb_launch, stress_comb_orbit,
                            stress_comb_dock, stress_comb_reentry, stress_comb_launchpad]
    R_load_ratio_global = R_calculation(sigma_global_loading)
    # print("Global Load ratio R = ", R_load_ratio_global)

    #R for thermal stress 
    sigma_thermal_load = [sigma_thermal_start, sigma_thermal_maxq, sigma_thermal_launch, sigma_thermal_orbit,
                          sigma_thermal_dock, sigma_thermal_reentry, sigma_thermal_launchpad]
    R_load_ratio_thermal = R_calculation(sigma_thermal_load)
    ## print("Maximum stress ratio R thermal: ", R_load_ratio_thermal)

    #R for mechanical stress
    sigma_mechanical_load = [sigma_axial_start, sigma_axial_maxq, sigma_axial_launch, 0, 0, sigma_axial_reentry,
                             sigma_axial_launchpad]
    R_load_ratio_mechanical = R_calculation(sigma_mechanical_load)
    ## print("Maximum stress ratio R mechanical: ", R_load_ratio_mechanical)

    #R for pressure stress
    sigma_pressure_load = [sigma_pressure_start, sigma_pressure_maxq, sigma_pressure_launch, sigma_pressure_orbit,
                           sigma_pressure_dock, sigma_pressure_reentry, 0]
    R_load_ratio_pressure = R_calculation(sigma_pressure_load)
    ## print("Maximum stress ratio R pressure: ", R_load_ratio_pressure)
    if plot is True:
        plt.figure(figsize=(8, 5))
        plt.plot(time_mission, np.array(sigma_pressure_load) / 1e6, marker='o', label='Pressure Stress')
        plt.plot(time_mission, np.array(sigma_thermal_load) / 1e6, marker='s', label='Thermal Stress')
        plt.plot(time_mission, np.array(sigma_mechanical_load) / 1e6, marker='^', label='Mechanical Stress')
        plt.plot(time_mission, np.array(sigma_global_loading) / 1e6, marker='x', label='Total Stress')
        plt.xlabel('Mission Time (hours)')
        plt.ylabel('Stress (MPa)')
        plt.title('Pressure, Thermal, Mechanical and Total Stress vs Mission Time')
        plt.legend()
        plt.grid(True)
        plt.tight_layout()
        plt.show()

    return sigma_global_loading


def rainflow_counting(loading_series):
    """
    Performs rainflow counting on a loading series with a time history.
    Returns:
        stress_range: array of stress ranges (MPa)
        stress_cycle: array of number of cycles for each range
    """
    # Use rainflow.extract_cycles to get cycles (assumes loading_series is in MPa)
    cycles = rainflow.extract_cycles(loading_series)
    stress_range = []
    stress_cycle = []

    # Each cycle is (range, mean, count, start_index, end_index)
    for rng, mean, count, istart, iend in cycles:
        stress_range.append(rng)
        stress_cycle.append(count)

    return stress_range, stress_cycle


def thickness_optimization_fatigue(phi,
                                   tank_radius,
                                   thickness_tank,
                                   material,
                                   fuel_reentry_LH2,
                                   dry_mass,
                                   launch_mass,
                                   payload_mass,
                                   g_reentry_force_ratio,
                                   g_launch_force_ratio,
                                   max_thrust2weight,
                                   number_of_launches=40,
                                   min_launches=25,
                                   safety_factor=2,
                                   time_mission=[0, 0.1, 0.3, 18, 21, 23, 24],
                                   a_crack_depth=0.001,
                                   pressure_launch_tank=3e5,
                                   pressure_dock_vent=5e5,
                                   pressure_vent=1e6,
                                   pressure_after_dock=4.5e5,
                                   T_lh2=20,
                                   T_space=4,
                                   T_ambient_earth=300,
                                   T_gh2=150,
                                   T_gh2_ext=200
                                   ):
    """

    Args:
        T_gh2_ext:
        T_gh2:
        T_ambient_earth:
        T_space:
        T_lh2:
        pressure_after_dock:
        pressure_vent:
        pressure_dock_vent:
        pressure_launch_tank:
        a_crack_depth:
        max_thrust2weight:
        g_launch_force_ratio:
        g_reentry_force_ratio:
        phi: from tank sizing
        tank_radius: from tank sizing
        thickness_tank: from tank sizing
        material:
        fuel_reentry_LH2: from mass integration
        dry_mass: updated from mass integration
        launch_mass: updated from mass integration
        payload_mass: updated from mass integration
        number_of_launches
        min_launches
        safety_factor


    Returns:

    """

    force_launch = g_launch_force_ratio * cs.g_0 * launch_mass
    thrust_engines = max_thrust2weight * (
            dry_mass + payload_mass) * cs.g_0  # N, thrust of engines, later import from dictionary in main branch.

    # Crack Growth Parameters
    da_dn = []  # initialize da/dn as an empty list
    a_crack = []  # to store crack growth
    a_crack_depth = a_crack_depth  # initial crack size in m
    count = 0
    Damage = 0

    delta_T_earth = T_ambient_earth - T_lh2
    delta_T_tank = T_gh2 - T_lh2
    delta_T_tank_extreme = T_gh2_ext - T_lh2

    # Paris law
    paris_coeff_C = 5.131e-17 / 1e3  # convert from mm/cycle to m/cycle
    paris_exp_m = 7.02  # dimensionless

    # MAT200 paper --- 1884, -0.1555
    const = 1884
    exp_coeff = -0.1555
    miner_m = -1 / exp_coeff
    miner_c = 10 ** (np.log10(const) / -exp_coeff)
    # print("The C coefficient for miner equation is", miner_c, "and m coefficient is", miner_m)
    # Conditions
    plot = False
    crack_cond = True

    sigma_global_loading = loading_phases(delta_T_earth,
                                          delta_T_tank,
                                          delta_T_tank_extreme,
                                          material,
                                          phi,
                                          thickness_tank,
                                          tank_radius,
                                          pressure_launch_tank,
                                          pressure_after_dock,
                                          launch_mass,
                                          cs,
                                          g_launch_force_ratio,
                                          thrust_engines,
                                          pressure_vent,
                                          pressure_dock_vent,
                                          g_reentry_force_ratio,
                                          dry_mass,
                                          time_mission,
                                          plot)
    stress_range, stress_cycle = rainflow_counting(cycle_launch(sigma_global_loading, 1e-6))  # Result in Mpa

    # with paris equation - crack growth
    count = fatigue_paris_estimation(
        a_crack_depth,
        thickness_tank,
        count,
        sigma_global_loading,
        stress_range,
        stress_cycle,
        number_of_launches,
        paris_coeff_C,
        paris_exp_m,
        plot,
        material,
        min_launches)

    while count < min_launches:
        # print("Failure predicted due to Paris equation")
        thickness_tank = thickness_tank + 0.001

        count = fatigue_paris_estimation(
            a_crack_depth,
            thickness_tank,
            count,
            sigma_global_loading,
            stress_range,
            stress_cycle,
            number_of_launches,
            paris_coeff_C,
            paris_exp_m,
            plot,
            material,
            min_launches)

    #applying miners rule 
    damage = miners(stress_range, miner_c, miner_m, stress_cycle, safety_factor, min_launches)[2]
    if damage >1:
        pass
        # print("Failure was predicted due to Miners Law. Thickness will be increased")
    while damage > 1:

        thickness_tank = thickness_tank + 0.0001
        sigma_global_loading = loading_phases(delta_T_earth,
                                              delta_T_tank,
                                              delta_T_tank_extreme,
                                              material,
                                              phi,
                                              thickness_tank,
                                              tank_radius,
                                              pressure_launch_tank,
                                              pressure_after_dock,
                                              launch_mass,
                                              cs,
                                              g_launch_force_ratio,
                                              thrust_engines,
                                              pressure_vent,
                                              pressure_dock_vent,
                                              g_reentry_force_ratio,
                                              dry_mass,
                                              time_mission,
                                              plot)
        stress_range, stress_cycle = rainflow_counting(cycle_launch(sigma_global_loading, 1e-6))
        damage = miners(stress_range, miner_c, miner_m, stress_cycle, safety_factor, min_launches)[2]

        if thickness_tank > 0.01:
            raise RuntimeError("failed on Miners equation fatigue calculation")
    damage,launches = fatigue_miner_estimation(stress_range,miner_c,miner_m,stress_cycle,safety_factor,min_launches,plot)
    # print("Final damage number is:",damage,"which means it will survive",launches,"Launches")
    return thickness_tank


def miner_coefficients(stress, cycles):
    log_stress = np.log10(stress)
    log_N = np.log10(cycles)

    # Linear regression
    slope, intercept, _, _, _ = linregress(log_stress, log_N)

    # Extract coefficients
    m = -slope
    C = 10 ** intercept

    return m, C

if __name__ == "__main__":
#--------------------------------------------------------
#     Main Input Parameters
# --------------------------------------------------------

    #Constants
    # launches = 40
    # safety_factor = 2
    # min_launches = 25

    # Geometry Properties
    phi = 6  # degrees, conical head angle, later import from tank sizing file in final sizing.
    tank_radius = 5  # 5 # m, tank radius, later import from tank sizing file in final sizing.
    thickness_tank = 0.004  # m, tank thickness later import from tank sizing file in final sizing.

    # Material Properties

    # Mat = material.Material(youngs_modulus=material.E,density=material.rho,thermal_expansion_coefficient=material.cte)
    material = mat.Material(
        density=7850,
        youngs_modulus=200e9,
        thermal_expansion_coefficient=1.5e-6,  # Changed from thermal_expansion_coefficient
        fracture_toughness=200,
        yield_strength=500e6)

    # # Mass Inputss
    fuel_reentry_LH2 = 3000  # kg, fuel mass during re-entry
    dry_mass = 20000  # kg, dry mass of the rocket, later import from tank sizing file in final sizing.
    launch_mass = 220_000  # kg, total launch mass to be corrected
    payload_mass = 15500

    # Forces and time inputs
    # time_mission = [0, 0.1, 0.3, 18, 21, 23, 24]  # hours, mission time points
    g_reentry_force_ratio = 6
    g_launch_force_ratio = 6
    max_thrust2weight = 4.3
    # force_launch = g_launch_force_ratio * cs.g_0 * launch_mass
    # thrust_engines = max_thrust2weight * (
    #             dry_mass + payload_mass) * cs.g_0  # N, thrust of engines, later import from dictionary in main branch.

    # Crack Growth Parameters
    # da_dn = []  # initialize da/dn as an empty list
    # a_crack = []  # to store crack growth
    # a_crack_depth = 0.001  # initial crack size in m
    # count = 0
    # Damage = 0

    # Pressure Inputs
    # pressure_launch_tank = 3e5
    # pressure_dock_vent = 5e5
    # pressure_vent = 1e6  # Pa, pressure at which tank is vented. Taken from boil-off file
    # pressure_after_dock = 4.5e5
    # Temperature Inputs
    # T_lh2 = 20  #K
    # T_space = 4  #K, temperature in space
    # T_ambient_earth = 300  #K
    # T_gh2 = 150  #K Temperature of gaseous hydrogen
    # T_gh2_ext = 200  #K temperature of gaseous hydrogen at extreme
    # delta_T_earth = T_ambient_earth - T_lh2
    # delta_T_tank = T_gh2 - T_lh2
    # delta_T_tank_extreme = T_gh2_ext - T_lh2

    #Extrapolated parameters
    # # Paris law
    # paris_coeff_C = 5.131e-17 / 1e3  # convert from mm/cycle to m/cycle
    # paris_exp_m = 7.02  # dimensionless

    # Example Miner's law coefficients for stainless steel 304L
    # Cryogenic paper miner coeff --- 408,-0.02
    # MAT200 paper --- 1884, -0.1555
    # const = 1884
    # exp_coeff = -0.1555
    # miner_m = -1 / exp_coeff
    # miner_c = 10 ** (np.log10(const) / -exp_coeff)
    # # print("The C coefficient for miner equation is", miner_c, "and m coefficient is", miner_m)
    # # Conditions
    # plot = True
    # crack_cond = True
    #--------------------------------------------------------

    thickness_minimum_tank = thickness_optimization_fatigue(phi,
                                       tank_radius,
                                       thickness_tank,
                                       material,
                                       fuel_reentry_LH2,
                                       dry_mass,
                                       launch_mass,
                                       payload_mass,
                                       g_reentry_force_ratio,
                                       g_launch_force_ratio,
                                       max_thrust2weight,
                                       number_of_launches=40,
                                       min_launches=25,
                                       safety_factor=2,
                                       time_mission=[0, 0.1, 0.3, 18, 21, 23, 24],
                                       a_crack_depth=0.001,
                                       pressure_launch_tank=3e5,
                                       pressure_dock_vent=5e5,
                                       pressure_vent=1e6,
                                       pressure_after_dock=4.5e5,
                                       T_lh2=20,
                                       T_space=4,
                                       T_ambient_earth=300,
                                       T_gh2=150,
                                       T_gh2_ext=200)
    # print("Min thickness of tank:", thickness_minimum_tank,"m")
