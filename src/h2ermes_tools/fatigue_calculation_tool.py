import numpy as np
import matplotlib.pyplot as plt

from data import constants as cs
from data import material as mat

#--------------------------------------------------------
#Functions for Stress and Crack Growth Calculations
#--------------------------------------------------------

def mechanical_stress(force, radius, thickness, phi):
    """
    Calculates the mechanical stress in the wall of the conical tank for an axial load.
    Takes into account force, radius, thickness and degree of curvature.
    """
    sigma_conical = force/(2 * np.pi * radius * thickness * np.cos(np.radians(phi)))
    return sigma_conical

def t_crit_mechanical(force, radius, sigma_fracture, phi):
    """Calculates the minimum conical tank wall thickness to not fail for a specified axial load."""
    t = force /(2 * np.pi * radius * sigma_fracture* np.cos(np.radians(phi)))
    return t

def thermal_stress(delta_T, young_mod, thermal_expansion_coeff, phi):
    """
    Calculates the thermal stress in the tank by assuming that the tank is clamped at both ends
    and then uses mechanical stress formula to see effect of thermal induced stress.
    """
    stress = delta_T * young_mod * thermal_expansion_coeff
    thermal_stress = stress/ np.cos(np.radians(phi))
    return thermal_stress

def pressure_stress(Pressure_vapor,radius_tank,thickness, phi):
    """
    Calculates the stress caused by the pressure inside a conical tank. 
    """
    conical_stress = Pressure_vapor /0.85*(radius_tank /thickness/np.cos(np.radians(phi))/0.85+1/6894.76) # Conical stress formula
    return conical_stress

def crack_growth(paris_coeff_C, paris_exp_m, sigma_max, sigma_min, Y_geometry_factor, a_crack_depth):
    """
    Uses the paris equation model crack growth for a sigma max and sigma min for a specified crack depth a. 
    """
    delta_K = Y_geometry_factor * (sigma_max - sigma_min)/10**6 * (np.pi*a_crack_depth)
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

def t_critical(yield_stress,pressure_vapor,phi,radius_tank):
    '''
    Calculates minimum conical tank wall thickness to not fail for a specified pressure load.
    '''
    thickness_crit = pressure_vapor*radius_tank/(np.cos(np.radians(phi))*(yield_stress*0.85 -pressure_vapor/6894.76)) 
    return thickness_crit

def y_calculation(a_crack_depth,thickness):
    """
    Calculates the Y factor that goes into paris equation. Essentially is a geometry factor.
    """
    return 1+1.5*(a_crack_depth/thickness)**2

def miners(stress_range, miner_coeff_C, miner_coeff_m, stress_cycle,safety_factor,launches):
    """
    Calculates if the material is safe from fracture using the Miner rule for fracture. Takes into account stress range and cycle calculated with rainflow.
    """
    sigma = cycle_launch(stress_range,safety_factor)
    cycle = cycle_launch(stress_cycle,launches)
    Di = []
    for i in range(len(sigma)):
        Di.append(cycle[i]*sigma[i]**miner_coeff_m/miner_coeff_C)
    damage = sum(Di)
    return damage

def cycle_launch(stress_cycle, launches):
    """
    Multiplies a given vector (stress_cycle) by a constant (launches)
    """
    cycle = np.zeros(len(stress_cycle))
    for i in range(len(stress_cycle)):
        cycle[i] = stress_cycle[i]*launches
    return cycle

def fatigue_paris_estimation(a_crack_depth, thickness, count, sigma_global_loading, paris_coeff_C, paris_exp_m, pressure_vent, phi, tank_radius, thrust_engines, delta_T_reentry, material,plot):
    """
    Calculates the fatigue life due to paris equation. Also takes into account the thickenss at which failure happens due to 3 separate loading conditions.
    """
    thickness_crit_pressure = t_critical(material.ys, pressure_vent, phi, tank_radius)
    thickness_crit_mech = t_crit_mechanical(thrust_engines,tank_radius, material.fs,phi)
    thickness_crit = max(thickness_crit_pressure,thickness_crit_mech)

    if thickness-thickness_crit<0:
        t = thickness_crit + 0.002

    while a<t-thickness_crit:
        Y_geometry_factor = y_calculation(a_crack_depth, thickness)  # Calculate Y factor based on current crack size
        da_dn = crack_growth(paris_coeff_C, paris_exp_m, max(sigma_global_loading), min(sigma_global_loading), Y_geometry_factor, a_crack_depth) 
        count += 1
        a += da_dn 
        a_crack.append(a_crack_depth) 
        if count > 100000:
            print('no point to worry')
            break
    print('Expected lifecycles: ',count)
    if plot == True:
        cycle = np.arange(0, count, 1)
        plt.figure(figsize=(8, 5))
        plt.plot(cycle, a_crack, marker='o')
        plt.xlabel('Cycle')
        plt.ylabel('Crack Size (m)')
        plt.title('Crack Growth vs Cycle')
        plt.grid(True)
        plt.tight_layout()
        plt.show()
    return a_crack,count

def fatigue_miner_estimation(stress_range, miner_c, miner_m, stress_cycle, safety_factor,launches):
    """
    Calculates the fatigue life due to Miners cycle. Takes into account a safety factor and a conservative number of launches.
    """
    damage = miners(stress_range, miner_c, miner_m, stress_cycle, safety_factor,launches)
    while damage>1:
        launches = launches - 1
        damage = miners(stress_range, miner_c, miner_m, stress_cycle, safety_factor,launches)
    print('Damage count is ', damage, launches)
    return damage,launches

def fatigue_check(a_crack_depth, thickness, count, paris_coeff_C, paris_exp_m, pressure_vent, phi, tank_radius, force_launch, delta_T_reentry,material,plot,min_launch,stress_range, miner_c, miner_m, stress_cycle, safety_fac,launches):
    """
    The main code is here. This is where both paris and miners cycle estimation is done. Returns true if fatigue is not expected. 
    """
    # with paris equation - crack growth
    sigma_global_loading = loading_phases()
    a_crack,count = fatigue_paris_estimation(a_crack_depth, thickness, count, sigma_global_loading, paris_coeff_C, paris_exp_m, pressure_vent, phi, tank_radius, thrust_engines, delta_T_reentry, material,plot)

    #applying miners rule 
    damage, launches = fatigue_miner_estimation(stress_range, miner_c, miner_m, stress_cycle, safety_fac,launches)
    if damage<1 and count>min_launch:
        return True
    else:
        return False

def loading_phases(): 
    """
    Calculates the loading conditions for different phases in the mission. Also calculates the individual load ratios for the different conditions (thermal, pressure, mechanical).
    """
    #------------------------------------------------------------
    # Different phases of mission and their loading conditions
    #------------------------------------------------------------
    # Before launch
    sigma_thermal_start = thermal_stress(delta_T_earth, material.E, material.cte,phi) # is valid for this.
    sigma_pressure_start = pressure_stress(pressure_launch_tank, tank_radius, thickness_tank, phi)
    sigma_axial_start = mechanical_stress(launch_mass*cs.g_0, tank_radius, thickness_tank, phi)
    stress_comb_start = sigma_thermal_start + sigma_pressure_start + sigma_axial_start 

    # 1st stage max q
    sigma_thermal_maxq = 0 #thermal_stress(delta_T_earth, material.E, material.cte,tank_radius, t, phi) # is valid for this.
    sigma_pressure_maxq = pressure_stress(pressure_launch_tank, tank_radius, thickness_tank, phi)
    sigma_axial_maxq = mechanical_stress(launch_mass*g_launch_force_ratio*cs.g_0, tank_radius, thickness_tank, phi)
    stress_comb_maxq = sigma_thermal_maxq + sigma_pressure_maxq + sigma_axial_maxq 

    # During second stage fire
    sigma_thermal_launch = 0 #thermal_stress(delta_T_earth, material.E, material.cte,tank_radius, t, phi)
    sigma_pressure_launch = pressure_stress(pressure_launch_tank, tank_radius, thickness_tank, phi)
    sigma_axial_launch = mechanical_stress(thrust_engines, tank_radius,thickness_tank, phi)
    stress_comb_launch = sigma_thermal_launch + sigma_pressure_launch + sigma_axial_launch 

    # Orbit - Right before docking
    sigma_thermal_orbit = thermal_stress(delta_T_tank, material.E, material.cte, phi)
    sigma_pressure_orbit = pressure_stress(pressure_vent, tank_radius, thickness_tank, phi)
    stress_comb_orbit = sigma_thermal_orbit + sigma_pressure_orbit

    # Orbit - After docking
    sigma_thermal_dock = thermal_stress(delta_T_tank_extreme, material.E, material.cte,phi)
    sigma_pressure_dock= pressure_stress(pressure_dock_vent, tank_radius, thickness_tank, phi)
    stress_comb_dock = sigma_thermal_dock + sigma_pressure_dock

    # Worst point during re-entry
    sigma_thermal_reentry = thermal_stress(delta_T_tank,material.E,material.cte, phi)
    sigma_axial_reentry = mechanical_stress(g_reentry_force_ratio *cs.g_0* dry_mass, tank_radius, thickness_tank, phi)
    sigma_pressure_reentry = pressure_stress(pressure_vent, tank_radius, thickness_tank, phi)
    stress_comb_reentry = sigma_thermal_reentry + sigma_pressure_reentry + sigma_axial_reentry

    # On launch pad
    sigma_thermal_launchpad = thermal_stress(delta_T_earth,material.E,material.cte, phi)
    sigma_axial_launchpad= mechanical_stress(dry_mass*cs.g_0, tank_radius, thickness_tank, phi)
    stress_comb_launchpad = sigma_thermal_launchpad + sigma_axial_launchpad
        
    #------------------------------------------------------------
    # Calculate load ratio's R
    #------------------------------------------------------------

    sigma_global_loading = [stress_comb_start, stress_comb_maxq, stress_comb_launch, stress_comb_orbit, stress_comb_dock, stress_comb_reentry, stress_comb_launchpad]
    R_load_ratio_global = R_calculation(sigma_global_loading)
    print("Global Load ratio R = ", R_load_ratio_global)

    #R for thermal stress 
    sigma_thermal_load = [sigma_thermal_start, sigma_thermal_maxq, sigma_thermal_launch, sigma_thermal_orbit, sigma_thermal_dock, sigma_thermal_reentry, sigma_thermal_launchpad]
    R_load_ratio_thermal = R_calculation(sigma_thermal_load)
    print("Maximum stress ratio R thermal: ", R_load_ratio_thermal)

    #R for mechanical stress
    sigma_mechanical_load = [sigma_axial_start, sigma_axial_maxq, sigma_axial_launch, 0, 0, sigma_axial_reentry, sigma_axial_launchpad]
    R_load_ratio_mechanical = R_calculation(sigma_mechanical_load)
    print("Maximum stress ratio R mechanical: ", R_load_ratio_mechanical)

    #R for pressure stress
    sigma_pressure_load = [sigma_pressure_start, sigma_pressure_maxq, sigma_pressure_launch, sigma_pressure_orbit, sigma_pressure_dock, sigma_pressure_reentry, 0]
    R_load_ratio_pressure = R_calculation(sigma_pressure_load)
    print("Maximum stress ratio R pressure: ", R_load_ratio_pressure)
    if plot == True:
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

if __name__ =="__main__":
    #--------------------------------------------------------
    #Main Input Parameters
    #--------------------------------------------------------

    #Constants 
    launches = 40
    safety_factor = 1.1
    min_launch = 25

    # Geometry Properties
    phi = 20 # degrees, conical head angle, later import from tank sizing file in final sizing.
    tank_radius = 5# 5 # m, tank radius, later import from tank sizing file in final sizing.
    thickness_tank = 0.008 # m, tank thickness later import from tank sizing file in final sizing.
    
    # Material Properties

    # # mat = material.Material(youngs_modulus=material.E,density=material.rho,thermal_expansion_coeffient=material.cte)
    material = mat.Material(
        density = 7850,
        youngs_modulus=200e9,
        fracture_strength=800e6,
        yield_strength=500e6,
        thermal_expansion_coeffient=1.5e-6
    )

    # Mass Inputs
    fuel_reentry_LH2 = 3000 #kg, fuel mass during re-entry
    dry_mass = 30000 # kg, dry mass of the rocket, later import from tank sizing file in final sizing.
    launch_mass = 200000# kg, total launch mass to be corrected
    payload_mass = 15500

    # Forces and time inputs
    time_mission = [0, 0.1, 0.3, 18, 21, 23, 24] # hours, mission time points
    g_reentry_force_ratio = 8 
    g_launch_force_ratio = 6 
    max_thrust2weight = 7.9
    force_launch = g_launch_force_ratio*cs.g_0*launch_mass 
    thrust_engines = max_thrust2weight * (dry_mass + payload_mass) * cs.g_0# N, thrust of engines, later import from dictionary in main branch.

    # Crack Growth Parameters
    da_dn = [] # initialize da/dn as an empty list
    a_crack = [] # to store crack growth
    a_crack_depth = 0.001 # initial crack size in m
    count = 0
    Damage= 0

    # Pressure Inputs
    pressure_launch_tank = 1e5
    pressure_dock_vent =3e5
    pressure_vent = 1e6 # Pa, pressure at which tank is vented. Taken from boil-off file

    # Temperature Inputs
    T_lh2 = 20 #K
    T_space=4 #K, temperature in space
    T_ambient_earth = 300 #K
    T_gh2 = 150 #K Temperature of gaseous hydrogen
    T_gh2_ext = 200 #K temperature of gaseous hydrogen at etreme
    delta_T_earth = T_ambient_earth - T_lh2 
    delta_T_tank = T_gh2 - T_lh2
    delta_T_tank_extreme = T_gh2_ext - T_lh2

    #Extrapolated parameters 
    # Paris law
    paris_coeff_C = 5.131e-11 # m/cycle
    paris_exp_m = 7.02 # dimensionless
    # Miners Law
    miner_c = 1.8e12 #mpa
    miner_m = 3.2
    # Rainflow
    stress_range = [80,300,360,430]#mpa
    stress_cycle = [1,1,0.5,0.5]

    # Conditions
    plot = False
    crack_cond = False
    #--------------------------------------------------------

    fatigue = fatigue_check(a_crack_depth,
                            thickness_tank,
                            count,
                            paris_coeff_C,
                            paris_exp_m,
                            pressure_vent,
                            phi,
                            tank_radius,
                            force_launch,
                            delta_T_reentry,
                            material,
                            plot,
                            min_launch,
                            stress_range,
                            miner_c,
                            miner_m,
                            stress_cycle,
                            safety_factor,
                            launches)
    
    if fatigue == True:
        print("Fatigue failure not expected due to crack propagation and fracture failure. Update stress range values for more precision.")
    else:
        print("!Warning! = Fatigue failure expected. Please increase thickness")