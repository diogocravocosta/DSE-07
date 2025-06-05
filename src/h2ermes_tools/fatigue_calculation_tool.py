import numpy as np
from pathlib import Path
import sys
import matplotlib.pyplot as plt
current_dir = Path(__file__).parent
sys.path.append(str(current_dir.parent))
from data import constants as cs
from scipy.interpolate import interp1d
from scipy.signal import find_peaks
import rainflow


#--------------------------------------------------------
#Main Input Parameters
#--------------------------------------------------------
#Constants 
g_0 = 9.81 # m/s^2, gravitational acceleration
T_lh2 = 20 #K
T_space=4 #K, temperature in space
T_ambient_earth = 300 #K
T_gh2 = 150
T_gh2_ext = 200
delta_T_earth = T_ambient_earth - T_lh2 
delta_T_reentry = 150-T_lh2 #K, temperature difference during re-entry, will be changed later for reentry heat.
delta_T_space = T_lh2 - T_space
delta_T_h2 = T_gh2 - T_lh2
delta_T_ext = T_gh2_ext - T_lh2
launches = 40
safety_fac = 1.4

# Geometry Properties
phi = 20 # degrees, conical head angle, later import from tank sizing file in final sizing.
r_small = 2.492 # m, radius of the small end of the tank, later import from tank sizing file in final sizing.
tank_radius = 5# 5 # m, tank radius, later import from tank sizing file in final sizing.
tank_length = 14 # m, later import from tank sizing file in final sizing.
t = 0.008 # m, tank thickness later import from tank sizing file in final sizing.
 
# Material Properties
ss_E = 200e9  # Pa, steel
ss_rho = 7850 # kg/m3, steel
ss_alpha_expansion = 1.5e-6 #1/K
ss_sigma_f = 1200e6 # Pa, fatigue strength of steel
ss_fatigue_m = 3
ss_yield = 500e6

# Mass Inputs
fuel_reentry_LH2 = 3000 #kg, fuel mass during re-entry
dry_mass = 30000 # kg, dry mass of the tank, later import from tank sizing file in final sizing.
launch_mass = 200000# kg, total launch mass
fuel_mass = 42000 # kg, later import from dictionary in main branch
m_payload = 15500

# Forces and time inputs
t_mission = [0,0.1,0.3,18,21,23,24] # hours, mission time points
thrust_engines = 7.9 * (dry_mass+m_payload)*g_0# N, thrust of engines, later import from dictionary in main branch.
g_reentry = 8
g_launch = 6 

# Crack Growth Parameters
da_dn = [] # initialize da/dn as an empty list
a_crack = [] # to store crack growth
a = 0.001 # initial crack size in m
count = 0
Dcount = 0

# Pressure Inputs
p_launch_tank = 10**5
p_dock_vent =3*10**5
p_vent = 10**6 # Pa, pressure at which tank is vented. Taken from boil-off file

#Extrapolated parameters to apply paris law - 
paris_coeff_C = 5.131*10**(-11) # m/cycle
paris_exp_m = 7.02 # dimensionless

miner_c = 1.8*10e12 #mpa
miner_m = 3.2
stress_range = [80,300,360,430]#mpa
stress_cycle = [1,1,0.5,0.5]

# Conditions
plot = False
crack_cond = False

#--------------------------------------------------------
#Functions for Stress and Crack Growth Calculations
#--------------------------------------------------------
def mechanical_stress(force, radius, t, phi):
    sigma_conical = force/(2 * np.pi * radius * t * np.cos(np.radians(phi)))
    return sigma_conical

def t_crit_mechanical(force,radius,sigma_fracture,phi):
    t = force /(2 * np.pi * radius * sigma_fracture* np.cos(np.radians(phi)))
    return t

def thermal_stress(delta_T, young_mod, thermal_expansion_coeff,r,t, phi):
    pressure = delta_T * young_mod * thermal_expansion_coeff
    force = pressure * 2*np.pi * r * t
    thermal_stress = mechanical_stress(force, r, t, phi)
    return thermal_stress

def t_crit_theraml(delta_T, young_mod, thermal_expansion_coeff,r,t, phi, sigma_fracture):
    pressure = delta_T * young_mod * thermal_expansion_coeff
    force = pressure * 2*np.pi * r * t
    t = t_crit_mechanical(force,r,sigma_fracture, phi)
    return t

def pressure_stress(P,r_big,t, phi):
    conical_stress = P /0.85*(r_big /t/np.cos(np.radians(phi))/0.85+1/6894.76) # Conical stress formula
    return conical_stress

def crack_growth(paris_coeff_C, paris_exp_m, sigma_max, sigma_min, Y, a):
    delta_K = Y * (sigma_max - sigma_min)/10**6 * (np.pi*a)
    da_dn = paris_coeff_C * (delta_K ** paris_exp_m)
    return da_dn  

def R_calculation(loading):
    sigma_max = max(loading)
    s = []
    for i in range(len(loading)):
        if loading[i] != 0:
            s.append(loading[i]) 
    sigma_min = min(s)
    R = sigma_min / sigma_max   
    return R

def t_critical(yield_stress,p,phi,r_outer):
    t_crit = p*r_outer/(np.cos(np.radians(phi))*(yield_stress*0.85 -p/6894.76)) # Critical thickness for crack growth
    return t_crit

def y_calculation(a,t):
    return 1+1.5*(a/t)**2

def miners(stress_range, C, m, stress_cycle,sf,launches):
    sigma = cycle_launch(stress_range,sf)
    cycle = cycle_launch(stress_cycle,launches)
    Di = []
    for i in range(len(sigma)):
        Di.append(cycle[i]*sigma[i]**m/C)
    Dcount = sum(Di)
    return Dcount

def cycle_launch(stress_cycle, launches):
    cycle = np.zeros(len(stress_cycle))
    for i in range(len(stress_cycle)):
        cycle[i] = stress_cycle[i]*launches
    return cycle

def fatigue_estimation(a, t, count, sigma_global_loading, paris_coeff_C, paris_exp_m, ss_sigma_f, p_vent, phi, tank_radius, thrust_engines, delta_T_reentry,ss_yield,plot):
    t_crit_pressure = t_critical(ss_yield, p_vent, phi, tank_radius)
    t_crit_mech = t_crit_mechanical(thrust_engines,tank_radius, ss_sigma_f,phi)
    t_crit_thermal = t_crit_theraml(delta_T_reentry, ss_E, ss_alpha_expansion, tank_radius, t, phi, ss_sigma_f)
    t_crit = max(t_crit_pressure,t_crit_mech, t_crit_thermal)

    if t-t_crit<0:
        t = t_crit + 0.002

    while a<t-t_crit:
        Y = y_calculation(a, t)  # Calculate Y factor based on current crack size
        da_dn = crack_growth(paris_coeff_C, paris_exp_m, max(sigma_global_loading), min(sigma_global_loading), Y, a) 
        count += 1
        a += da_dn 
        a_crack.append(a) 
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

def fatigue_miner_estimation(stress_range, miner_c, miner_m, stress_cycle, sf,launches):
    Dcount = miners(stress_range, miner_c, miner_m, stress_cycle, sf,launches)
    while Dcount>1:
        launches = launches - 1
        Dcount = miners(stress_range, miner_c, miner_m, stress_cycle, sf,launches)
    print('Dcount minimum is ', Dcount, launches)
    return Dcount,launches

#------------------------------------------------------------
# Different phases of mission and their loading conditions
#------------------------------------------------------------

# Before Launch
sigma_thermal_start = thermal_stress(delta_T_earth, ss_E, ss_alpha_expansion,tank_radius, t, phi) # is valid for this.
sigma_pressure_start = pressure_stress(p_launch_tank, tank_radius, t, phi)
sigma_axial_start = mechanical_stress(launch_mass*g_0, tank_radius, t, phi)
stress_comb_start = sigma_thermal_start + sigma_pressure_start + sigma_axial_start 

# 1st stage max q
sigma_thermal_maxq = 0 #thermal_stress(delta_T_earth, ss_E, ss_alpha_expansion,tank_radius, t, phi) # is valid for this.
sigma_pressure_maxq = pressure_stress(p_launch_tank, tank_radius, t, phi)
sigma_axial_maxq = mechanical_stress(launch_mass*g_launch*g_0, tank_radius, t, phi)
stress_comb_maxq = sigma_thermal_maxq + sigma_pressure_maxq + sigma_axial_maxq 

# During second stage fire
sigma_thermal_launch = 0 #thermal_stress(delta_T_earth, ss_E, ss_alpha_expansion,tank_radius, t, phi)
sigma_pressure_launch = pressure_stress(p_launch_tank, tank_radius, t, phi)
sigma_axial_launch = mechanical_stress(thrust_engines, tank_radius, t, phi)
stress_comb_launch = sigma_thermal_launch + sigma_pressure_launch + sigma_axial_launch 

# Orbit - Right before docking
sigma_thermal_orbit = thermal_stress(delta_T_h2, ss_E, ss_alpha_expansion,tank_radius, t, phi)
sigma_pressure_orbit = pressure_stress(p_vent, tank_radius, t, phi)
stress_comb_orbit = sigma_thermal_orbit + sigma_pressure_orbit

# Orbit - After docking
sigma_thermal_dock = thermal_stress(delta_T_ext, ss_E, ss_alpha_expansion,tank_radius, t, phi)
sigma_pressure_dock= pressure_stress(p_dock_vent, tank_radius, t, phi)
stress_comb_dock = sigma_thermal_dock + sigma_pressure_dock

# Worst point during re-entry
sigma_thermal_reentry = thermal_stress(delta_T_h2,ss_E,ss_alpha_expansion,tank_radius, t, phi) 
sigma_axial_reentry = mechanical_stress(g_reentry *g_0* dry_mass, tank_radius, t, phi)
sigma_pressure_reentry = pressure_stress(p_vent, tank_radius, t, phi)
stress_comb_reentry = sigma_thermal_reentry + sigma_pressure_reentry + sigma_axial_reentry

# On launch pad
sigma_thermal_launchpad = thermal_stress(delta_T_earth,ss_E,ss_alpha_expansion,tank_radius, t, phi) 
sigma_axial_launchpad= mechanical_stress(dry_mass*g_0, tank_radius, t, phi)
stress_comb_launchpad = sigma_thermal_launchpad + sigma_axial_launchpad

#------------------------------------------------------------
# Calculate load ratio's R
#------------------------------------------------------------
# Global load ratio R
sigma_global_loading = [stress_comb_start, stress_comb_maxq, stress_comb_launch, stress_comb_orbit, stress_comb_dock, stress_comb_reentry, stress_comb_launchpad]
R_global = R_calculation(sigma_global_loading)
print("Global Load ratio R = ", R_global)

#R for thermal stress 
sigma_thermal_load = [sigma_thermal_start, sigma_thermal_maxq, sigma_thermal_launch, sigma_thermal_orbit, sigma_thermal_dock, sigma_thermal_reentry, sigma_thermal_launchpad]
R_thermal = R_calculation(sigma_thermal_load)
print("Maximum stress ratio R thermal: ", R_thermal)

#R for mechanical stress
sigma_mechanical_load = [sigma_axial_start, sigma_axial_maxq, sigma_axial_launch, 0, 0, sigma_axial_reentry, sigma_axial_launchpad]
R_mechanical = R_calculation(sigma_mechanical_load)
print("Maximum stress ratio R mechanical: ", R_mechanical)

#R for pressure stress
sigma_pressure_load = [sigma_pressure_start, sigma_pressure_maxq, sigma_pressure_launch, sigma_pressure_orbit, sigma_pressure_dock, sigma_pressure_reentry, 0]
R_pressure = R_calculation(sigma_pressure_load)
print("Maximum stress ratio R pressure: ", R_pressure)

#Plots
if plot == True:
    plt.figure(figsize=(8, 5))
    plt.plot(t_mission, np.array(sigma_pressure_load)/1e6, marker='o', label='Pressure Stress')
    plt.plot(t_mission, np.array(sigma_thermal_load)/1e6, marker='s', label='Thermal Stress')
    plt.plot(t_mission, np.array(sigma_mechanical_load)/1e6, marker='^', label='Mechanical Stress')
    plt.plot(t_mission, np.array(sigma_global_loading)/1e6, marker='x', label='Total Stress')
    plt.xlabel('Mission Time (hours)')
    plt.ylabel('Stress (MPa)')
    plt.title('Pressure, Thermal, Mechanical and Total Stress vs Mission Time')
    plt.legend()
    plt.grid(True)
    plt.tight_layout()
    plt.show()

#------------------------------------------------------------
# Crack Growth Calculation
#------------------------------------------------------------

# with paris equation - crack growth
a_crack,count = fatigue_estimation(a, t, count, sigma_global_loading, paris_coeff_C, paris_exp_m, ss_sigma_f, p_vent, phi, tank_radius, g_launch*g_0*launch_mass, delta_T_reentry,ss_yield,plot)

#applying miners rule 
Dcount, launches = fatigue_miner_estimation(stress_range, miner_c, miner_m, stress_cycle, safety_fac,launches)

