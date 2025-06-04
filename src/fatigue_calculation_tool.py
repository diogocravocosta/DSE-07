import numpy as np
from pathlib import Path
import sys
import matplotlib.pyplot as plt
current_dir = Path(__file__).parent
sys.path.append(str(current_dir.parent))
from data import constants as cs
from scipy.interpolate import interp1d
import rainflow

#--------------------------------------------------------
#Main Input Parameters
#--------------------------------------------------------
g_0 = 9.81
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
p_vent = 10**6 # Pa, pressure at which tank is vented. Taken from boil-off file
fuel_mass = 42000 # kg, later import from dictionary in main branch
t = 0.007 # m, tank thickness later import from tank sizing file in final sizing.
tank_length = 14 # m, later import from tank sizing file in final sizing.
thrust_engines = 2.13*10**6 # N, thrust of engines, later import from dictionary in main branch.
fuel_reentry_LH2 = 3000 #kg, fuel mass during re-entry
dry_mass = 30000 # kg, dry mass of the tank, later import from tank sizing file in final sizing.
launch_mass = 200000# kg, total launch mass
thrust_engines = 7.9 * (dry_mass+15000)*g_0
da_dn = []
a_crack = [] # to store crack growth
a = 0.002 # initial crack size in m
count = 0
ss_E = 200e9  # Pa, steel
ss_rho = 7850 # kg/m3, steel
ss_alpha_expansion = 1.5e-6 #1/K
ss_sigma_f = 600e6 # Pa, fatigue strength of steel
ss_fatigue_m = 3
g_reentry = 8
g_launch = 6
p_launch_tank = 10**5
p_dock_vent =3*10**5
phi = 20 # degrees, conical head angle, later import from tank sizing file in final sizing.
r_small = 2.492 # m, radius of the small end of the tank, later import from tank sizing file in final sizing.
tank_radius =r_small# 5 # m, tank radius, later import from tank sizing file in final sizing.
plot = True
crack_cond = True
#--------------------------------------------------------
#Extrapolated parameters to apply paris law - 
#--------------------------------------------------------
paris_coeff_C = 5.131*10**(-14) # m/cycle
paris_exp_m = 7.02 # dimensionless

#--------------------------------------------------------
#Functions for Stress and Crack Growth Calculations
#--------------------------------------------------------
def mechanical_stress(force, radius, t, phi):
    sigma_conical = force/(2 * np.pi * radius * t * np.cos(np.radians(phi)))
    return sigma_conical

def thermal_stress(delta_T, young_mod, thermal_expansion_coeff,r,t, phi):
    pressure = delta_T * young_mod * thermal_expansion_coeff
    force = pressure * 2*np.pi * r * t
    thermal_stress = mechanical_stress(force, r, t, phi)
    return thermal_stress

def pressure_stress(P,r_big,t, phi):
    conical_stress = P *(r_big /t/np.cos(np.radians(phi))/0.85+1/0.85) # Conical stress formula
    return conical_stress

def crack_growth(paris_coeff_C, paris_exp_m, sigma_max, sigma_min, Y, a):
    delta_K = Y * (sigma_max - sigma_min)/10**6 * (np.pi*a)
    da_dn = paris_coeff_C * (delta_K ** paris_exp_m)
    return da_dn  

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
sigma_max = max(sigma_global_loading)
sigma_min = min(sigma_global_loading)
R_global = sigma_min/sigma_max
print("Global Load ratio R = ", R_global)

#R for thermal stress 
sigma_thermal_load = [sigma_thermal_start, sigma_thermal_maxq, sigma_thermal_launch, sigma_thermal_orbit, sigma_thermal_dock, sigma_thermal_reentry, sigma_thermal_launchpad]
t_mission = [0,0.2,0.5,18,21,23.5,24]
sigma_thermal_max = max(sigma_thermal_load)
s = []
for i in range(len(sigma_thermal_load)):
    if sigma_thermal_load[i] != 0:
        s.append(sigma_thermal_load[i]) 
sigma_thermal_min = min(s)
R_thermal = sigma_thermal_min / sigma_thermal_max
print("Maximum stress ratio R thermal: ", R_thermal)

#R for mechanical stress
sigma_mechanical_load = [sigma_axial_start, sigma_axial_maxq, sigma_axial_launch, 0, 0, sigma_axial_reentry, sigma_axial_launchpad]
t_mission =[0,0.2,0.5,18,21,23.5,24]
s = []
for i in range(len(sigma_mechanical_load)):
    if sigma_mechanical_load[i] != 0:
        s.append(sigma_mechanical_load[i]) 
sigma_mechanical_min = min(s)
sigma_mechanical_max = max(sigma_mechanical_load)
R_mechanical = sigma_mechanical_min / sigma_mechanical_max
print("Maximum stress ratio R mechanical: ", R_mechanical)

#R for pressure stress
sigma_pressure_load = [sigma_pressure_start, sigma_pressure_maxq, sigma_pressure_launch, sigma_pressure_orbit, sigma_pressure_dock, sigma_pressure_reentry, 0]
t_mission =[0,0.2,0.5,18,21,23.5,24]
s= []
for i in range(len(sigma_thermal_load)):
    if sigma_pressure_load[i] != 0:
        s.append(sigma_pressure_load[i]) 
sigma_pressure_min = min(s)
sigma_pressure_max = max(sigma_pressure_load)
R_pressure = sigma_pressure_min / sigma_pressure_max
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

while crack_cond == True:
    da_dn = crack_growth(paris_coeff_C, paris_exp_m, sigma_thermal_max, sigma_thermal_min, 1.12, a) 
    count += 1
    a += da_dn 
    a_crack.append(a) 
    if count>15000:
        print('Fatigue failure not expected')
        break
    if a > t:  # Stop when crack size exceeds 0.1 m
        cycle = np.arange(0, count, 1)
        plt.figure(figsize=(8, 5))
        plt.plot(cycle, a_crack, marker='o')
        plt.xlabel('Cycle')
        plt.ylabel('Crack Size (m)')
        plt.title('Crack Growth vs Cycle')
        plt.grid(True)
        plt.tight_layout()
        plt.show()
        break
