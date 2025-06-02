import numpy as np

from pathlib import Path
import sys
import matplotlib.pyplot as plt

current_dir = Path(__file__).parent
sys.path.append(str(current_dir.parent))

from data import constants as cs
from data import material_prop as mp
#from final_sizing import tank_sizing as ts

#Functions
def thermal_stress(delta_T, young_mod, thermal_expansion_coeff):
    return delta_T * young_mod * thermal_expansion_coeff

def pressure_stress(P,r,t):
    hoop_stress = P * r / t
    head_stress = 0.5*P*(2*r)/t # FOR 2:1 ELLIPTICAL HEADS
    return max(hoop_stress,head_stress)

def mechanical_stress(thrust_engines, radius, t, fuel_mass, tank_length):
    # Calculate axial stress due to the launch acceleration
    sigma_axial = thrust_engines / (2 * np.pi * (radius) * t)
    # Calculate bending stress due to lateral loads
    sigma_bend = (3 * 9.81 * fuel_mass * tank_length / 2) / (2 * np.pi * (radius)**2 * t)
    return sigma_axial, sigma_bend

def crack_growth(paris_coeff_C, paris_exp_m, sigma_max, sigma_min, Y, a):
    delta_K = Y * (sigma_max - sigma_min)/10**6 * (np.pi*a)
    da_dn = paris_coeff_C * (delta_K ** paris_exp_m)
    return da_dn  

#sigma_cr = gamma * 0.605 * E * t / radius


#--------------------------------------------------------
#Main Input Parameters
#--------------------------------------------------------
T_lh2 = 20 #K
T_ambientn_earth = 300 #K
delta_T = T_ambientn_earth - T_lh2
p_vent = 10**6 # Pa, pressure at which tank is vented. Taken from boil-off file
tank_radius = 3.5 # m, tank radius, later import from tank sizing file in final sizing.
fuel_mass = 42000 # kg, later import from dictionary in main branch
t = 0.01 # m, tank thickness later import from tank sizing file in final sizing.
tank_length = 20 # m, later import from tank sizing file in final sizing.
thrust_engines = 2.13*10**6 # N, thrust of engines, later import from dictionary in main branch.
force_axial_reentry = thrust_engines # N, axial force during re-entry, for now taken as thrust of engines but will be taken from aero people.
fuel_reentry_LH2 = 3000 #kg, fuel mass during re-entry
dry_mass = 40000 # kg, dry mass of the tank, later import from tank sizing file in final sizing.
da_dn = []
a_crack = [] # to store crack growth
a = 0.001 # initial crack size in m
count = 0

#Extrapolated parameters to apply paris law - 
paris_coeff_C = 3.6e-8#5.131*10**(-14) # m/cycle
paris_exp_m = 3.44 # dimensionless


sigma_thermal = thermal_stress(delta_T, mp.ss_E, mp.ss_alpha_expansion)
sigma_pressure = pressure_stress(p_vent, tank_radius, t)
sigma_axial, sigma_bend = mechanical_stress(thrust_engines, tank_radius, t, fuel_mass, tank_length)
print("Thermal stress: " + str(sigma_thermal/10**6) + " MPa")
print("Pressure stress on tank wall: " + str(sigma_pressure/10**6) + " MPa")
print("Axial stress: " + str(sigma_axial/10**6) + " MPa")
print("Bending stress: " + str(sigma_bend/10**6) + " MPa")
#------------------------------------------------------------
# Different phases of mission and their loading conditions
#------------------------------------------------------------
# During launch, orbit, docking, re-entry and on launch pad
# Note: The stresses are calculated for the worst case scenario in each phase.
# During launch
p_launch_tank = 10**5
sigma_thermal = sigma_thermal # is valid for this.
sigma_pressure_launch = pressure_stress(p_launch_tank, tank_radius, t)
sigma_axial, sigma_bend = mechanical_stress(thrust_engines, tank_radius, t, fuel_mass, tank_length)
stress_comb_launch = sigma_thermal + sigma_pressure_launch + sigma_axial + sigma_bend 

# During orbit
sigma_pressure_orbit = pressure_stress(p_vent, tank_radius, t)
stress_comb_orbit = sigma_thermal + sigma_pressure_orbit

# During docking
p_dock_vent =3*10**5
sigma_pressure_dock= pressure_stress(p_dock_vent, tank_radius, t)
stress_comb_dock = sigma_pressure_dock

# During re-entry
sigma_thermal = sigma_thermal # will always be taken as its worst case.
sigma_axial, sigma_bend = mechanical_stress(force_axial_reentry, tank_radius, t, fuel_reentry_LH2, tank_length)
sigma_pressure_reentry = pressure_stress(p_vent, tank_radius, t)
stress_comb_reentry = sigma_thermal + sigma_pressure_reentry + sigma_axial + sigma_bend    

# On launch pad
sigma_thermal = sigma_thermal # will always be taken as its worst case.
sigma_axial, sigma_bend = mechanical_stress(dry_mass*cs.g_0, tank_radius, t, fuel_reentry_LH2, tank_length)
stress_comb_launchpad = sigma_thermal + sigma_axial + sigma_bend

sigma_max = max(stress_comb_launch, stress_comb_orbit, stress_comb_dock, stress_comb_reentry, stress_comb_launchpad)
sigma_min = min(stress_comb_launch, stress_comb_orbit, stress_comb_dock, stress_comb_reentry, stress_comb_launchpad)
R = sigma_min/sigma_max
print("Maximum stress during mission: ", sigma_max)
print("Load ratio R = ", R)
sigma_plot = [stress_comb_launch, stress_comb_orbit, stress_comb_dock, stress_comb_reentry, stress_comb_launchpad]
time_mission = [0.1,1.5,1.7,1.9,2.0]
plt.figure(figsize=(8, 5))
plt.plot(time_mission, sigma_plot, marker='o')
plt.xlabel('Mission Time (day)')
plt.ylabel('Combined Stress (Pa)')
plt.title('Combined Stress vs Mission Time')
plt.grid(True)
plt.tight_layout()
plt.show()

while True:
    da_dn = crack_growth(paris_coeff_C, paris_exp_m, sigma_max, sigma_min, 1.12, a) 
    count += 1
    a += da_dn 
    a_crack.append(a) 

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
