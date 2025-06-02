import numpy as np

from pathlib import Path
import sys

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
    return hoop_stress, head_stress

def crack_growth(paris_coeff_C, paris_exp_m, sigma_max, sigma_min, Y, a):
    delta_K = Y * (sigma_max -sigma_min) * (np.pi*a)
    da_dn = paris_coeff_C * (delta_K ** paris_exp_m)
    return da_dn  

        #axial stress due to the launch acceleration
def mechanical_stress(thrust_engines, radius, t, fuel_mass, tank_length):
    # Calculate axial stress due to the launch acceleration
    sigma_axial = thrust_engines / (2 * np.pi * (radius) * t)
    # Calculate bending stress due to lateral loads
    sigma_bend = (3 * 9.81 * fuel_mass * tank_length / 2) / (2 * np.pi * (radius)**2 * t)
    return sigma_axial, sigma_bend

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
t = 0.001 # m, tank thickness later import from tank sizing file in final sizing.
tank_length = 20 # m, later import from tank sizing file in final sizing.
thrust_engines = 2.13*10**6 # N, thrust of engines, later import from dictionary in main branch.
force_axial_reentry = thrust_engines # N, axial force during re-entry, for now taken as thrust of engines but will be taken from aero people.

#Extrapolated parameters to apply paris law - 
R = 0.1
paris_coeff_C = 10**(-12)*5.131 # m/cycle
paris_exp_m = 7.02 # dimensionless


while True:

    sigma_thermal = thermal_stress(delta_T, mp.ss_E, mp.ss_alpha_expansion)
    sigma_pressure_wall, sigma_pressure_head = pressure_stress(p_vent, tank_radius, t)
    sigma_axial, sigma_bend = mechanical_stress(thrust_engines, tank_radius, t, fuel_mass, tank_length)
    print("Thermal stress: " + str(sigma_thermal/10**6) + " MPa")
    print("Pressure stress wall: " + str(sigma_pressure_wall/10**6) + " MPa")
    print("Pressure stress head: " + str(sigma_pressure_head/10**6) + " MPa")
    print("Axial stress: " + str(sigma_axial/10**6) + " MPa")
    print("Bending stress: " + str(sigma_bend/10**6) + " MPa")
    #------------------------------------------------------------
    # Different phases of mission and their loading conditions
    #------------------------------------------------------------
    # During launch
    p_launch_tank = 10**5
    sigma_thermal = sigma_thermal # is valid for this.
    sigma_pressure_wall, sigma_pressure_head = pressure_stress(p_launch_tank, tank_radius, t)
    sigma_axial, sigma_bend = mechanical_stress(thrust_engines, tank_radius, t, fuel_mass, tank_length)

    # During orbit
    sigma_thermal = sigma_thermal # will always be taken as its worst case.
    sigma_pressure_wall, sigma_pressure_head = pressure_stress(p_vent, tank_radius, t)

    # During docking
    p_dock_vent =3*10**5
    sigma_pressure_wall, sigma_pressure_head = pressure_stress(p_dock_vent, tank_radius, t)
    
    #During re-entry
    sigma_axial, sigma_bend = mechanical_stress(force_axial_reentry, tank_radius, t, fuel_mass, tank_length)
    break