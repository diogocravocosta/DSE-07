import numpy as np

from pathlib import Path
import sys

current_dir = Path(__file__).parent
sys.path.append(str(current_dir.parent))

from data import constants as cs
from data import material_prop as mp


#Functions
def thermal_stress(delta_T, young_mod, thermal_expansion_coeff):
    return delta_T * young_mod * thermal_expansion_coeff

def pressure_stress(P,r,t,D):
    hoop_stress = P * r / t
    head_stress = 0.5*P*2*r/t # FOR 2:1 ELLIPTICAL HEADS
    return hoop_stress, head_stress

def crack_growth(paris_coeff_C, paris_exp_m, sigma_max, sigma_min, Y, a):
    delta_K = Y * (sigma_max -sigma_min) * (np.pi*a)
    da_dn = paris_coeff_C * (delta_K ** paris_exp_m)
    return da_dn  

sigma_thermal = thermal_stress(273, mp.ss_E, mp.ss_alpha_expansion)
sigma_launch = 100*10**6
sigma_pressure_wall, sigma_pressure_head = pressure_stress(10**7, 3.5, 0.01, 1)