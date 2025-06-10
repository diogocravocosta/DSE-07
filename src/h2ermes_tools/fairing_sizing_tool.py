import numpy as np
from pathlib import Path
import sys
import matplotlib.pyplot as plt
current_dir = Path(__file__).parent
sys.path.append(str(current_dir.parent))
from data import constants as cs
from data import material as mat


#spherical end cap

def volume_spherical_endcap_sheet(radius,thickness):
    surface_area = 2*np.pi*radius**2
    volume_sheet = surface_area*thickness
    return volume_sheet

def mass_endcap(material,radius,thickness):
    mass = material.rho*volume_spherical_endcap_sheet(radius,thickness)
    return mass

def chapman_stagnation_heat_flux(radius_nose,material,velocity_max,air_rho,n,m,c_star,v_c):
    qs = 1.63e-4*(air_rho/radius_nose)**0.5*velocity_max**3
    c1 = c_star * (1 / np.sqrt(air_rho)) * (1 / v_c ** 3)
    qc_max = c1 * (1 / radius_nose^n) * (air_rho)^(1-n) * (velocity_max)^m
    return qs,qc_max

def 

def chapman_stagnation_heat_flux(self):
    c1 = self.c_star * (1 / np.sqrt(self.rho_0)) * (1 / self.v_c ** 3)
    c2 = c1 * (1 / self.nose_radius ** self.n) * (2 * self.lift_parameter / self.v_c ** 2) ** (
                1 - self.n) * self.v_c ** self.constant.m

    rho_max = 4 * (self.lift_parameter / (self.v_c ** 2)) * ((1 - self.n) / (self.constant.m + 2 * (self.n - 1)))
    V_max = np.sqrt((self.constant.m - 2 * (1 - self.n)) / self.constant.m) * self.v_c
    qc_max = c1 * (1 / self.nose_radius ** self.n) * (rho_max) ** (1 - self.n) * (V_max) ** self.constant.m
    qc = c2 * (((1 / self.normalized_v_ratio ** 2) - 1) ** (
                1 - self.n)) * self.normalized_v_ratio ** self.constant.m
    normalised_heat_flux = qc/ qc_max


    return qc*1e-3, qc_max*1e-3, normalised_heat_flux

if __name__ =="__main__":
    radius = 2.5
    thickness = 0.003
    radius_nose = 0.5
    n_coefficient = 0.5
    m_coefficient = 3
    c_star = 1.18e8
    velocity_max = 
    air_rho = 
    V_c
    material = mat.Material(density = 7850,
                            youngs_modulus=200e9,
                            fracture_strength=800e6,
                            yield_strength=500e6,
                            thermal_conductivity = 15)
    mass_cap = mass_endcap(material,radius,thickness)
    print("Mass of cap",mass_cap)

    qs , qc_max= chapman_stagnation_heat_flux(radius_nose,material,velocity_max,air_rho,n_coefficient,m_coefficient,c_star,v_c)
    print(qs, qc_max)