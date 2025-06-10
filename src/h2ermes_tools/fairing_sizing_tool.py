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

if __name__ =="__main__":
    radius = 2.5
    thickness = 0.003
    material = mat.Material(density = 7850,
                            youngs_modulus=200e9,
                            fracture_strength=800e6,
                            yield_strength=500e6,
                            thermal_conductivity = 15)
    mass_cap = mass_endcap(material,radius,thickness)
    print("Mass of cap",mass_cap)