import numpy as np
from pathlib import Path
import sys
import matplotlib.pyplot as plt

current_dir = Path(__file__).parent
sys.path.append(str(current_dir.parent))
from data import constants as cs
from data import material as mat


def radius_tank(m_h2_header, rho_h2, liquid_fraction):
    volume_tank = m_h2_header / rho_h2 / liquid_fraction
    radius = (volume_tank * 3 / 4 / np.pi) ** (1 / 3)
    return radius


def mass_tank(radius, material, thickness):
    surface_area = 4 * np.pi * radius ** 2
    mass = surface_area * thickness * material.rho
    return mass


def heat_load_tank(deltaT_wall, thickness, thermal_conductivity, radius_header):
    """
    Calculates the heat load into a tank wall based on Fourier's law.

    Parameters:
        surface_area (float): Surface area of the tank (m^2)
        thickness (float): Thickness of the tank wall (m)
        thermal_conductivity (float): Thermal conductivity of the tank material (W/m·K)
        incident_heat_load (float): Incident heat load on the outer surface (W/m^2)

    Returns:
        float: Total heat load conducted into the tank (W)
    """
    # Assuming steady-state, 1D conduction, and that the inside is at 0 K for max load
    r1 = radius_header
    r2 = r1 + thickness
    heat_flux = thermal_conductivity * deltaT_wall / (r1 ** 2 * (1 / r1 - 1 / r2))

    return heat_flux


def stress_sphere_tank(pressure, radius, thickness):
    sigma = pressure * radius / (2 * thickness)
    return sigma


def thickness_optimization(material, pressure, radius, thickness, safety_factor):
    sigma_critical = material.ys / safety_factor
    sigma = stress_sphere_tank(pressure, radius, thickness)
    if sigma > sigma_critical:
        print("Thickness must be increased as failure is predicted")
    else:
        print("Thickness is good. Sigma in walls is", sigma / 1e6, "MPa and sigma critical with safety factor of",
              safety_factor, "is", sigma_critical / 1e6, "Pa")
    while sigma > sigma_critical:
        thickness = thickness + 0.0001
        sigma = stress_sphere_tank(pressure, radius, thickness)
    print("Final thickness is:", thickness, "m")
    return thickness


def mass_dimensions_header(m_h2_header, rho_h2, liquid_fraction, m_o2_header, rho_o2, material, pressure,
                           thickness_h2_header, thickness_o2_header, safety_factor):
    radius_h2_header = radius_tank(m_h2_header, rho_h2, liquid_fraction)
    radius_o2_header = radius_tank(m_o2_header, rho_o2, liquid_fraction)

    thickness_h2_header = thickness_optimization(material, pressure, radius_h2_header, thickness_h2_header,
                                                 safety_factor)
    thickness_o2_header = thickness_optimization(material, pressure, radius_o2_header, thickness_o2_header,
                                                 safety_factor)

    radius_h2_header = radius_tank(m_h2_header, rho_h2, liquid_fraction)
    mass_h2_header = mass_tank(radius_h2_header, material, thickness_h2_header)
    radius_o2_header = radius_tank(m_o2_header, rho_o2, liquid_fraction)
    mass_o2_header = mass_tank(radius_o2_header, material, thickness_o2_header)
    return radius_h2_header, radius_o2_header, mass_h2_header, mass_o2_header

def size_header_tank(m_h2_header,
                     m_o2_header,
                     material,
                     rho_h2=71,
                     rho_o2 = 1141,
                     liquid_fraction = 0.9,
                     thickness_h2_header = 0.0024,
                     thickness_o2_header = 0.003,
                     pressure = 1e6,
                     safety_factor = 2,
                     ):

    (radius_h2_header,
     radius_o2_header,
     mass_h2_header,
     mass_o2_header) = mass_dimensions_header(
        m_h2_header,
        rho_h2,
        liquid_fraction,
        m_o2_header,
        rho_o2,
        material,
        pressure,
        thickness_h2_header,
        thickness_o2_header,
        safety_factor
    )

    return mass_o2_header + mass_h2_header



if __name__ == '__main__':
    m_propellant = 3000
    OF_Ratio = 6
    m_h2_header = m_propellant / OF_Ratio
    m_o2_header = m_propellant - m_h2_header
    rho_h2 = 71  # kg/m3
    rho_o2 = 1141  # kg/m3
    liquid_fraction = 0.9
    thickness_h2_header = 0.0024
    thickness_o2_header = 0.003

    pressure = 1e6
    safety_factor = 2
    thermal_conductivity_insulation = 3.8e-3
    material = mat.Material(density=7850,
                            youngs_modulus=200e9,
                            fracture_strength=800e6,
                            yield_strength=500e6,
                            thermal_conductivity=0.1,
                            specific_heat=500)

    radius_h2_header, radius_o2_header, mass_h2_header, mass_o2_header = mass_dimensions_header(m_h2_header, rho_h2,
                                                                                                liquid_fraction,
                                                                                                m_o2_header, rho_o2,
                                                                                                material, pressure,
                                                                                                thickness_h2_header,
                                                                                                thickness_o2_header,
                                                                                                safety_factor)
    print("Radius of LH2 tank", radius_h2_header, "m with mass", mass_h2_header, "kg.")
    print("Radius of LOX tank", radius_o2_header, "m with mass", mass_o2_header, "kg.")

