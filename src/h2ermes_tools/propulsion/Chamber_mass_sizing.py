import scipy as sp
import numpy as np
import matplotlib.pyplot as plt

#using Inconel 718


def chamber_hoop_stress(pressure, r_c, thickness):
    """
    Calculate the hoop stress in the combustion chamber wall.

    Args:
        pressure (float): Pressure in Pascals.
        r_c (float): Inner radius of the combustion chamber in meters.
        thickness (float): Thickness of the combustion chamber wall in meters.

    Returns:
        float: Hoop stress in Pascals.
    """
    return pressure * r_c / thickness

#ultimate tensile strength of Inconel 718 is the allowable stress, with safety factor of 2.

def obtain_chamber_thickness(pressure, r_c, safety_factor, UTS_strength):
    """
    Calculate the required thickness of the combustion chamber wall.

    Args:
        pressure (float): Pressure in Pascals.
        r_c (float): Inner radius of the combustion chamber in meters.
        safety_factor (float, optional): Safety factor to apply. Defaults to 2.
        UTS_strength (float): Ultimate tensile strength of the material in Pascals.

    Returns:
        float: Required thickness of the combustion chamber wall in meters.
    """
    p_b = safety_factor * pressure #factored burst pressure
    t_w = p_b * r_c / UTS_strength  #thickness of the combustion chamber wall

    return t_w*2 #accounting for stress concentrations and other stresses not modeled in hopp-stress analysis (SPAD)


#linear taper assumption conservative (spad)

def combustion_chamber_mass(r_throat, r_chamber, t_wall, theta_cc, l_cc, rho_material):
    """
    Calculate the mass of the combustion chamber.

    Args:
        r_throat (float): Throat radius in meters.
        r_chamber (float): Chamber radius in meters.
        t_wall (float): Wall thickness in meters.
        theta_cc (float): Half angle of the combustion chamber in radians.
        l_cc (float): Length of the combustion chamber in meters.
        rho_material (float): Density of the material in kg/m^3. 

    Returns:
        float: Mass of the combustion chamber in kg.
    """
    m_cc = np.pi*rho_material*t_wall*(2*r_chamber*l_cc+((r_chamber**2 - r_throat**2)/(np.tan(theta_cc))))
    
    return m_cc

def nozzle_mass(rho_material, t_wall, l_nozzle, theta_cn, r_exit, r_throat):
    """
    Calculate the mass of the nozzle.

    Args:
        rho_material (float): Density of the material in kg/m^3.
        t_wall (float): Wall thickness in meters.
        l_nozzle (float): Length of the nozzle in meters.
        theta_cn (float): Half angle of the nozzle in radians.
        r_exit (float): Exit radius in meters.
        r_throat (float): Throat radius in meters.

    Returns:
        float: Mass of the nozzle in kg.
    """
    m_n = np.pi * rho_material* t_wall* l_nozzle*(r_exit+r_throat)
    
    return m_n

def engine_mass(r_throat, r_chamber, t_wall, theta_cc, l_cc, rho_material, l_nozzle, theta_cn, r_exit):
    """
    Calculate the total mass of the engine, including the combustion chamber and nozzle.

    Args:
        r_throat (float): Throat radius in meters.
        r_chamber (float): Chamber radius in meters.
        t_wall (float): Wall thickness in meters.
        theta_cc (float): Half angle of the combustion chamber in radians.
        l_cc (float): Length of the combustion chamber in meters.
        rho_material (float): Density of the material in kg/m^3.
        l_nozzle (float): Length of the nozzle in meters.
        theta_cn (float): Half angle of the nozzle in radians.
        r_exit (float): Exit radius in meters.

    Returns:
        float: Total mass of the engine in kg.
    """
    m_cc = combustion_chamber_mass(r_throat, r_chamber, t_wall, theta_cc, l_cc, rho_material)
    m_n = nozzle_mass(rho_material, t_wall, l_nozzle, theta_cn, r_exit, r_throat)
    
    return m_cc + m_n

if __name__ == "__main__":
    #for vaccuum
    pressure = 6100000  # Pa
    strength = 1151424468 # Pa (Inconel 718 UTS)
    rho_material = 8220.9316989 # kg/m^3 (Inconel 718 density)
    r_throat = (129.71*1e-3) /2  
    r_chamber = (241.97 * 1e-3) /2
    t_wall = obtain_chamber_thickness(pressure, r_chamber, UTS_strength=strength, safety_factor=2)  # m
    theta_cc = np.deg2rad(30)  # Half angle of the combustion chamber in radians
    l_cc = 185.42*1e-3 #lcyl
    l_nozzle = 1448.83 * 1e-3 #le
    theta_cn = np.deg2rad(38.41)  # Tn Half angle of the nozzle in radians
    r_exit = 1160.18 * 1e-3 / 2  # de/2
    # Calculate the mass of the combustion chamber and nozzle
    m_eng = engine_mass(r_throat, r_chamber, t_wall, theta_cc, l_cc, rho_material, l_nozzle, theta_cn, r_exit)
    print("The mass of the engine for vaccuum is ",m_eng)

    #for sea level
    r_throat = (129.71*1e-3) /2  
    r_chamber = (241.97 * 1e-3) /2
    t_wall = obtain_chamber_thickness(pressure, r_chamber, UTS_strength=strength, safety_factor=2)  # m
    theta_cc = np.deg2rad(30)  # Half angle of the combustion chamber in radians
    l_cc = 185.42*1e-3 #lcyl
    l_nozzle = 334.48 * 1e-3 #le
    theta_cn = np.deg2rad(23.23)  # Tn Half angle of the nozzle in radians
    r_exit = 366.88 * 1e-3 / 2  # de/2
    # Calculate the mass of the combustion chamber and nozzle
    m_eng = engine_mass(r_throat, r_chamber, t_wall, theta_cc, l_cc, rho_material, l_nozzle, theta_cn, r_exit)
    print("The mass of the engine for sea level is ",m_eng)
