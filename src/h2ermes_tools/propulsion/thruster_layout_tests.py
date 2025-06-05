import numpy as np
import matplotlib.pyplot as plt
import math
import 

# Data
g0 = 9.80665
#expansion_ratios = np.array([40, 50, 60, 70, 80, 85, 90, 95, 100, 110, 120, 130, 140])
#vacuum_Isp = np.array([435.5235, 439.7381, 443.0285, 445.669, 447.8579, 448.8227, 449.7163, 450.5474, 451.3234, 452.7236, 453.9782, 455.0922, 456.0993])
#sea_level_Isp = np.array([274.4114, 238.3587, 201.3354, 163.7011, 125.6168, 106.4455, 87.2033, 67.8989, 48.5397, 9.6704,-29.7567, -68.9614, -108.2737])
#sigma = 0.121029372  # Structural mass ratio
#structural_mass = 20642  # Structural mass in kg
#delta_V_vacuum = 7264.29 # Delta V in vacuum in m/s
#delta_V_sea_level = 250  # Delta V in m/s for landing burn
#structural_mass = 21 #assuming structural mass is constant throughout the mission
#payload_mass = 15000
#REVISE THIS VALUE

def mprop_calc(delta_V, Isp, mf):
    """
    Calculate the propellant mass required for a given delta V, specific impulse, and structural mass.

    Args:
        delta_V (float): Delta V in m/s.
        Isp (float): Specific impulse in seconds.
        mf (float): Final mass in kg.

    Returns:
        float: Propellant mass required in kg.
    """
    return mf * (math.exp(delta_V / (9.80665 * Isp)) - 1)


def sl_vac_thruster_perf(n_sl, n_vac, Isp_sl, Isp_vac, t_sl, t_vac, mdot_sl, mdot_vac, m_struct, payload, delta_v_sl, delta_v_vac, delta_v_deorbit):
    """
    Calculate the performance of a thruster in sea level and vacuum conditions.

    Args:
        n_sl (int): Number of thrusters in sea level.
        n_vac (int): Number of thrusters in vacuum.
        Isp_sl (float): Specific impulse in sea level.
        Isp_vac (float): Specific impulse in vacuum.
        t_sl (float): Thrust in sea level.
        t_vac (float): Thrust in vacuum.
        mdot_sl (float): Mass flow rate in sea level.
        mdot_vac (float): Mass flow rate in vacuum.
        m_struct (float): Structural mass in kg.
        payload (float): Payload mass in kg.
        delta_v_sl (float): Delta V in sea level in m/s.
        delta_v_vac (float): Delta V in vacuum in m/s.
        delta_v_deorbit (float): Delta V for deorbit burn in m/s.

    Returns:
        
    """
    
