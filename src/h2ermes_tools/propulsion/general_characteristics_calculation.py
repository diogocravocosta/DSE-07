import numpy as np
import matplotlib.pyplot as plt


def total_thrust_to_individual_chamber_thrust(t_total, Isp_sl = 393.3471, Isp_vac = 447.9481, sea_level=8, total = 24):
    
    """Calculate the individual chamber thrusts based on total thrust and specific impulses.
    Parameters:
    t_total (float): Total thrust in Newtons.
    Isp_sl (array): Specific impulse at sea level for the sea level [0] and vacuum [1] optimized chambers in seconds.
    Isp_vac (array): Specific impulse in vacuum for sea level [0] and vacuum [1] optimized chambers in seconds.
    sea_level (int): Number of sea level chambers.
    total (int): Total number of chambers.
    Returns:
    tuple: Individual thrusts at sea level and in vacuum.
    """
    
    #t_sl_opt = t_total/(sea_level + ((total - sea_level)*(Isp_sl[0]/Isp_sl[1])))
    #t_vac_non_opt = t_sl_opt* (Isp_sl[0]/Isp_sl[1])
    t_sl_non_opt = t_total/(sea_level + ((total - sea_level)*(Isp_vac/Isp_sl)))
    t_vac_opt = t_sl_non_opt * (Isp_vac/Isp_sl)
    return t_sl_non_opt, t_vac_opt

def mass_flow_rate(thrust, specific_impulse):
    """
    Calculate the mass flow rate based on thrust and specific impulse.
    Parameters:
    thrust (float): Thrust in Newtons.
    specific_impulse (float): Specific impulse in seconds.
    Returns:
    float: Mass flow rate in kg/s.

    """
    g0 = 9.81  # Standard gravity in m/s^2
    return thrust / (specific_impulse * g0)

if __name__ == "__main__":
    thrust_sl = total_thrust_to_individual_chamber_thrust(2168252)[0]
    thrust_space = total_thrust_to_individual_chamber_thrust(2168252)[1]

    print(mass_flow_rate(thrust_sl, 393.3471))
    print(mass_flow_rate(thrust_space, 447.9481))
    #its right