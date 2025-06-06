import numpy as np
import matplotlib.pyplot as plt
import math

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

#m_props.append(sl_vac_thruster_perf(sea_level[i], vacuum[i], Isp_sl, Isp_vac, t_sl, t_vac, mdot_sl, mdot_vac, m_struct, payload, delta_v_sl, delta_v_vac, delta_v_deorbit)[-1])
def sl_vac_thruster_perf(n_sl, n_vac, Isp_sl, Isp_vac, t_sl, t_vac, m_struct, payload, delta_v_sl, delta_v_vac, delta_v_deorbit):
    """
    Calculate the performance of the launch vehicle in sea level and vacuum conditions for a given thruster config.

    Args:
        n_sl (int): Number of thrusters in sea level.
        n_vac (int): Number of thrusters in vacuum.
        Isp_sl (array): Specific impulse in sea level. Isp_sl[0] for the sea level optimized chambers and Isp_sl[1] for the vacuum optimized chambers.
        Isp_vac (array): Specific impulse in vacuum. Isp_vac[0] for the sea level optimized chambers and Isp_vac[1] for the vacuum optimized chambers.
        t_sl (array): Thrust for each thruster at sea level. t_sl[0] for the sea level optimized chambers and t_sl[1] for the vacuum optimized chambers.
        t_vac (array): Thrust for each thruster in vacuum. t_vac[0] for the sea level optimized chambers and t_vac[1] for the vacuum optimized chambers.
        mdot_sl (float): Mass flow rate for the sea level thrusters.
        mdot_vac (float): Mass flow rate  for the vacuum thrusters.
        m_struct (float): Structural mass in kg.
        payload (float): Payload mass in kg.
        delta_v_sl (float): Delta V in sea level in m/s. (landing)
        delta_v_vac (float): Delta V in vacuum in m/s.
        delta_v_deorbit (float): Delta V for deorbit burn in m/s.

    Returns:
        
    """
    mdot_sl = t_vac[0]/(Isp_vac[0] * 9.80665)  # kg/s
    print("mdot_sl is", mdot_sl)
    mdot_vac = t_vac[1]/(Isp_vac[1] * 9.80665)  # kg/s
    print("mdot_vac is", mdot_vac)
    Isp_compound_sl = ((t_sl[0]*n_sl) + (t_sl[1]*n_vac))/(((mdot_vac*n_vac)+(mdot_sl*n_sl))*9.80665)
    Isp_compound_vac = ((t_vac[0]*n_sl) + (t_vac[1]*n_vac))/(((mdot_vac*n_vac)+(mdot_sl*n_sl))*9.80665)
    # calculate the mass of propellant necessary given a numver of thrusters
    m_prop_sl = mprop_calc(delta_v_sl, Isp_compound_sl, m_struct)
    m_prop_deorbit = mprop_calc(delta_v_deorbit, Isp_compound_vac, m_struct+m_prop_sl)
    m_prop_vac = mprop_calc(delta_v_vac, Isp_compound_vac, m_struct+m_prop_sl+m_prop_deorbit)
    m_final = m_struct + m_prop_sl + m_prop_deorbit + m_prop_vac + payload
    m_prop_total = m_prop_sl + m_prop_deorbit + m_prop_vac
    return m_prop_sl,m_prop_deorbit, m_prop_vac, m_prop_total

if __name__ == "__main__":
    sea_level = [1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,23]
    vacuum = [23,22,21,20,19,18,17,16,15,14,13,12,11,10,9,8,7,6,5,4,3,2,1]
    Isp_sl = np.array([361.6451, 130.9458])
    Isp_vac = np.array([393.3471, 447.9481])
    t_sl = np.array([61364.9, 19510.9])
    t_vac = np.array([66744.2 , 66744.2])
    #mdot_vac = 15.1937 # kg/s 15.193756344148227
    #mdot_sl = 17.3028 # kg/s 17.30282055269797
    m_struct = 20642.21346 # kg
    payload = 18000 # kg
    delta_v_sl = 500 # m/s
    delta_v_vac = 5800+114+113+30 # m/s
    delta_v_deorbit = 160 # m/s
    m_props = []
    for i in range(len(sea_level)):
        m_props.append(sl_vac_thruster_perf(sea_level[i], vacuum[i], Isp_sl, Isp_vac, t_sl, t_vac, m_struct, payload, delta_v_sl, delta_v_vac, delta_v_deorbit)[-1])
    
    print("The smallest propellant mass is achieved with ", sea_level[np.argmin(m_props)], "sea level thrusters and", vacuum[np.argmin(m_props)], "vacuum thrusters.")
    print("The total propellant mass is ", min(m_props), "kg.")
    print("The propellant mass for the sea level burn is ", sl_vac_thruster_perf(sea_level[np.argmin(m_props)], vacuum[np.argmin(m_props)], Isp_sl, Isp_vac, t_sl, t_vac, m_struct, payload, delta_v_sl, delta_v_vac, delta_v_deorbit)[0], "kg.")
    print("The propellant mass for the deorbit burn is ", sl_vac_thruster_perf(sea_level[np.argmin(m_props)], vacuum[np.argmin(m_props)], Isp_sl, Isp_vac, t_sl, t_vac, m_struct, payload, delta_v_sl, delta_v_vac, delta_v_deorbit)[1], "kg.")
    print("The propellant mass for the vacuum burn is ", sl_vac_thruster_perf(sea_level[np.argmin(m_props)], vacuum[np.argmin(m_props)], Isp_sl, Isp_vac, t_sl, t_vac, m_struct, payload, delta_v_sl, delta_v_vac, delta_v_deorbit)[2], "kg.")
    print("The total propellant mass is ", sl_vac_thruster_perf(sea_level[np.argmin(m_props)], vacuum[np.argmin(m_props)], Isp_sl, Isp_vac, t_sl, t_vac, m_struct, payload, delta_v_sl, delta_v_vac, delta_v_deorbit)[-1], "kg.")

    # Plotting the results
    print("Total Massflow is ", str(15.193756344148227* vacuum[np.argmin(m_props)] + 17.30282055269797*sea_level[np.argmin(m_props)]), "kg/s")

    plt.plot(sea_level, m_props, label='Total Propellant Mass')
    plt.xlabel('Number of Sea Level Thrusters')
    plt.ylabel('Total Propellant Mass (kg)')
    plt.title('Total Propellant Mass vs Number of Sea Level Thrusters')
    plt.grid()
    plt.legend()
    plt.show()
