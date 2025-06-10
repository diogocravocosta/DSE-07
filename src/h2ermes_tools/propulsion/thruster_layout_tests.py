import numpy as np
import matplotlib.pyplot as plt
import math
from h2ermes_tools.delta_v.landing_burn import landing_burn_no_drag, landing_burn_with_drag
from h2ermes_tools.delta_v.helpers import delta_v_from_final_mass

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
    #print("mdot_sl is", mdot_sl)
    mdot_vac = t_vac[1]/(Isp_vac[1] * 9.80665)  # kg/s
    #print("mdot_vac is", mdot_vac)
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
    t_w_list_landing = np.arange(1.1, 5, 0.1)
    #print(len(t_w_list_landing)) #cd = 0.2
    #delta_v_list_landing = 
    delta_v_list_landing = np.array([311.8233405404127, 253.1417922255299, 219.84431134863922, 198.49630304083158, 183.70535870152332, 172.825038764721, 164.4735159527933, 157.7408310615044, 152.48469649800793, 147.99294491417857, 144.1707497167712, 141.06697615781263, 138.22083983369478, 135.78456999288196, 133.65548132588245, 131.95516764375782, 130.25567020881962, 128.61800461405184, 127.39787685763969, 126.07654615649157, 124.99951831459842, 123.90250608242643, 122.81598649534395, 122.11551048882053, 121.15131372170298, 120.65409381582843, 119.9337562014674, 118.99043705805146, 118.61520792674995, 118.05750838239992, 117.3174246518772, 116.8207428852476, 116.59776557296384, 115.7870573505504, 115.72626171740927, 115.07783991727594, 114.76380329824582, 114.32825023093629, 114.26747964431526])
    t_w_list_vacuum = np.arange(0.6, 1.5, 0.1)
    print(len(t_w_list_vacuum))
    delta_v_list_vacuum = np.array([6210, 5937, 5810, 5739, 5673, 5645, 5629, 5618, 5610])+112+114
    sea_level = np.arange(0, 25, 2)  # Number of sea level thrusters
    vacuum = np.arange(24, -1, -2)  # Number of vacuum thrusters
    Isp_sl = np.array([361.6451, 130.9458])
    Isp_vac = np.array([393.3471, 447.9481])
    m_struct = 20642.21346 # kg
    payload = 18000 # kg
    t_sl_opt = 76026.9 #N t_w_list_landing*m_struct # np.array([61364.9, 19510.9])
    t_sl_non_opt = 82691.5 #N 
    t_vac_opt = 94170.0 #N #t_w_list_vacuum*(m_struct+payload+(185555.4038-20642.21346)) #np.array([66744.2 , 66744.2])
    t_vac_non_opt = 27528.1 #N
    assumed_mtot = 185555.4038 #* 9.81 #N
    #mdot_vac = 15.1937 # kg/s 15.193756344148227
    #mdot_sl = 17.3028 # kg/s 17.30282055269797
    #delta_v_sl = 500 # m/s
    #delta_v_vac = 5800+114+113+30 # m/s
    #delta_v_deorbit = 160 # m/s
    m_props_both = []
    m_props_one = []
    lst_n = []
    #for i in range(len(sea_level)):
    #    m_props.append(sl_vac_thruster_perf(sea_level[i], vacuum[i], Isp_sl, Isp_vac, t_sl, t_vac, m_struct, payload, delta_v_sl, delta_v_vac, delta_v_deorbit)[-1])
    for i in range(len(sea_level)):
        #calculate the thrust to weight ratio of given configuration at sea level and vacuum
        t_sl_all = t_sl_opt * sea_level[i] #+ t_vac_non_opt * vacuum[i]
        t_vac_all = t_sl_non_opt * sea_level[i] + t_vac_opt * vacuum[i]
        t_sl_only = t_sl_opt * sea_level[i]
        t_vac_only = t_vac_opt * vacuum[i]

        t_w_sl_all = float(np.round(t_sl_all / ((m_struct + payload)*9.81),1))
        t_w_vac_all = float(np.round(t_vac_all / ((assumed_mtot)*9.81),1))
        t_w_sl_only = float(np.round(t_sl_only / ((m_struct + payload)*9.81),1))
        t_w_vac_only = float(np.round(t_vac_only / ((assumed_mtot)*9.81),1))
        if sea_level[i] == 8:
            print(t_w_sl_all, t_w_vac_all, t_w_sl_only, t_w_vac_only)
        if np.any(np.isclose(t_w_sl_all, t_w_list_landing, atol=1e-08)) and np.any(np.isclose(t_w_vac_all, t_w_list_vacuum, atol=1e-08)) and np.any(np.isclose(t_w_sl_only, t_w_list_landing, atol=1e-08)) and np.any(np.isclose(t_w_vac_only, t_w_list_vacuum, atol=1e-08)):
            lst_n.append(sea_level[i])
            delta_v_sl_all = delta_v_list_landing[np.argmin(np.abs(t_w_list_landing - t_w_sl_all))]
            delta_v_vac_all = delta_v_list_vacuum[np.argmin(np.abs(t_w_list_vacuum - t_w_vac_all))]
            delta_v_sl_only = delta_v_list_landing[np.argmin(np.abs(t_w_list_landing - t_w_sl_only))]
            delta_v_vac_only = delta_v_list_vacuum[np.argmin(np.abs(t_w_list_vacuum - t_w_vac_only))]
            delta_v_deorbit = 160 # m/s, deorbit burn delta V
            #print(delta_v_sl_all, delta_v_vac_all, delta_v_sl_only, delta_v_vac_only)
            both_props = sl_vac_thruster_perf(sea_level[i], vacuum[i], Isp_sl, Isp_vac, [t_sl_opt, t_vac_non_opt], [t_sl_non_opt, t_vac_opt], m_struct, payload, delta_v_sl_all, delta_v_vac_all, delta_v_deorbit)
            m_props_both.append(both_props[-1] - both_props[0]+ mprop_calc(delta_v_sl_only, Isp_sl[0], m_struct))  # Total propellant mass for both thrusters
            m_props_one.append(mprop_calc(delta_v_sl_only, Isp_sl[0], m_struct) + mprop_calc(delta_v_vac_only, Isp_vac[1], m_struct + mprop_calc(delta_v_sl_only, Isp_sl[0], m_struct)))
            
            
        #print(t_w_sl_all, t_w_vac_all, t_w_sl_only, t_w_vac_only)

    #print the best configuration with either one or both thrusters which is included in lst_n
    
    print(lst_n)

    
    #graph the results
    print(m_props_both)
    print(m_props_one)
    plt.plot(lst_n, m_props_both, label='Both Types of Thrusters in Sea Level and Vacuum')
    plt.plot(lst_n, m_props_one, label='Both Types of Thrusters in Vacuum only')
    plt.xlabel('Number of Sea Level Thrusters')
    plt.ylabel('Total Propellant Mass (kg)')
    plt.title('Total Propellant Mass vs Number of Sea Level Thrusters')
    plt.grid()
    plt.legend()
    plt.show()


    #print("The smallest propellant mass is achieved with ", sea_level[np.argmin(m_props)], "sea level thrusters and", vacuum[np.argmin(m_props)], "vacuum thrusters.")
    #print("The total propellant mass is ", min(m_props), "kg.")
    #print("The propellant mass for the sea level burn is ", sl_vac_thruster_perf(sea_level[np.argmin(m_props)], vacuum[np.argmin(m_props)], Isp_sl, Isp_vac, t_sl, t_vac, m_struct, payload, delta_v_sl, delta_v_vac, delta_v_deorbit)[0], "kg.")
    #print("The propellant mass for the deorbit burn is ", sl_vac_thruster_perf(sea_level[np.argmin(m_props)], vacuum[np.argmin(m_props)], Isp_sl, Isp_vac, t_sl, t_vac, m_struct, payload, delta_v_sl, delta_v_vac, delta_v_deorbit)[1], "kg.")
    #print("The propellant mass for the vacuum burn is ", sl_vac_thruster_perf(sea_level[np.argmin(m_props)], vacuum[np.argmin(m_props)], Isp_sl, Isp_vac, t_sl, t_vac, m_struct, payload, delta_v_sl, delta_v_vac, delta_v_deorbit)[2], "kg.")
    #print("The total propellant mass is ", sl_vac_thruster_perf(sea_level[np.argmin(m_props)], vacuum[np.argmin(m_props)], Isp_sl, Isp_vac, t_sl, t_vac, m_struct, payload, delta_v_sl, delta_v_vac, delta_v_deorbit)[-1], "kg.")

    # Plotting the results
    #print("Total Massflow is ", str(15.193756344148227* vacuum[np.argmin(m_props)] + 17.30282055269797*sea_level[np.argmin(m_props)]), "kg/s")

    #plt.plot(sea_level, m_props, label='Total Propellant Mass')
    #plt.xlabel('Number of Sea Level Thrusters')
    #plt.ylabel('Total Propellant Mass (kg)')
    #plt.title('Total Propellant Mass vs Number of Sea Level Thrusters')
    #plt.grid()
    #plt.legend()
    
    # Add annotation
    #min_index = np.argmin(m_props)
    #plt.annotate(f'Min Prop Mass:\n{sea_level[min_index]} SL, {vacuum[min_index]} Vac',
    #             xy=(sea_level[min_index], m_props[min_index]),
    #             xytext=(sea_level[min_index] + 1, m_props[min_index] + 1000),  # Adjust text position as needed
    #             arrowprops=dict(facecolor='black', shrink=0.05),
    #             )

    #plt.show()
