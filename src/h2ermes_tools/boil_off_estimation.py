import numpy as np
from sklearn.linear_model import LinearRegression
import matplotlib.pyplot as plt
import pytest
from data import material as mat
import data.constants as cn

#------------------------------------------------
#Functions
#------------------------------------------------

def sa_cone(ro,ri,h):
    l = ((ro-ri)**2 + h**2)**0.5
    return np.pi *( (ro + ri) * l +ri**2 + ro**2)

def volume_cone(h,ro,ri):
    volume = (1 / 3) * np.pi * h * (ro ** 2 + ro * ri + ri ** 2) + 2/3*np.pi * ri**3/4 + 2/3*np.pi * ro**3/4
    area_proj_cylinder = (ro + ri)*h
    area_proj = area_proj_cylinder + np.pi*ri**2/4+np.pi*ro**2/4
    return volume, area_proj

def height_calc_cone(volume):
    h = volume*3/np.pi/(ro ** 2 + ro * ri + ri ** 2)
    return h
# Calculate the new inner radius based on ullage height
def calculate_cone_param(ro, ri, h, mass,rho_lh2):
    vol_cone = volume_cone(h,ro,ri)[0]
    vol_ullage = vol_cone - mass / rho_lh2
    h_ull = 3
    ro= ri + h_ull * np.tan(phi)
    while volume_cone(h_ull, ro, ri)[0] < vol_ullage:
        ro= ri + h_ull * np.tan(phi)
        h_ull = h_ull + 0.1  # Increment ullage height until the volume condition is met
    return ro, h_ull

def rad_load(T_tank, T_lh2, material,area_gh2):
    q = material.eps *boltzman*(T_tank**4-T_lh2**4)
    q_load = q * area_gh2
    return q_load

def heat_load(solar_power, planetary_power, albedo_power,area, material,rho_lh2_20k, m_payload):
    # Geometry
    incident_area = 7 * 9 + np.pi * 0.875 * 3.5
    planetary_flux = planetary_power / incident_area * material.eps
    solar_flux = solar_power / incident_area * material.abs
    albedo_flux = albedo_power / incident_area * material.abs

    ro_gh2, h_gh2 = calculate_cone_param(ro, ri, h, m_payload,rho_lh2_20k)  # Calculate new outer radius based on ullage height 
    area_gh2 = sa_cone(ro_gh2, ri, h_gh2)  # Calculate the surface area of the cone with the new outer radius
    q_load = rad_load(150, 20, material,area_gh2)  # Example temperatures in Kelvin
    heat_load = (solar_flux + planetary_flux + albedo_flux) * area + q_load
    return heat_load

def linear_regression(heat_load_data, boil_off_mass, heat_load_h2go):
    """
    Perform linear regression to estimate boil-off mass based on heat load data.
    
    Parameters:
    heat_load_data (list): List of heat load values in W.
    boil_off_mass (list): List of corresponding boil-off mass values in kg.
    
    Returns:
    tuple: Coefficients (a, b) of the linear regression model.
    """
    # Convert to numpy arrays
    X = np.array(heat_load_data).reshape(-1, 1)
    y = np.array(boil_off_mass)

    # Perform linear regression
    model = LinearRegression()
    model.fit(X, y)

    a = model.coef_[0]
    b = model.intercept_
    boil_off_specific = a*heat_load_h2go + b
    return boil_off_specific

def vanderwaals(P, V, R, T, a, b, h2_nm):
    n = 600000
    f = ((P + a * (n / V) ** 2) * (V - n * b)) / (n * R * T)
    iter_count = 0
    while not (0.99 <= f <= 1.01):
        if f < 1:
            n = n -500
        else:
            n =n +500
        f = ((P + a * (n / V) ** 2) * (V - n * b)) / (n * R * T)
        m_gh2 = n * h2_nm / 1000
        iter_count += 1
        if iter_count > 1000000:
            raise RuntimeError("Van der Waals solver did not converge")
    return m_gh2

def pres_vanderwaals(n, V, R, T, a, b):
    P = (n * R * T)/(V - n * b) - a * (n / V) ** 2
    return P

def boil_off_launch(P1,T1,R,m_h2_tot,m_pl,ro,ri,h,rho_lh2_20k):
    V1_vapor = 0.1*volume_cone(h,ro,ri)[0]#volume_cone(calculate_cone_param(ro,ri,h,m_h2_tot,rho_lh2_20k)[1], calculate_cone_param(ro,ri,h,m_h2_tot,rho_lh2_20k)[0], ri)[0]  # Volume of the cone
    V2_vapor =  volume_cone(calculate_cone_param(ro,ri,h,m_pl,rho_lh2_20k)[1], calculate_cone_param(ro,ri,h,m_pl,rho_lh2_20k)[0], ri)[0]  #m3
    n1 = P1*V1_vapor/T1/R
    m_vap_h2_launchstart = vanderwaals(P1, V1_vapor, R, T1, a, b, h2_nm)

    m_vap_h2_launchend = vanderwaals(P1, V2_vapor, R, T1, a, b, h2_nm)

    mass_boil_off_launch = m_vap_h2_launchend - m_vap_h2_launchstart 
    # print("Boil off during launch (due to rapid vaporization): ", mass_boil_off_launch, "kg for the start volume of: ",volume_cone(h,ro,ri)[0],"m3 and projected area of",volume_cone(h,ro,ri)[1])
    return mass_boil_off_launch

def orbit_boil_off(h,ro,ri, solar_power, planetary_power, albedo_power, material, heat_load_data, boil_off_mass,rho_lh2_20k, m_payload):
    area_proj = volume_cone(h, ro, ri)[1]
    heat_load_h2go = heat_load(solar_power, planetary_power, albedo_power,area_proj, material,rho_lh2_20k, m_payload)
    boil_off_specific = linear_regression(heat_load_data, boil_off_mass, heat_load_h2go)
    # print("Boil off during orbit (due to external heat sources): ", boil_off_specific, "kg for the given heat load of: ", heat_load_h2go,"W")
    return boil_off_specific

def boil_off_refueling(p_vent, T_vapor,T_vapor_refuel, a,b,R,h2_nm,rho_lh2_20k,m_h2_reentry,m_h2_dock):
    V1 = volume_cone(calculate_cone_param(ro,ri,h,m_h2_dock,rho_lh2_20k)[1], calculate_cone_param(ro,ri,h,m_h2_dock,rho_lh2_20k)[0], ri)[0]  # Volume of the cone
    V2 = volume_cone(calculate_cone_param(ro,ri,h,m_h2_reentry,rho_lh2_20k)[1], calculate_cone_param(ro,ri,h,m_h2_reentry,rho_lh2_20k)[0], ri)[0]  #m3
    P1 = p_vent
    T1 = T_vapor # K, temperature before refueling (temperature of gh2 during venting. should be ideally reset every iteration)
    T2 = T_vapor_refuel# K, temperature after refueling (temperature of gh2 after long period of venting. should be ideally reset every iteration)
    m_gh2_orbit = vanderwaals(P1,V1, R, T1, a, b,h2_nm)
    nh2 = m_gh2_orbit / h2_nm*1000
    p = pres_vanderwaals(nh2, V2, R, T2, a, b)
    m_gh2_refuel = vanderwaals(p,V2, R, T2,a, b,h2_nm)

    print("No boil off is expected in this region as pressure will drop from ",P1/10e5, "bar to ",p/10e5,"bar")
    m_boiloff_worst_case = vanderwaals(P1,V2, R, T2, a, b,h2_nm) - m_gh2_orbit

    print('Worst case boil off if pressure is held constant ',m_boiloff_worst_case,'kg at pressure: ', P1/10e5, 'bar. The change in mass is ',m_gh2_refuel/m_gh2_orbit)

    return m_boiloff_worst_case

def boiloff_reentry(ro,ri,h,m_h2_reentry,material,T_skin_reentry,T_lh2,rho_lh2_30k):
    ro_gh2, h_gh2 = calculate_cone_param(ro, ri, h, m_h2_reentry,rho_lh2_30k)  # Calculate new outer radius based on ullage height 
    area_gh2 = sa_cone(ro_gh2, ri, h_gh2)  
    radiation_load = rad_load(T_skin_reentry, T_lh2, material,area_gh2)  # Example temperatures in Kelvin
    m_boil_off_reentry = linear_regression(heat_load_data, boil_off_mass, radiation_load)

    print('Boil off in reentry due to tank wall heating up: ',m_boil_off_reentry,'kg for radiation load: ',radiation_load,"W")
    return m_boil_off_reentry

def total_boil_off_h2(m_h2_tot,
                      m_payload,
                      m_coolant,
                      m_boil_off,
                      middle_radius,
                      top_radius,
                      height,
                      material,
                      T_gh2_launch=20,
                      heat_load_data=[44000,40000,50000,47000],
                      boil_off_mass_data=[16381-14200, 16381-15040,16381-13000,16381-13604],
                      p_vent=10e6,
                      T_vapor=75,
                      T_vapor_refuel=50,
                      a=0.2453e-6,
                      b=0.02651e-3,
                      h2_nm=2,
                      solar_power=135311.68,
                      planetary_power=25795.63,
                      albedo_power=13604.74,
                      rho_lh2_20k=71,
                      worst_case=False,
                      P_min_vapor=3e5,
                      ):
    m_payload_coolant_boil_off = m_payload + m_coolant + m_boil_off

    # During launch
    total_boil_off = 0
    mass_boil_off_launch = boil_off_launch(
        P_min_vapor,
        T_gh2_launch,
        cn.R_star/1000,
        m_h2_tot,
        m_payload_coolant_boil_off,
        middle_radius,
        top_radius,
        height,
        rho_lh2_20k)
    total_boil_off += mass_boil_off_launch

    # During orbit
    boil_off_specific = orbit_boil_off(
        height,
        middle_radius,
        top_radius,
        solar_power,
        planetary_power,
        albedo_power,
        material,
        heat_load_data,
        boil_off_mass_data,
        rho_lh2_20k,
        m_payload_coolant_boil_off
    )
    total_boil_off += boil_off_specific  # kg, total boil-off mass during orbit

    transfer_boil_off = 0.05/0.95*m_payload  # kg, transfer boil-off mass during orbit
    total_boil_off += transfer_boil_off
    # During Refueling
    if worst_case is True:
        m_boiloff_worst_case = boil_off_refueling(
            p_vent,
            T_vapor,
            T_vapor_refuel,
            a,
            b,
            R,
            h2_nm,
            rho_lh2_20k,
            m_coolant,
            m_payload + m_coolant
        )
        total_boil_off +=  m_boiloff_worst_case
    # print('Total boil off of LH2 is: ',total_boil_off,"kg")
    # print('New payload mass is: ', h2_depot + total_boil_off + m_h2_coolant, "kg")
    return total_boil_off
#------------------------------------------------
# Calculations
#------------------------------------------------
if __name__ =='__main__':
    #------------------------------------------------
    #Input parameters
    #------------------------------------------------

    #Hydrogen Parameters
    rho_lh2_20k = 71
    rho_lh2_30K = 50

    # Mass paramters
    h2_nm = 2.016 #g/mol
    total_boil_off = 0

    m_payload = 10000
    m_coolant = 3000
    m_boil_off = 2500

    m_prop_h2 = 150000/7
    m_h2_tot = m_prop_h2 + m_payload + m_coolant + m_boil_off

    #Geometry parameters
    ro = 4.92
    ri = 2.46
    volume = m_h2_tot / rho_lh2_20k/0.9
    phi = 10
    h = 12-ro/4-ri/4
    # Pressure Parameters
    p_vent = 10e6 #pa
    P_min_vapor= 3e5 #pa

    #Temperature parameters
    T_vapor = 75
    T_vapor_refuel = 50
    T_skin_reentry = 200
    T_lh2 = 20
    T_gh2_launch = 20

    # Material properties
    material = mat.Material(absorptivity=0.2,
                            emissivity=0.08)

    #------------------------------------------------
    # Constants
    #------------------------------------------------
    boltzman = 5.67e-8
    solar_power = 135311.68  # W
    planetary_power = 25795.63  # W
    albedo_power = 13604.74  # W
    R = 8.314
    a = 0.2453e-6
    b = 0.02651e-3
    h2_nm = 2 #g/mol

    # Boil-off data40000, 42500, 45000, 47500, 50000, 52500, 55000: 861, 1470, 1978, 2434, 2855, 3253
    heat_load_data = [44000,40000,50000,47000]  # W
    boil_off_mass = [16381-14200, 16381-15040,16381-13000,16381-13604] #kg

    # Conditions
    worst_case = False

    total_boil_off_mass = total_boil_off_h2(
        m_h2_tot,
        m_payload,
        m_coolant,
        m_boil_off,
        ro,
        ri,
        h,
        material,
    )

    print("Total boil off mass of LH2 is: ", total_boil_off_mass, "kg")
