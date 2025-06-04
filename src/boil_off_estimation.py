import numpy as np
from sklearn.linear_model import LinearRegression
from sklearn.preprocessing import PolynomialFeatures
import matplotlib.pyplot as plt

#------------------------------------------------
#Input parameters
#------------------------------------------------
# Boil-off data
heat_load_data = [40000, 42500, 45000, 47500, 50000, 52500, 55000]  # W
boil_off_mass = [605, 861, 1470, 1978, 2434, 2855, 3253] #kg

#Geometry parameters
ro = 4.92
ri = 2.46
h = 13.95
phi = np.arctan((ro - ri) / h)  # angle in radians

# Constants
boltzman = 5.67e-8
solar_power = 135311.68  # W
planetary_power = 25795.63  # W
albedo_power = 13604.74  # W

# Material properties
emissivity_ss = 0.8
absorptivity_ss = 0.2

# Mass paramters
m_payload = 15000
m_h2_reentry = 3000
h2_nm = 2.016 #g/mol

#------------------------------------------------
#Functions
#------------------------------------------------

def sa_cone(ro,ri,h):
    l = ((ro-ri)**2 + h**2)**0.5
    return np.pi *( (ro + ri) * l +ri**2 + ro**2)

def volume_cone(h,ro,ri):
    volume = (1 / 3) * np.pi * h * (ro ** 2 + ro * ri + ri ** 2)
    area_proj = (ro + ri)*h/2
    return volume, area_proj

# Calculate the new inner radius based on ullage height
def calculate_cone_param(ro, ri, h, mass):
    vol_cone = volume_cone(h,ro,ri)[0]
    vol_ullage = vol_cone - mass / 71
    h_ull = 3
    ro= ri + h_ull * np.tan(phi)
    while volume_cone(h_ull, ro, ri)[0] < vol_ullage:
        ro= ri + h_ull * np.tan(phi)
        h_ull = h_ull + 0.1  # Increment ullage height until the volume condition is met
    return ro, h_ull

def rad_load(T_tank, T_lh2, emissivity,area_gh2):
    q = emissivity *boltzman*(T_tank**4-T_lh2**4)
    q_load = q * area_gh2
    return q_load

def heat_load(solar_power, planetary_power, albedo_power,area, emissivity, absorptivity):
    # Geometry
    incident_area = 7 * 9 + np.pi / 2 * 0.875 * 3.5
    planetary_flux = planetary_power / incident_area * emissivity
    solar_flux = solar_power / incident_area * absorptivity
    albedo_flux = albedo_power / incident_area * absorptivity

    ro_gh2, h_gh2 = calculate_cone_param(ro, ri, h, m_payload)  # Calculate new outer radius based on ullage height 
    area_gh2 = sa_cone(ro_gh2, ri, h_gh2)  # Calculate the surface area of the cone with the new outer radius
    q_load = rad_load(150, 20, emissivity_ss,area_gh2)  # Example temperatures in Kelvin
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
    while not (0.999 <= f <= 1.001):
        if f < 1:
            n = n -500
        else:
            n =n +500
        f = ((P + a * (n / V) ** 2) * (V - n * b)) / (n * R * T)
        m_gh2 = n * h2_nm / 1000
        iter_count += 1
        if iter_count > 10000:
            raise RuntimeError("Van der Waals solver did not converge")
    return m_gh2

def pres_vanderwaals(n, V, R, T, a, b):
    P = (n * R * T)/(V - n * b) - a * (n / V) ** 2
    return P
#------------------------------------------------
# Calculations
#------------------------------------------------
# During 1st stage of launch is assumed to be negligible.

# During Launch
P1= 10**5
V1 = volume_cone(calculate_cone_param(ro,ri,h,36000)[1], calculate_cone_param(ro,ri,h,40000)[0], ri)[0]  # Volume of the cone
V2 =  volume_cone(calculate_cone_param(ro,ri,h,15000)[1], calculate_cone_param(ro,ri,h,15000)[0], ri)[0]  #m3
P2 = P1 *V1/V2 #pa
T1 = 20
T2 = V2/V1 * T1 #K, temperature after first stage
rhp_gh2 = 0.33 #kg/m3
R = 8.314  # J/(mol*K), universal gas constant
n1 = P1*V1/T1/R
m_vap_h2_1 = n1 * h2_nm/1000  # kg, mass of vaporized hydrogen
n2 = V2*n1/V1
m_vap_h2_2 = n2 * h2_nm/1000  # kg, mass of vaporized hydrogen after first stage

mass_boil_off_launch = m_vap_h2_2 - m_vap_h2_1  # kg, mass of hydrogen vaporized during launch
print(mass_boil_off_launch)
print("Mass of vaporized hydrogen during launch: ", mass_boil_off_launch, "kg")

#During orbit
vol_cone, area_proj = volume_cone(h, ro, ri)
heat_load_h2go = heat_load(solar_power, planetary_power, albedo_power,area_proj, emissivity_ss, absorptivity_ss)
boil_off_specific = linear_regression(heat_load_data, boil_off_mass, heat_load_h2go)
print("Boil off mass: ", boil_off_specific, "kg for the given heat load of: ", heat_load_h2go)

# During Refueling
V1 = volume_cone(calculate_cone_param(ro,ri,h,13500)[1], calculate_cone_param(ro,ri,h,13500)[0], ri)[0]  # Volume of the cone
V2 = volume_cone(calculate_cone_param(ro,ri,h,3000)[1], calculate_cone_param(ro,ri,h,3000)[0], ri)[0]  #m3
P1 = 10**6
T1 = 75 # K, temperature before refueling (temperature of gh2 during venting. should be ideally reset every iteration)
T2 = 50 # K, temperature after refueling (temperature of gh2 after long period of venting. should be ideally reset every iteration)
m_gh2_orbit = vanderwaals(P1,V1, R, T1, 0.2453e-6, 0.02651e-3,h2_nm)
nh2 = m_gh2_orbit / h2_nm*1000
p = pres_vanderwaals(nh2, V2, R, T2, 0.2453e-6, 0.02651e-3)
m_gh2_refuel = vanderwaals(p,V2, R, T2, 0.2453e-6, 0.02651e-3,h2_nm)
if m_gh2_refuel/m_gh2_orbit<1.01:
    print("The pressure will decrease which will not cause boil-off. The vapor mass is: ", m_gh2_refuel, "kg")
    m_boiloff_worst_case = vanderwaals(P1,V2, R, T2, 0.2453e-6, 0.02651e-3,h2_nm) - m_gh2_orbit
    print('Worst case boil off if pressure is held constant: ',m_boiloff_worst_case,'kg at pressure: ', P1/10e5, 'bar')
else:
    print('run complex calc on boilfast to get initial and final pressure and temperature')


#During re-entry
ro_gh2, h_gh2 = calculate_cone_param(ro, ri, h, m_h2_reentry)  # Calculate new outer radius based on ullage height 
area_gh2 = sa_cone(ro_gh2, ri, h_gh2)  
radiation_load = rad_load(200, 20, emissivity_ss,area_gh2)  # Example temperatures in Kelvin
m_boil_off_reentry = linear_regression(heat_load_data, boil_off_mass, radiation_load)
# can be left as is and can be ignored.