import scipy as sp
import numpy as np
from rocketcea.cea_obj_w_units import CEA_Obj
import matplotlib.pyplot as plt

C = CEA_Obj( oxName='LOX', fuelName='LH2', pressure_units='Pa', cstar_units='m/s', temperature_units='K', sonic_velocity_units='m/s', enthalpy_units='J/kg', density_units='kg/m^3', specific_heat_units= 'J/kg-K', viscosity_units='millipoise', thermal_cond_units='mcal/cm-K-s')


#p_e = 15000  # Pa (page 205 SPAD)

#frozen flow approximation (appendix B SPAD)

#T_f = 3400  # K hydrolox combustion flame temperature

#mol = 13.5 # kg/kmol hydrolox combustion product molecular weight (frozen at the chamber M=0)	

#gamma = 1.2  # exhaust product isentropic parameter (frozen at the throat M=1)

#O_F = 6 # oxidizer to fuel ratio 

#------------------------------------------Mixture ratio exploration -----------------------------------------------------------------------------------------

# area_ratios = np.linspace(20, 140, 10)  # Area ratios from 1 to 100

# Pc = 345000

# for e in area_ratios:
#     ispArr = []
#     MR = 1.1
#     mrArr = []
#     while MR < 8:
#         ispArr.append( C(Pc, MR, e ))
#         mrArr.append(MR)
#         MR += 0.05
#     plt.plot(mrArr, ispArr, label='AreaRatio %g, Max ISP %g At MR %g'%(e, np.max(ispArr), np.max(mrArr[np.argmax(ispArr)])) )

# plt.legend(loc='best')
# plt.grid(True)
# plt.title( C.desc )
# plt.xlabel( 'Mixture Ratio' )
# plt.ylabel( 'Isp ODE (sec)' )

# plt.show()

#------------------------------------------ Assuming a Certain O/F Ratio  -----------------------------------------------------------------------------------------


#print(C.get_Cstar(Pc=6000000, MR=6.0))

#chamber_pressures = np.linspace(4*10**6, 10*10**6, 10000)  # Chamber pressures from 1 to 6 MPa


# c_star_arr = []
# MR = 6.0
# for Pc in chamber_pressures:
#     c_star_arr.append(C.get_Cstar(Pc=Pc, MR=MR))

# plt.plot(chamber_pressures, c_star_arr, label='Cstar at MR=%.2f' % MR)
# plt.xlabel('Chamber Pressure (Pa)')
# plt.ylabel('Cstar (m/s)')
# plt.title('Cstar vs Chamber Pressure for O/F Ratio %.2f' % MR)
# plt.grid(True)
# plt.legend()
# plt.show()

def eps_for_isp(MR=6, pc = 6.1*10**6, Isp_desired=450):
    # Initial Parameters
    # MR is the mixture ratio (oxidizer to fuel ratio) in mass
    # Pc is the Chamber pressure in Pa  (default the one for vinci) 
    # Isp_desired Desired specific impulse in seconds

    #parameter initialization

    IspArr = []
    Engine_arr = [['RS-25',20.64*10**6, 6.03, 78, 452.3],['RL10', 4.412*10**6, 5.88, 280, 462],['Vinci', 6.1*10**6, 5.8, 240, 457],['YF-75D', 4.1*10**6, 6.0, 80, 442.6],['YF-79', 4.1*10**6, 6.0, 160, 455.2],['LE-5B', 3.58*10**6, 5, 110, 447], ['LE-9', 10.0*10**6, 5.9, 37, 426]] #database of engines to try and make a correction factor (calculated is ideal isp, RPA gives real ISP)
    area_ratios = np.linspace(20, 200, 1000)
    #print(C.get_Isp(Pc=20.64*10**6, MR=6.03, eps=78))  # Example Isp calculation for a specific area ratio
    #print(450*C.get_Isp(Pc=20.64*10**6, MR=6.03, eps=78)/452.3) #corrected Isp so that if it deviates as much as the RS 25 we can get 450 in reality
    
    #getting the correction factor for the Isp of the engines in the database

    for i in range(len(Engine_arr)):
        Engine_arr[i].append(C.get_Isp(Pc=Engine_arr[i][1], MR=Engine_arr[i][2], eps=Engine_arr[i][3]))
        Engine_arr[i].append((Engine_arr[i][-2]-Engine_arr[i][-1])/Engine_arr[i][-2]*100)
    correction_factor = np.average([Engine_arr[i][-1] for i in range(len(Engine_arr))])
    
    #Calculating the Isp for every area Ratio

    for e in area_ratios:
        isp = C.get_Isp(Pc=pc, MR=MR, eps=e)
        IspArr.append(isp)
    
    #Getting The Engine Area Ratio for the desired Isp

    necesarry_eps = np.interp(desired_isp*(1-correction_factor/100), IspArr, area_ratios) #assume we design for 450 s
    isp_theoretical = C.get_Isp(Pc=pc, MR=MR, eps=necesarry_eps)
    return necesarry_eps, isp_theoretical, correction_factor


def calculate_propmass(struct_ratio, Isp, deltaV, payload):
    #struct_ratio is the structural mass ratio (0.121029372 for the Stoke Space Nova Second Stage)
    #Isp in seconds
    #deltaV in m/s
    #payload in kg

    mu = np.exp(deltaV/(Isp*g0))
    prop_mass = payload*((mu-1)*(1-struct_ratio))/(1-mu*struct_ratio)*1.0126
    return prop_mass

def Mixture_Ratio_Optimizer(pc=6.1*10**6, struct_ratio = 0.121029372, deltaV = 7264.29, payload = 15000, OF_MIN= 4, OF_MAX = 8):
    #Pc in Pa
    #struct_ratio is the structural mass ratio (0.121029372 for the Stoke Space Nova Second Stage)
    #deltaV in m/s
    #payload in kg
    #OF_MIN and OF_MAX are the minimum and maximum mixture ratios to explore (MASS MIXTURE RATIOS)

    g0 = 9.81

    mixture_ratios = np.linspace(OF_MIN, OF_MAX, 10000) 
    Isp_list = []
    prop_masses = []
    for i in range(len(mixture_ratios)):
        Isp_new = C.get_Isp(Pc=pc, MR=mixture_ratios[i], eps=necesarry_eps)/((1-correction_factor/100))
        Isp_list.append(Isp_new)
        prop_masses.append(calculate_propmass(struct_ratio=struct_ratio, Isp=Isp_new, deltaV=deltaV, payload=payload))


    optimal_index = np.argmin(prop_masses)
    optimal_mixture_ratio = mixture_ratios[optimal_index]
    optimal_Isp = Isp_list[optimal_index]
    opt_propellant_mass = prop_masses[optimal_index]

    return optimal_mixture_ratio, optimal_Isp, opt_propellant_mass
    #print("Optimal MR is: "+str(optimal_mixture_ratio)+" with Isp "+str(optimal_Isp)+" and propellant mass "+str(prop_masses[optimal_index]))


    #plt.plot(mixture_ratios, prop_masses)
    #plt.xlabel('Mixture Ratio')
    #plt.ylabel('Propellant Mass')
    #plt.title('Propelant Mass vs Mixture Ratio')
    #plt.grid(True)
    #plt.show()

