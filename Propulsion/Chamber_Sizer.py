import scipy as sp
import numpy as np
from rocketcea.cea_obj_w_units import CEA_Obj
import rocketcea 
import matplotlib.pyplot as plt



C = CEA_Obj( oxName='LOX', fuelName='LH2', pressure_units='Pa', cstar_units='m/s', temperature_units='K', sonic_velocity_units='m/s', enthalpy_units='J/kg', density_units='kg/m^3', specific_heat_units= 'J/kg-K', viscosity_units='millipoise', thermal_cond_units='mcal/cm-K-s')

Pc = 6000000
MR = 6.0
eps = 80

Isp = C.get_Isp(Pc, MR, eps=eps)
thrust = 66700 #in N

def obtain_cstar(Pc, MR):
    """
    Obtain the characteristic velocity (C*) for given chamber pressure and mixture ratio.
    
    Parameters:
    Pc (float): Chamber pressure in Pascals.
    MR (float): Mixture ratio (mass of oxidizer to mass of fuel).
    
    Returns:
    float: Characteristic velocity (C*) in m/s.
    """
    return C.get_Cstar(Pc, MR)


def obtain_exit_pressure(Pc, MR, eps, frozen=0, frozenAtThroat=0): 
    pressure_ratio = C.get_PcOvPe(Pc= Pc, MR=MR, eps=eps, frozen=0, frozenAtThroat=0)
    P_exit = Pc / pressure_ratio
    return P_exit


#print(Isp)

def throat_geometry(thrust, Isp = Isp, Pc= Pc, MR=MR):
    g_0=9.81
    mass_flow = thrust / (g_0 * Isp)
    area_throat = (mass_flow * C.get_Cstar(Pc, MR)) / Pc
    diameter_throat = (4 * area_throat / np.pi) ** 0.5
    return area_throat, diameter_throat

def exit_geometry(thrust, Isp=Isp, Pc=Pc, MR=MR, eps=eps):
    area = throat_geometry(thrust, Isp=Isp, Pc=Pc, MR=MR)[0] * eps
    diameter = (4 * area / np.pi) ** 0.5
    return area, diameter

print("exit pressure is", obtain_exit_pressure(Pc, MR, eps))
print("exit area is", throat_geometry(thrust, Isp=Isp, Pc=Pc, MR=MR)[0]*eps)
print("exit diameter is", (((area_throat(thrust, Isp=Isp, Pc=Pc, MR=MR)*eps))/np.pi)**0.5 * 2)
print("area of the throat", area_throat(thrust, Isp=Isp, Pc=Pc, MR=MR))
print("throat diameter is", (((area_throat(thrust, Isp=Isp, Pc=Pc, MR=MR)))/np.pi)**0.5 * 2)

def chamber_geometry(thrust, Isp=Isp, Pc=Pc, MR=MR, eps=eps, L_star=0.76):
    #combustion chamber sizing
    #historical characteristic length between 0.76 to 1.02 m, lower more aggressive, higher more conservative (not really physical limits)
    #L_star = 0.76-1.02  # m, characteristic length of the combustion chamber, historical value high conservative, low less conservative, author from SPAD says low can be taken due to progress of comb chambers
        
    a_t = throat_geometry(thrust, Isp=Isp, Pc=Pc, MR=MR)[0]
    D_t = throat_geometry(thrust, Isp=Isp, Pc=Pc, MR=MR)[1]
    
    v_c = L_star*a_t

    #mach in the chamber between 0.1 and 0.6, 0.1 to be conservative

    gamma = C.get_Chamber_MolWt_gamma(Pc=Pc, MR=MR, eps=eps)[-1]
    chamber_contraction_ratio = 8*(D_t*100)**-0.6+1.25 #statistical relation #1/mach *((2/(gamma+1)*(1+(gamma-1)/2*mach**2))**((gamma+1)/(2*(gamma-1))))
    combustion_length = L_star / chamber_contraction_ratio
    combustion_area = a_t * chamber_contraction_ratio


    #lc/dc ratios from 0.5 to 2.5

    D_c = np.sqrt(4*combustion_area/np.pi)
    #print(D_c)
    #print(combustion_area)
    #print(combustion_length)
    #print(combustion_length/ (np.sqrt(4*combustion_area/np.pi)))

    D_e = np.sqrt(4*eps*a_t/np.pi)

    return combustion_length, D_c, a_t, D_t, D_e, chamber_contraction_ratio, combustion_area



