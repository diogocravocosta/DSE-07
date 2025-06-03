import scipy as sp
import numpy as np
from rocketcea.cea_obj_w_units import CEA_Obj
import rocketcea 
import matplotlib.pyplot as plt



C = CEA_Obj( oxName='LOX', fuelName='LH2', pressure_units='Pa', cstar_units='m/s', temperature_units='K', sonic_velocity_units='m/s', enthalpy_units='J/kg', density_units='kg/m^3', specific_heat_units= 'J/kg-K', viscosity_units='millipoise', thermal_cond_units='mcal/cm-K-s')

Pc = 6000000
MR = 6.0
eps = 120

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
    pressure_ratio = C.get_PcOvPe(Pc= 6000000, MR=6, eps=eps, frozen=0, frozenAtThroat=0)
    P_exit = Pc / pressure_ratio
    return P_exit

Pc = 6000000
print(Isp)

def area_throat(thrust, Isp = Isp, g_0=9.81, Pc= Pc, MR=MR):
    mass_flow = thrust / (g_0 * Isp)
    area_throat = (mass_flow * C.get_Cstar(Pc, MR)) / Pc
    return area_throat

print("exit presssure is", obtain_exit_pressure(Pc, MR, eps))
print("exit area is", area_throat(thrust, Isp=Isp, Pc=Pc, MR=MR)*eps)
print("exit diameter is", (((area_throat(thrust, Isp=Isp, Pc=Pc, MR=MR)*eps))/np.pi)**0.5 * 2)
print("area of the throat", area_throat(thrust, Isp=Isp, Pc=Pc, MR=MR))
