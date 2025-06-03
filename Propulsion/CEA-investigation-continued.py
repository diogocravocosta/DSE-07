import scipy as sp
import numpy as np
from rocketcea.cea_obj_w_units import CEA_Obj
import matplotlib.pyplot as plt

C = CEA_Obj( oxName='LOX', fuelName='LH2', pressure_units='Pa', cstar_units='m/s', temperature_units='K', sonic_velocity_units='m/s', enthalpy_units='J/kg', density_units='kg/m^3', specific_heat_units= 'J/kg-K', viscosity_units='millipoise', thermal_cond_units='mcal/cm-K-s')

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
    pressure_ratio = C.get_PcOvPePc(Pc= 6000000, MR=6, eps=, frozen=0, frozenAtThroat=0)
    P_exit = Pc / pressure_ratio
    return P_exit


def area_throat(thrust, Isp, g_0=9.81, Pc):
    mass_flow = thrust / (g_0 * Isp)
    area_throat = (mass_flow * c_star) / Pc
    return area_throat
