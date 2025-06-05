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

area_ratios = np.linspace(20, 100, 10)  # Area ratios from 1 to 100

Pc = 600000

for e in area_ratios:
    ispArr = []
    MR = 1.1
    mrArr = []
    while MR < 8:
        ispArr.append( C(Pc, MR, e ))
        mrArr.append(MR)
        MR += 0.05
    plt.plot(mrArr, ispArr, label='AreaRatio %g, Max ISP %g At MR %g'%(e, np.max(ispArr), np.max(mrArr[np.argmax(ispArr)])) )

plt.legend(loc='best')
plt.grid(True)
plt.title( C.desc )
plt.xlabel( 'Mixture Ratio' )
plt.ylabel( 'Isp ODE (sec)' )

plt.show()