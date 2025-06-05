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

# Current Inputs

MR = 6

pc = 6.1*10**6  # Chamber pressure in Pa  (vinci) 

desired_isp = 450  # Desired specific impulse in seconds

IspArr = []

Engine_arr = [['RS-25',20.64*10**6, 6.03, 78, 452.3],['RL10', 4.412*10**6, 5.88, 280, 462],['Vinci', 6.1*10**6, 5.8, 240, 457],['YF-75D', 4.1*10**6, 6.0, 80, 442.6],['YF-79', 4.1*10**6, 6.0, 160, 455.2],['LE-5B', 3.58*10**6, 5, 110, 447], ['LE-9', 10.0*10**6, 5.9, 37, 426]]

#calculations

area_ratios = np.linspace(20, 200, 1000)

#print(C.get_Isp(Pc=20.64*10**6, MR=6.03, eps=78))  # Example Isp calculation for a specific area ratio
#print(450*C.get_Isp(Pc=20.64*10**6, MR=6.03, eps=78)/452.3) #corrected Isp so that if it deviates as much as the RS 25 we can get 450 in reality

for i in range(len(Engine_arr)):
    Engine_arr[i].append(C.get_Isp(Pc=Engine_arr[i][1], MR=Engine_arr[i][2], eps=Engine_arr[i][3]))
    Engine_arr[i].append((Engine_arr[i][-2]-Engine_arr[i][-1])/Engine_arr[i][-2]*100)

correction_factor = np.average([Engine_arr[i][-1] for i in range(len(Engine_arr))])
#print(correction_factor)

for e in area_ratios:
    isp = C.get_Isp(Pc=pc, MR=MR, eps=e)
    IspArr.append(isp)

# plt.plot(area_ratios, IspArr, label='Isp at MR=%.2f' % MR)
# plt.xlabel('Area Ratio')
# plt.ylabel('Isp (sec)')
# plt.title('Isp vs Area Ratio for O/F Ratio %.2f' % MR)
# plt.grid(True)
# plt.legend()
# plt.show()


#print("The necessary area ratio to achieve 450 seconds of Isp is: %.2f" % (necesarry_eps))
#Isp_chosen = 450*(1-correction_factor/100)

necesarry_eps = np.interp(desired_isp*(1-correction_factor/100), IspArr, area_ratios) #assume we design for 450 s
isp_theoretical = C.get_Isp(Pc=pc, MR=MR, eps=necesarry_eps)



