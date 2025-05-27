import numpy as np
from pyfluids import Fluid, FluidsList, Input
import matplotlib.pyplot as plt

d_tank = 7 #m
rad_tank = d_tank/2
len_tank = 10.54 #m
A_bestcase = np.pi*(d_tank/2)**2  #m2
A_worstcase = d_tank * len_tank #+ np.pi*rad_tank*rad_tank/4
ellipsoid_cap_SA = 43.809  # m for 0.25*rad tank cap height
cylinder_SA = 2 * np.pi * rad_tank * len_tank  # m
tot_SA = cylinder_SA +2*ellipsoid_cap_SA # m^2
absorptivity_ss = 0.4
emissivity = 0.11
Lv_lh2 = 461000#J/kg
sigma = 5.67*10**(-8)

view_sun = 1-(np.arcsin(6371/(6371+600)))/(np.pi)
flux_sun = 1361 #W/m2
q_solar = flux_sun*A_worstcase*absorptivity_ss*view_sun #W

albedo = 0.3
F = (6371/(6371+600))**2
flux = 1361 #W/m2
q_albedo = flux*albedo*F*A_worstcase *absorptivity_ss

E_ir = 5.67*10**(-8)*255**4 *F#W/m28
q_IR= absorptivity_ss*E_ir* A_worstcase #w

T_internal = 20#353*view_sun + 123*(1-view_sun) #K
q_dissipated = 0#emissivity * 5.67*10**(-8)*(tot_SA)*(T_internal**4)

q_in = q_solar+q_albedo+q_IR
q_total = q_in-q_dissipated
T_ss_plate = (q_in/(emissivity * 5.67*10**(-8)*2*A_worstcase))**0.25
print('Temp of plate', T_ss_plate)

print('q_solar:', q_solar)
print('q_albedo:', q_albedo)
print('q_IR:', q_IR)
print('q_dissipated',q_dissipated)
print('total q is: ',q_total)
 
# define the fluid
_hydrogen = Fluid(FluidsList.Hydrogen).with_state(Input.pressure(1e5), Input.temperature(13.8))  # 1 bar, 13.8 K


# assume the hydrogen is heated to 20 °C
_initial_energy = _hydrogen.internal_energy  # internal energy at 13.8 K

time = 24*3600  # s
print('Tank mass',8000*0.003*(2*np.pi*3.5*13 + 87.62))
print('latent heat boil off',q_total*time/Lv_lh2)

T = np.arange(150,293.15,1)
_boiled_off_mass = []
for i in range(len(T)):
    _final_energy = _hydrogen.heating_to_temperature(temperature=T[i])  # internal energy at 20 °C
    _heating_energy = _final_energy.internal_energy - _initial_energy  # J/ 

    # calculate the boiled-off mass

    _total_heat = q_total * time  # J

    _boiled_off_mass.append(_total_heat / _heating_energy)  # kg

print('Boiled off mass:', _boiled_off_mass[-1])  # kg

q_abs = 10 #w/m2
SA= 385

k_sofi = 0.02 #W/mK
t_sofi = k_sofi*(20-4)/q_abs
rho_sofi = 50 #kg/m3
print('thickness of sofi:', t_sofi,'and mass of sofi:',t_sofi*rho_sofi*SA)


k_mli = 0.001 #W/mK
layer_mli = 40 #number of layers
t_mli = k_mli*(20-4)/q_abs
rho_mli = 0.012928571 #kg/layer/m2
print('thickness of mli:', t_mli,' and mass:',layer_mli*rho_mli*SA)

k_mlivcs = 0.001 #W/mK
layer_mlivcs = 16 #number of layers
t_mlivcs = k_mlivcs*(20-4)/q_abs
rho_mli = 0.012928571 #kg/layer/m2
m_vcs = 250/(2*np.pi*3.5*10) #kg
print('thickness of mli-vcs:', layer_mlivcs*0.0005+,' and mass of mli-vcs:',layer_mlivcs*rho_mli*SA)
print('mass of inner shell' ,2700*0.005*(2*np.pi*3.5*10)*1.1)

'''
plt.plot(T, _boiled_off_mass, label="Boiled-off Mass")
plt.xlabel("Temperature (K)")
plt.ylabel("Boiled-off Mass (kg)")
plt.title("Boiled-off Mass vs Temperature")
plt.legend()
plt.grid()
plt.show()
'''


