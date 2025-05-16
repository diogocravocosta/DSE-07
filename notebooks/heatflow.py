import numpy as np
from pyfluids import Fluid, FluidsList, Input
import matplotlib.pyplot as plt

d_tank = 7 #m
len_tank = 10 #m
A_bestcase = np.pi*(d_tank/2)**2  #m2
A_worstcase = d_tank * len_tank
absorptivity = 0.6
view_sun = 1-(np.arcsin(6371/(6371+600)))/(np.pi)
flux_sun = 1361 #W/m2
q_solar = flux_sun*A_worstcase*absorptivity#*view_sun #W

albedo = 0.3
F = (6371/(6371+600))**2
flux = 1361 #W/m2
q_albedo = flux*albedo*A_worstcase *F

emissivity = 1
E_ir = 237 *(6371/(6371+600))**2#W/m28
q_IR= emissivity * A_worstcase*E_ir #w

T_internal = 353#353*view_sun + 123*(1-view_sun) #K

emissivity = 0.2
T_space = 4 #K
q_dissipated = emissivity * 5.67*10**(-8)*(np.pi*3.5*10)*(T_internal**4)

q_total = q_solar+q_albedo+q_IR-q_dissipated

print('q_solar:', q_solar)
print('q_albedo:', q_albedo)
print('q_IR:', q_IR)
print('q_dissipated',q_dissipated)
print('total q is: ',q_total)

# define the fluid
_hydrogen = Fluid(FluidsList.Hydrogen).with_state(Input.pressure(1e5), Input.temperature(13.8))  # 1 bar, 13.8 K

# assume the hydrogen is heated to 20 °C
_initial_energy = _hydrogen.internal_energy  # internal energy at 13.8 K

len_tank = 10.54  # m
rad_tank = 3.5
ellipsoid_cap_SA = 43.809  # m for 0.25*rad tank cap height
cylinder_SA = 2 * np.pi * rad_tank * len_tank  # m
_surface_area = cylinder_SA  # m^2
ellipsoid_vol = 22.45
print('total_SA:', _surface_area) 

_total_heat_flux =  q_total/_surface_area
_time = 24*3600  # s
print('tank mass',8000*0.003*(2*np.pi*3.5*13 + 87.62))
T = np.arange(150,293.15,1)
_boiled_off_mass = []
for i in range(len(T)):
    _final_energy = _hydrogen.heating_to_temperature(temperature=T[i])  # internal energy at 20 °C
    _heating_energy = _final_energy.internal_energy - _initial_energy  # J/ 

    # calculate the boiled-off mass

    _total_heat = _total_heat_flux * _surface_area * _time  # J

    _boiled_off_mass.append(_total_heat / _heating_energy)  # kg

print('Boiled off mass:', _boiled_off_mass[-1])  # kg


plt.plot(T, _boiled_off_mass, label="Boiled-off Mass")
plt.xlabel("Temperature (K)")
plt.ylabel("Boiled-off Mass (kg)")
plt.title("Boiled-off Mass vs Temperature")
plt.legend()
plt.grid()
plt.show()