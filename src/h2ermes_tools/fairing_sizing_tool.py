import numpy as np
import matplotlib.pyplot as plt

from data import constants as cs
from data import material as mat
from h2ermes_tools.reentry.atmosphere import ExponentialAtmosphere

#spherical end cap

def volume_spherical_endcap_sheet(radius,thickness):
    surface_area = 2*np.pi*radius**2
    volume_sheet = surface_area*thickness
    return volume_sheet

def surface_area_spherical_end_cap(radius):
    return 2*np.pi*radius**2

def mass_spherical_endcap(material,radius,thickness):
    mass = material.rho*volume_spherical_endcap_sheet(radius,thickness)
    return mass

def chapman_stagnation_heat_flux(radius_nose,velocity_max,air_rho,n,m,c_star,v_c):
    qs = 1.63e-4*(air_rho/radius_nose)**0.5*velocity_max**3
    c1 = c_star * (1 / np.sqrt(air_rho)) * (1 / v_c ** 3)
    qc_max = c1 * (1 / radius_nose^n) * (air_rho)^(1-n) * (velocity_max)^m
    return qs,qc_max

def chapman_simplified_stagnation_heat_flux(radius_nose,velocity_max,air_rho):
    qs = 1.63e-4*(air_rho/radius_nose)**0.5*velocity_max**3
    return qs

def launcher_velocity(time):
    velocity = 0.09333*time**2 + 0.66667*time
    return velocity

def altitude_km(time):
    altitude_vec = 3.61/2*time**2 - 0.0285/3*time**3 + 1.83e-3/4*time**4 - 9.31e-6/5*time**5
    return altitude_vec/1e3

def heat_flux_internal(material, qs, time, thickness, dt):
    alpha = material.k / (material.rho * material.Cp)
    tau = thickness ** 2 / (alpha * np.pi ** 2)
    q_internal = [0] * len(time)

    for i in range(1, len(time)):
        q_internal[i] = q_internal[i-1] + (dt / tau) * (qs[i] - q_internal[i-1])
    return q_internal



def heat_flux_stagnation(time):
    altitude_vec=[]
    velocity_vec = []
    air_densities = []
    qs = []
    for i in range(len(time)):
        alt = altitude_km(time[i])
        altitude_vec.append(alt)
        velocity_vec.append(launcher_velocity(time[i]))
        atm = ExponentialAtmosphere(alt)
        air_densities.append(atm.rho_exp)  
        qs.append(chapman_simplified_stagnation_heat_flux(radius_nose,velocity_vec[i],air_densities[i]))
    print("Max stagnation heat flux is",max(qs),"W/m²")
    return qs

def total_heat_load(q_internal,dt,surface_area):
    heat_load = []
    total_heat_load = 0
    for i in range(len(q_internal)):
        heat_load_instantaneous = q_internal[i]*dt*surface_area
        total_heat_load = total_heat_load + heat_load_instantaneous
        heat_load.append(heat_load_instantaneous)
    return heat_load, total_heat_load

def heat_radiation(T,area,emissivity,sigma_boltzman):
    return sigma_boltzman*emissivity * area * T**4

def T_black_box(q_internal,absorptivity, dt, area,mass,Cp,T_ambient,emissivity,sigma_boltzman,heat_produced):
    q_absorbed = q_internal * absorptivity
    heat_load = total_heat_load(q_absorbed,dt,area)[1]

    heat_load = heat_load + heat_produced

    heat_in = 0
    T_comp = T_ambient

    while heat_in<heat_load:
        T_comp = T_comp +1
        dT = T_comp - T_ambient
        heat_conductive = (dT)*mass*Cp
        heat_rad = heat_radiation(T_comp,area,emissivity,sigma_boltzman)
        heat_in = heat_conductive+heat_rad
    return T_comp

def thickness_optimization(T_blackbox,T_max_operating,T_ambient,qs,time,dt,thickness_insulation,absorptivity_blackbox,surface_area_blackbox,mass_blackbox,Cp_blackbox,emissivity_blackbox,sigma_boltzman,q_in):
    while T_blackbox>T_max_operating:
        thickness_insulation = thickness_insulation+0.001
        q_internal = heat_flux_internal(material,qs,time,thickness_insulation,dt)
        T_blackbox = T_black_box(q_internal,absorptivity_blackbox,dt,surface_area_blackbox,mass_blackbox,Cp_blackbox,T_ambient,emissivity_blackbox,sigma_boltzman,q_in)

    print("New thickness of insulation is: ",thickness_insulation)
    print("Temperature of blackbox:",T_blackbox-273,"Celcius")
    return T_blackbox,thickness_insulation

def mass_insulation(rho_insulation,radius,thickness_insulation):
    surface_area = surface_area_spherical_end_cap(radius)
    mass = surface_area*thickness_insulation*rho_insulation
    return mass

def total_mass_heat_flux_calculation(material,rho_insulation,radius,thickness_insulation,thickness_wall,T_max_operating,T_ambient,time,dt,absorptivity_blackbox,surface_area_blackbox,mass_blackbox,Cp_blackbox,emissivity_blackbox,sigma_boltzman,T_surface_fuel_cell , area_fuel_cell,emissivity_fuel_cell,):
    qs = heat_flux_stagnation(time)
    q_internal = heat_flux_internal(material,qs,time,thickness_insulation,dt)
    q_in = heat_radiation(T_surface_fuel_cell , area_fuel_cell,emissivity_fuel_cell,sigma_boltzman)
    T_blackbox= T_black_box(q_internal,absorptivity_blackbox,dt,surface_area_blackbox,mass_blackbox,Cp_blackbox,T_ambient,emissivity_blackbox,sigma_boltzman,q_in)
    thickness_insulation = thickness_optimization(T_blackbox,T_max_operating,T_ambient,qs,time,dt,thickness_insulation,absorptivity_blackbox,surface_area_blackbox,mass_blackbox,Cp_blackbox,emissivity_blackbox,sigma_boltzman,q_in)[1]
    mass_insulation_endcap = mass_insulation(rho_insulation,radius,thickness_insulation)
    mass_endcap = mass_spherical_endcap(material,radius,thickness_wall)
    mass_total = mass_insulation_endcap + mass_endcap 
    return mass_total,q_internal,qs

def heatflux

if __name__ =="__main__":
    radius = 2.5
    thickness_wall = 0.002
    thickness_insulation = 0.01
    radius_nose = 0.5
    n_coefficient = 0.5
    m_coefficient = 3

    plot = True
    dt = 0.1 #s
    time_count = 170/dt
    time = np.linspace(0,170,int(time_count))
    air_densities = []
    qs = []
    altitude_vec = []
    velocity_vec = []
    T_ambient = 300 #K
    emissivity_blackbox = 0.5
    absorptivity_blackbox = 1
    sigma_boltzman = 5.67e-8
    T_surface_fuel_cell = 670
    area_fuel_cell = 2*0.4**2+4*0.4*1.2
    emissivity_fuel_cell = 0.7
    surface_area_blackbox =  6*0.4**2
    Cp_blackbox = 500 #J/kg*K
    T_max_operating = 273+50
    mass_blackbox = 5 #kg
    rho_insulation = 20 #kg/m3


    material = mat.Material(density = 7850,
                            youngs_modulus=200e9,
                            fracture_strength=800e6,
                            yield_strength=500e6,
                            thermal_conductivity = 0.1,
                            specific_heat = 500)


    mass_total,q_internal,qs = total_mass_heat_flux_calculation(material,rho_insulation,radius,thickness_insulation,thickness_wall,T_max_operating,T_ambient,time,dt,absorptivity_blackbox,surface_area_blackbox,mass_blackbox,Cp_blackbox,emissivity_blackbox,sigma_boltzman,T_surface_fuel_cell , area_fuel_cell,emissivity_fuel_cell,)

    print("Total mass is: ",mass_total,"kg")
    print("Max internal heat flux is:",max(q_internal),"W/m^2")
    if plot == True:
        plt.figure(figsize=(8, 5))
        plt.plot(time, q_internal, label="Internal Heat Flux")
        plt.plot(time, qs, label="Stagnation Heat Flux (qs)")
        plt.xlabel("Time (s)")
        plt.ylabel("Heat Flux (W/m²)")
        plt.title("Internal and Stagnation Heat Flux vs Time")
        plt.legend()
        plt.grid(True)
        plt.show()
    
