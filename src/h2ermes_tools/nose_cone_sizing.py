import numpy as np
import matplotlib.pyplot as plt

from data import constants as cs
from data import material as mat
from h2ermes_tools.reentry.atmosphere import ExponentialAtmosphere

def volume_spherical_endcap_sheet(radius,thickness):
    surface_area = 2*np.pi*radius**2
    volume_sheet = surface_area*thickness
    return volume_sheet

def surface_area_spherical_end_cap(radius):
    return 2*np.pi*radius**2

def mass_spherical_endcap(material,radius,thickness):
    mass = material.rho*volume_spherical_endcap_sheet(radius,thickness)
    return mass

def chapman_simplified_stagnation_heat_flux(radius_nose,velocity_max,air_rho):
    qs = 1.63e-4*(air_rho/radius_nose)**0.5*velocity_max**3
    return qs

def launcher_velocity(time):
    velocity = 0.09333*time**2 + 0.66667*time
    return velocity

def altitude_km(time):
    altitude_vec = 3.61/2*time**2 - 0.0285/3*time**3 + 1.83e-3/4*time**4 - 9.31e-6/5*time**5
    return altitude_vec/1e3

def heat_radiation(T,area,emissivity,sigma_boltzman):
    return sigma_boltzman*emissivity * area * T**4

def mass_mli(density_insulation,radius,thickness_insulation):
    surface_area = surface_area_spherical_end_cap(radius)
    mass = surface_area*density_insulation*thickness_insulation
    return mass

def projected_area_frontal_and_top(radius):
    proj_area_frontal = np.pi*radius**2/2
    proj_area_top = np.pi*radius**2
    return proj_area_frontal,proj_area_top

def heat_flux_orbital(solar_power, planetary_power, albedo_power, material_steel):
    # Geometry
    solar_power = 135311.68  # W
    planetary_power = 25795.63  # W
    albedo_power = 13604.74  # W

    incident_area = 7 * 9 + np.pi * 0.875 * 3.5
    planetary_flux = planetary_power / incident_area * material_steel.eps
    solar_flux = solar_power / incident_area * material_steel.abs
    albedo_flux = albedo_power / incident_area * material_steel.abs


    #total_heat_flux = (solar_flux + planetary_flux + albedo_flux)*projected_area_frontal_and_top(radius)[0]
    #total_heat_flux_out = sigma_boltzman*material_steel.ems*projected_area_frontal_and_top(radius)[0]*T_wall**4
    return planetary_flux+albedo_flux,solar_flux

def heat_flux_stagnation(time,radius):
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
        qs.append(chapman_simplified_stagnation_heat_flux(radius,velocity_vec[i],air_densities[i]))

    #force=  aerodynamic_load(air_densities,velocity_vec,radius,time)

    return qs

def radiative_temperature(total_heat_flux,time,material,sigma_boltzman):

    temp_wall = []

    for i in range(len(time)):
        temp_wall.append((total_heat_flux[i]/(material.eps*sigma_boltzman)+300**4)**0.25)

    return temp_wall

def temperature_blackbox_radiation(heat_tot,area_source,surface_area_blackbox,mass_blackbox,Cp_blackbox):
    heat_flux_source = heat_tot/area_source
    temp_wall = ((heat_flux_source*surface_area_blackbox)/(mass_blackbox*Cp_blackbox))
    return temp_wall

def fdm(thickness_wall,
        thickness_insulation,
        n_node_steel,
        n_node_mli,
        material_steel,
        dt,
        t_end,
        T_initial,
        time,
        radius,
        sigma_boltzman,
        emissivity_mli,
        Cp_insulation,
        density_insulation,
        thermal_conductivity_insulation):

    dx_steel = thickness_wall /(n_node_steel-1)
    k_combined = (thickness_wall+thickness_insulation)/(thickness_wall/material_steel.k+thickness_insulation/thermal_conductivity_insulation)

    alpha = k_combined/(material_steel.Cp*material_steel.rho)
    #alpha = material_steel.k/(material_steel.Cp*material_steel.rho)
    qs = heat_flux_stagnation(time,radius)

    T_wall_steel = radiative_temperature(qs,time,material_steel,sigma_boltzman)
    surface_area = surface_area_spherical_end_cap(radius)

    Fo_combined = alpha*dt/dx_steel**2
    T_steel = np.ones(n_node_steel)*T_initial

    T_inner_record_steel = []
    heat_radiation_propagation = []
    T_inner_record_mli = []
    length_time = len(time)
    T_steel[-1] = T_wall_steel[0]


    dx_mli = thickness_insulation/(n_node_mli-1)
    alpha_mli = alpha = thermal_conductivity_insulation/(Cp_insulation*density_insulation)
    Fo_mli = alpha_mli*dt/dx_mli**2
    T_mli = np.ones(n_node_mli)*T_initial

    if n_node_mli != 0:
        T_mli[-1] = T_steel[0]
        dx_mli = thickness_insulation/(n_node_mli-1)
        alpha_mli = alpha = thermal_conductivity_insulation/(Cp_insulation*density_insulation)
        Fo_mli = alpha_mli*dt/dx_mli**2
        T_mli = np.ones(n_node_mli)*T_initial


    if Fo_mli>0.5 and n_node_mli!=0:
        assert False, f"Fo mli is too high: Fo={Fo_mli} (must be < 0.5 for stability)"
    elif Fo_combined>0.5:
        assert False, f"Fo is too high: Fo={Fo_combined} (must be < 0.5 for stability)"

    for i in range(int(length_time)):
        T_new_steel = T_steel.copy()
        for j in range(1,n_node_steel-1):
            T_new_steel[j] = Fo_combined * (T_steel[j-1] + T_steel[j+1]) + (1 - 2*Fo_combined) * T_steel[j]

        if i <int(length_time)-1:
            T_new_steel[-1] = T_wall_steel[i+1]
        else:
            T_new_steel[-1] = T_wall_steel[i-1]

        T_new_steel[0] = T_new_steel[1]
        T_steel = T_new_steel
        T_inner_record_steel.append(T_new_steel[0])

        if n_node_mli !=0:
            T_new_mli = T_mli.copy()
            for j in range(1,n_node_mli-1):
                T_new_mli[j] = Fo_mli * (T_mli[j-1] + T_mli[j+1]) + (1 - 2*Fo_mli) * T_mli[j]

            T_new_mli[-1] = T_new_steel[0]


            T_new_mli[0] = T_new_mli[1]
            T_mli = T_new_mli
            T_inner_record_mli.append(T_new_mli[0])

        if n_node_mli ==0:
            heat_radiation_propagation.append(heat_radiation(T_new_steel[0],surface_area,emissivity_mli,sigma_boltzman))
        else:
            heat_radiation_propagation.append(heat_radiation(T_new_mli[0],surface_area,emissivity_mli,sigma_boltzman))
    
    return time,T_inner_record_steel,T_inner_record_mli, heat_radiation_propagation,max(T_inner_record_steel)

def total_mass(density_insulation,radius,thickness_insulation,thickness_wall,material_steel):
    mass_insulation = mass_mli(density_insulation,radius,thickness_insulation)
    mass_tank = surface_area_spherical_end_cap(radius)*thickness_wall*material_steel.rho
    return mass_insulation,mass_tank,

def thickness_optimization(thickness_wall,
        thickness_insulation,
        n_node_steel,
        n_node_mli,
        material_steel,
        dt,
        t_end,
        T_initial,
        time_mission,
        radius,
        sigma_boltzman,
        emissivity_mli,
        Cp_insulation,
        density_insulation,
        thermal_conductivity_insulation,
        safety_factor,
        max_operating_temperature,
        mass_per_area_insulation):
    max_dynamic_pressure = aerodynamic_load(radius,time)
    thickness_wall = critical_buckling_load(max_dynamic_pressure,material_steel,v=0.3,radius=2.5,phi = 80)

    max_temperature = fdm(thickness_wall,
        thickness_insulation,
        n_node_steel,
        n_node_mli,
        material_steel,
        dt,
        t_end,
        T_initial,
        time_mission,
        radius,
        sigma_boltzman,
        emissivity_mli,
        Cp_insulation,
        density_insulation,
        thermal_conductivity_insulation)[4]
    max_temperature_list= []
    thickness_insulation_list = []
    while max_temperature >max_operating_temperature:
        thickness_insulation = thickness_insulation + 0.0005

        max_temperature = fdm(thickness_wall,
        thickness_insulation,
        n_node_steel,
        n_node_mli,
        material_steel,
        dt,
        t_end,
        T_initial,
        time_mission,
        radius,
        sigma_boltzman,
        emissivity_mli,
        Cp_insulation,
        density_insulation,
        thermal_conductivity_insulation)[4]
        max_temperature_list.append(max_temperature)
        thickness_insulation_list.append(thickness_insulation)

        if thickness_insulation > 0.04:
            break

    thickness_insulation_list = [thick * 1000 for thick in thickness_insulation_list]
    max_temperature_list = [temp - 273.15 for temp in max_temperature_list]
    plt.figure(figsize=(8, 5))
    plt.plot(thickness_insulation_list, max_temperature_list, label="Max Temperature vs Nosecone Insulation thickness")
    plt.xlabel("Insulation Thickness [mm]")
    plt.ylabel("Max Temperature [Celcius]")
    plt.title("Max Temperature vs Nosecone Insulation Thickness")
    plt.legend()
    plt.grid(True)
    plt.show()


    thickness_insulation_final = thickness_insulation * safety_factor
    time_mission,T_inner_record_steel,T_inner_record_mli, heat_radiation_propagation,max_temperature = fdm(thickness_wall,
        thickness_insulation_final,
        n_node_steel,
        n_node_mli,
        material_steel,
        dt,
        t_end,
        T_initial,
        time_mission,
        radius,
        sigma_boltzman,
        emissivity_mli,
        Cp_insulation,
        density_insulation,
        thermal_conductivity_insulation)
    mass_mli,mass_tank = total_mass(density_insulation,radius,thickness_insulation_final,thickness_wall,material_steel)
    print("Total mass of spherical nosecone: ",mass_mli+mass_tank,"kg. Mass of insulation:",mass_mli,"The final thickness of wall is",thickness_wall,"and the thickness of the insulation is: ",thickness_insulation_final,"m.")
    mass_total = mass_mli+mass_tank
    return mass_total, T_inner_record_steel,T_inner_record_mli,heat_radiation_propagation

def aerodynamic_load(radius,time):
    surface_area = np.pi*radius**2
    Cd = 0.2
    dynamic_pressure = []
    force= []
    air_densities = []
    altitude_vec = []
    velocity_vec = []
    for i in range(len(time)):
        alt = altitude_km(time[i])
        altitude_vec.append(alt)
        velocity_vec.append(launcher_velocity(time[i]))
        atm = ExponentialAtmosphere(alt)
        air_densities.append(atm.rho_exp)  

        pressure_dynamic = (0.5*air_densities[i]*velocity_vec[i]**2)
        dynamic_pressure.append(pressure_dynamic)
        
        force.append(pressure_dynamic*surface_area*Cd)
    return max(dynamic_pressure)

def critical_buckling_load(max_dynamic_pressure,material,v,radius,phi):
    thickness = 0.00001
    pressure_record = []
    thickness_record = []

    pressure_critical = 2*material.E/((3*(1-v**2))**0.5)*(thickness/radius)**2
    lambda_factor =  (12*(1-v**2))**0.25 * (radius/thickness)**0.5 * 2*np.sin(np.radians(phi)/2)
    q_ebf = 0.693/((1-v)**0.2*lambda_factor**0.4)
    max_pressure = pressure_critical*q_ebf

    thickness_record.append(thickness)
    pressure_record.append(max_pressure)

    
    while max_pressure <max_dynamic_pressure:
        thickness += 0.00001
        pressure_critical = 2*material.E/((3*(1-v**2))**0.5)*(thickness/radius)**2
        lambda_factor =  (12*(1-v**2))**0.25 * (radius/thickness)**0.5 * 2*np.sin(np.radians(phi)/2)
        q_ebf = 0.693/((1-v)**0.2*lambda_factor**0.4)
        max_pressure = pressure_critical*q_ebf
        pressure_record.append(max_pressure)
        thickness_record.append(thickness)
    


    plt.figure(figsize=(8, 5))
    plt.plot(thickness_record, pressure_record, label="Critical Buckling Pressure vs Thickness")
    plt.xlabel("Thickness [m]")
    plt.ylabel("Critical Buckling Pressure [Pa]")
    plt.title("Critical Buckling Pressure vs Thickness")
    plt.legend()
    plt.grid(True)
    plt.show()
    return thickness


def total_flux_mission(dt,t_end,solar_power, planetary_power, albedo_power,radius, material_steel,cs,altitude_sc):
    time_ascent=np.linspace(0,t_end,int(t_end/dt))
    flux_stagnation = heat_flux_stagnation(time_ascent,radius)
    time_mission = np.linspace(0,24*3600,int(24*3600/dt))
    flux_orbit_constant,flux_orbit_solar = heat_flux_orbital(solar_power, planetary_power, albedo_power,material_steel)
    heat_flux_orbit_constant = np.ones(int(24*3600/dt-t_end/dt))*flux_orbit_constant
    total_heat_flux = []

    orbital_fraction_ellipse = np.arcsin(cs.earth_radius/(cs.earth_radius+altitude_sc))/(np.pi)

    orbital_period = (cs.gravitational_parameter/(cs.earth_radius+altitude_sc))**0.5

    heat_flux_orbit_varying = np.ones(int(orbital_period/dt))*flux_orbit_solar
    

    for i in range(int(orbital_fraction_ellipse*orbital_period/dt)):
        
        heat_flux_orbit_varying[i+int((1-orbital_fraction_ellipse)*orbital_period/dt)] = 0
    
    # Repeat the heat_flux_orbit_varying array to fill the 24-hour mission duration
    n_orbits = int(np.ceil((24*3600) / orbital_period))
    heat_flux_orbit_varying_full = np.tile(heat_flux_orbit_varying, n_orbits)
    # Trim to match the length of time_mission minus t_end/dt (as heat_flux_orbit_constant)
    heat_flux_orbit_varying_full = heat_flux_orbit_varying_full[:len(heat_flux_orbit_constant)]

    # Concatenate flux_stagnation (array of 1700 values) and heat_flux_orbit (array of ~86000 values)
    flux_stagnation = np.array(flux_stagnation)
    heat_flux_orbit_constant = np.array(heat_flux_orbit_constant)
    heat_flux_orbit_varying_full = np.array(heat_flux_orbit_varying_full)
    heat_flux_orbit =heat_flux_orbit_constant+ heat_flux_orbit_varying_full
    total_heat_flux = np.concatenate([flux_stagnation, heat_flux_orbit])

    return time_mission,total_heat_flux


if __name__ =="__main__":
    radius = 2.5
    thickness_wall = 0.002


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
    max_operating_temperature = 273+80
    mass_blackbox = 8000*surface_area_blackbox*0.003
    # Inputs for MLI
    margin_insulation = 1.05
    thickness_insulation = 0.00005
    mass_per_area_insulation = 0.429*margin_insulation #kg/m2
    heat_transfer_coefficient_insulation = 0.2 # W/m2/K
    thermal_conductivity_insulation = 4e-3#heat_transfer_coefficient_insulation*thickness_insulation*10

    emissivity_mli = 0.05
    Cp_insulation = 1
    density_insulation = mass_per_area_insulation/0.0004
    solar_power = 135311.68  # W
    planetary_power = 25795.63  # W
    albedo_power = 13604.74  # W

    # FDM inputs
    n_node_steel = 8
    n_node_mli = 0
    t_end = 170
    T_initial = 300
    safety_factor = 1.3
    heatflux = heat_flux_stagnation(time,radius)
    material_steel = mat.Material(density = 7850,
                            youngs_modulus=200e9,
                            fracture_strength=800e6,
                            yield_strength=500e6,
                            thermal_conductivity = 16,
                            specific_heat = 500,
                            emissivity=0.5,
                            absorptivity = 0.5)
    
    material_insulation = mat.Material(density = 7850,
                            thermal_conductivity = 3.7e-3,
                            emissivity = 0.5,
                            specific_heat = 500)
    #print(radiative_temperature([1180,230],radius,[1,1.5],material_steel,sigma_boltzman))
    time_mission,total_heat_flux = total_flux_mission(dt,t_end,solar_power,planetary_power,albedo_power,radius,material_steel,cs,600000)

    #mass_total,q_internal,qs = total_mass_heat_flux_calculation(material_steel,rho_insulation,radius,thickness_insulation,thickness_wall,T_max_operating,T_ambient,time,dt,absorptivity_blackbox,surface_area_blackbox,mass_blackbox,Cp_blackbox,emissivity_blackbox,sigma_boltzman,T_surface_fuel_cell , area_fuel_cell,emissivity_fuel_cell,)
    #T_internal = temperature_wall(thickness_insulation,material_insulation,qs)
    #print("Total mass is: ",mass_total,"kg")
    #print("Max internal heat flux is:",max(q_internal),"W/m^2")
    #T_internal, heat_flux_radiation = T_rise(time,thickness_insulation,material_insulation,radius,qs,sigma_boltzman)

    mass_total, T_inner_record_steel,T_inner_record_mli,heat_radiation_propagation = thickness_optimization(thickness_wall,
        thickness_insulation,
        n_node_steel,
        n_node_mli,
        material_steel,
        dt,
        t_end,
        T_initial,
        time,
        radius,
        sigma_boltzman,
        emissivity_mli,
        Cp_insulation,
        density_insulation,
        thermal_conductivity_insulation,
        safety_factor,
        max_operating_temperature,
        mass_per_area_insulation)
    

    plt.figure(figsize=(8, 5))
    plt.plot(time, T_inner_record_steel, label="FDM Inner Temperature (Steel)")
    #plt.plot(time, T_inner_record_mli, label="FDM Inner Temperature (MLI)")
    plt.xlabel("Time (s)")
    plt.ylabel("Temperature [K]")
    plt.title("FDM Inner Temperature vs Time")
    plt.legend()
    plt.grid(True)
    plt.show()

    # Plot heat_radiation_propagation against time
    plt.figure(figsize=(8, 5))
    plt.plot(time, heat_radiation_propagation, label="Heat Radiation Propagation")
    plt.xlabel("Time (s)")
    plt.ylabel("Heat Radiation [W]")
    plt.title("Heat Radiation Propagation vs Time")
    plt.legend()
    plt.grid(True)
    plt.show()

    


