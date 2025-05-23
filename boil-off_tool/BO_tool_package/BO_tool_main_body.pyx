import csv
import itertools as it
from math import *
from pathlib import Path
cimport numpy as cnp
import numpy as np
import pandas as pd
import Boiloff_tool as pf
from tqdm.auto import tqdm
import matplotlib.pyplot as plt
import CoolProp.CoolProp as CP
import os

cdef int months = pf.months

with open('mass_and_boil-off_v15-%sm.csv' %months, 'w') as csvfile:
    filewriter = csv.writer(csvfile, delimiter= ',', quotechar = '|', quoting=csv.QUOTE_MINIMAL)
    filewriter.writerow(['layer density MLI [layers/cm]', 'layers of MLI', 'MLI mass [kg]', 'cooler power [W]', 'Cooler mass [kg]', 't_shell [m]', 'Shell mass [kg]', 'Boil-off mass [kg]', 'Vent mass [kg]', 'Total mass [kg]', 'Boil-off rate [%/month]', 'Time start boil [days]', 'Time to boil-off [days]', 'T_liquid', 'T_vapour','Boil-off rate from start of month [%]'])
csvfile.close()

if months > 1:
    with open('mass_and_boil-off_v15-1m.csv', 'w') as csvfile:
        filewriter = csv.writer(csvfile, delimiter= ',', quotechar = '|', quoting=csv.QUOTE_MINIMAL)
        filewriter.writerow(['layer density MLI [layers/cm]', 'layers of MLI', 'MLI mass [kg]', 'cooler power [W]', 'Cooler mass [kg]', 't_shell [m]', 'Shell mass [kg]', 'Boil-off mass [kg]', 'Vent mass [kg]', 'Total mass [kg]', 'Boil-off rate [%/month]', 'Time start boil [days]', 'Time to boil-off [days]', 'T_liquid', 'T_vapour','Boil-off rate from start of month [%]'])
    csvfile.close()

if months > 2:
    with open('mass_and_boil-off_v15-2m.csv', 'w') as csvfile:
        filewriter = csv.writer(csvfile, delimiter= ',', quotechar = '|', quoting=csv.QUOTE_MINIMAL)
        filewriter.writerow(['layer density MLI [layers/cm]', 'layers of MLI', 'MLI mass [kg]', 'cooler power [W]', 'Cooler mass [kg]', 't_shell [m]', 'Shell mass [kg]', 'Boil-off mass [kg]', 'Vent mass [kg]', 'Total mass [kg]', 'Boil-off rate [%/month]', 'Time start boil [days]', 'Time to boil-off [days]', 'T_liquid', 'T_vapour','Boil-off rate from start of month [%]'])
    csvfile.close()

if months > 3:
    with open('mass_and_boil-off_v15-3m.csv', 'w') as csvfile:
        filewriter = csv.writer(csvfile, delimiter= ',', quotechar = '|', quoting=csv.QUOTE_MINIMAL)
        filewriter.writerow(['layer density MLI [layers/cm]', 'layers of MLI', 'MLI mass [kg]', 'cooler power [W]', 'Cooler mass [kg]', 't_shell [m]', 'Shell mass [kg]', 'Boil-off mass [kg]', 'Vent mass [kg]', 'Total mass [kg]', 'Boil-off rate [%/month]', 'Time start boil [days]', 'Time to boil-off [days]', 'T_liquid', 'T_vapour','Boil-off rate from start of month [%]'])
    csvfile.close()

if months > 4:
    with open('mass_and_boil-off_v15-4m.csv', 'w') as csvfile:
        filewriter = csv.writer(csvfile, delimiter= ',', quotechar = '|', quoting=csv.QUOTE_MINIMAL)
        filewriter.writerow(['layer density MLI [layers/cm]', 'layers of MLI', 'MLI mass [kg]', 'cooler power [W]', 'Cooler mass [kg]', 't_shell [m]', 'Shell mass [kg]', 'Boil-off mass [kg]', 'Vent mass [kg]', 'Total mass [kg]', 'Boil-off rate [%/month]', 'Time start boil [days]', 'Time to boil-off [days]', 'T_liquid', 'T_vapour','Boil-off rate from start of month [%]'])
    csvfile.close()

if months > 5:
    with open('mass_and_boil-off_v15-5m.csv', 'w') as csvfile:
        filewriter = csv.writer(csvfile, delimiter= ',', quotechar = '|', quoting=csv.QUOTE_MINIMAL)
        filewriter.writerow(['layer density MLI [layers/cm]', 'layers of MLI', 'MLI mass [kg]', 'cooler power [W]', 'Cooler mass [kg]', 't_shell [m]', 'Shell mass [kg]', 'Boil-off mass [kg]', 'Vent mass [kg]', 'Total mass [kg]', 'Boil-off rate [%/month]', 'Time start boil [days]', 'Time to boil-off [days]', 'T_liquid', 'T_vapour','Boil-off rate from start of month [%]']) 
    csvfile.close() 

if months > 6:
    with open('mass_and_boil-off_v15-6m.csv', 'w') as csvfile:
        filewriter = csv.writer(csvfile, delimiter= ',', quotechar = '|', quoting=csv.QUOTE_MINIMAL)
        filewriter.writerow(['layer density MLI [layers/cm]', 'layers of MLI', 'MLI mass [kg]', 'cooler power [W]', 'Cooler mass [kg]', 't_shell [m]', 'Shell mass [kg]', 'Boil-off mass [kg]', 'Vent mass [kg]', 'Total mass [kg]', 'Boil-off rate [%/month]', 'Time start boil [days]', 'Time to boil-off [days]', 'T_liquid', 'T_vapour','Boil-off rate from start of month [%]']) 
    csvfile.close() 

cdef conduction(k,DT,L):
    cdef dq = -k*DT/L
    return dq

cdef radiation(epsilon, T1, T2):
    cdef sigma = 5.670373e-8
    cdef dq = epsilon*sigma*-(T1**4.-T2**4.)
    return dq

cdef pressure_increase(T_s,L,m,V,Q):
    cdef R_G = 8.314469848  # J/K/mol (universal gas constant)
    cdef dpv = R_G*T_s/(L*m*V)*Q
    return dpv

cdef Rayleigh(L_c, g, beta, alpha, kinematic_viscosity, T_s, T_inf):
    cdef Ra = (g*beta)/(kinematic_viscosity*alpha)*abs(((T_s+T_inf)/2.-T_inf))*L_c**3.
    return Ra

cdef h(k,L_c,Ra):
    cdef h
    cdef Nu
    if Ra < 0:
        Ra = -Ra
        if Ra <= 1.e7:
            Nu = 0.642*Ra**(1./6.)
        elif 1.e7 < Ra <= 1.e10:
            Nu = 0.167*Ra**(1./4.)
        elif 1.e10 < Ra <= 5.*1.e13:
            Nu = 0.00053*Ra**(1./2.)
        else:
            Nu = 0.00053*Ra**(1./2.)
        h = -k/L_c*Nu
    else:
        if Ra <= 1.e7:
            Nu = 0.642*Ra**(1./6.)
        elif 1.e7 < Ra <= 1.e10:
            Nu = 0.167*Ra**(1./4.)
        elif 1.e10 < Ra <= 5.*1.e13:
            Nu = 0.00053*Ra**(1./2.)
        else:
            Nu = 0.00053*Ra**(1./2.)
        h = k/L_c*Nu

    return h

cdef h_vap(k,L_c,Ra):
    cdef h
    cdef Nu
    if Ra <= 1.e7:
        Nu = 4.5
    elif 1.e7 < Ra <= 1.e12:
        Nu = 0.08*Ra**(1./4.)
    else:
        Nu = 0.08*Ra**(1./4.)

    h = k/L_c*Nu
    return h


cdef NusseltfromRe (Reynolds,Prandtl):
    if Reynolds <= 1960:
        Nusselt = 3.66
    elif 1960 < Reynolds < 6420:
        Nusselt = 0.116*(Reynolds**(2./3.)-125.)*Prandtl**(0.4)
    else:
        Nusselt = 0.023*Reynolds**(0.8)*Prandtl**(0.4)
    return Nusselt

cdef Lockheed(T_H, T_C, emissivity, N_star, N_s):
    ### NON MODIFIED ###
    # cdef A = 7.30e-8
    # cdef n = 2.63
    # cdef T_m = (T_H+T_C)/2.
    # cdef B = 7.07e-10
    ### MODIFIED ### for perforated aluminized shields and Dacron net spacers
    cdef A = 2.4*10.**(-4)
    cdef n = 2.63
    cdef T_m = (0.017+7.*10.**(-6)*(800.-(T_H+T_C)/2.)+0.0228*log((T_H+T_C)/2.))
    cdef B = 4.944*10.**(-10)
    cdef C = 1.46e4
    cdef p_star = 1.33e-5 * 0.0075 #interstitial pressure in torr (multiply by 0.0075 to convert from Pa to torr)
    cdef m = -0.48
    cdef dq_solidconduction = A*N_star**n*T_m*-(T_H-T_C)/N_s
    cdef dq_radiation = B*emissivity*-(T_H**4.67-T_C**4.67)/N_s
    cdef dq_gasconduction = C*p_star*-(T_H**(m+1.)-T_C**(m+1.))/N_s
    cdef dq_total = dq_solidconduction+dq_radiation+dq_gasconduction
    return dq_total, dq_solidconduction, dq_radiation, dq_gasconduction

cdef Cryocooler_mass(Q_cooler, T_C, T_H, eta):
    cdef carnot = T_C/(T_H-T_C)
    cdef P_in_cooler = Q_cooler/eta*(1./carnot)
    cdef mass_cooler = 0.1422*P_in_cooler**(0.905)+0.325*P_in_cooler
    return mass_cooler

cdef efficiency(Q_cooler):
    cdef eta = 10. ** (min(
    (-1.26281 + 0.45936 * log10(Q_cooler) - 0.08743 * log10(Q_cooler) ** 2.),
    (-0.92237 + 0.07763 * log10(1 + Q_cooler)))) #FROM CHAI AND WILHITE
    return eta

class structure_sphere:
    def __init__(self, height, radius, thickness, sphere_parts, cylinder_parts, wall_density):
        self.height = height
        self.radius = radius
        self.thickness = thickness
        self.sphere_parts = float(sphere_parts)
        self.cylinder_parts = float(cylinder_parts)
        self.wall_density = wall_density
        self.wall_mass = (pi*((radius+thickness)**2.-radius**2.)*height+4./3.*pi*((radius+thickness)**3.-radius**3.))*wall_density

        # increment: data referring to the single tank section that is being analysed (corresponding to one node)
        #Calculations for endcaps (half-spheres)
        self.endcap_increment_area = 4.*pi*self.radius**2./self.sphere_parts
        self.endcap_internalvolume = (4/3.)*pi*(self.radius)**3.
        self.endcap_increment_wallvolume = (4 / 3.) * pi * ((self.radius + self.thickness) ** 3. - (self.radius) ** 3.)/self.sphere_parts
        self.endcap_increment_mass = self.endcap_increment_wallvolume * self.wall_density
        self.endcap_increment_length = 2.*pi*self.radius/self.sphere_parts

        #Calculations for cylinder
        self.increment_area = (2*pi*self.radius)*self.height/self.cylinder_parts
        self.internalvolume = (pi*(self.radius)**2.)*self.height
        self.increment_wallvolume = pi*((self.radius+self.thickness)**2. - (self.radius)**2.)*self.height/self.cylinder_parts
        self.increment_mass = self.increment_wallvolume * self.wall_density
        self.increment_length = height/self.cylinder_parts

class structure_ellipse:
    def __init__(self, height, radius, thickness, sphere_parts, cylinder_parts, wall_density, height_cap):
        self.height = height
        self.radius = radius
        self.thickness = thickness
        self.sphere_parts = float(sphere_parts)
        self.cylinder_parts = float(cylinder_parts)
        self.wall_density = wall_density

        self.wall_volume_cyl = pi*((radius+thickness)**2.-radius**2.)*height
        self.wall_volume_ellipse = 4./3.*pi*((height_cap + thickness)*(radius+thickness)**2.-(height_cap)*radius**2.)
        self.wall_mass = (self.wall_volume_cyl + self.wall_volume_ellipse)*wall_density 

        # increment: data referring to the single tank section that is being analysed (corresponding to one node)
        #Calculations for endcaps (oblate spheroids)
        self.ecc = 1 - height_cap**2./self.radius**2. #eccentricity of oblate spheroid
        self.endcap_increment_area = (2.*pi*self.radius**2. + pi*height_cap**2./self.ecc*log((1+self.ecc)/(1-self.ecc)))/self.sphere_parts #oblate spheroid
        self.endcap_internalvolume = (4/3.)*pi*(self.radius)**2.*height_cap
        self.endcap_increment_wallvolume = self.wall_volume_ellipse/self.sphere_parts
        self.endcap_increment_mass = self.endcap_increment_wallvolume * self.wall_density
        self.endcap_increment_length = 2.*pi*self.radius/self.sphere_parts # this is for conduction to opposite node, hence it is still a circle because the endcap is a spheroid

        #Calculations for cylinder
        self.increment_area = (2*pi*self.radius)*self.height/self.cylinder_parts
        self.internalvolume = (pi*(self.radius)**2.)*self.height
        self.increment_wallvolume = self.wall_volume_cyl/self.cylinder_parts
        self.increment_mass = self.increment_wallvolume * self.wall_density
        self.increment_length = height/self.cylinder_parts

# inputs are set only for initialization and will be changed later in the code:
cdef volume_calc(mass = 0., fill_ratio = 0., max_radius = 0., propellant_type = 'LH2', propellant_density = 71.541, tank_head_shape = 'sphere', ellipse_height = 0., volume = 0.):
    cdef double radius
    cdef double height

    if volume == 0.:
        volume = mass/propellant_density

    volume_tank = volume/fill_ratio
    #IF THERE IS A MAX RADIUS, THIS RADIUS IS ALSO USED. FOR THE HIGHER THE RADIUS, THE BIGGER VOLUME OVER AREA
    if tank_head_shape == 'sphere' or ellipse_height == max_radius:
        if max_radius != 0.:
            volume_sphere = (4/3.)*pi*max_radius**3.
            if volume_sphere > volume_tank:
                radius = (volume_tank/(4./3.*pi))**(1./3.)
                height = 0.
            else:
                volume_cylinder = volume_tank-volume_sphere
                height = volume_cylinder/(pi*max_radius**2.)
                radius = max_radius
        else:
            radius = (volume_tank/(4./3.*pi))**(1./3.)
            height = 0.

    elif tank_head_shape == 'ellipse':
        if max_radius != 0.:
            volume_ellipse = (4/3.)*pi*max_radius**2.*ellipse_height
            if volume_ellipse > volume_tank:
                radius = (volume_tank/(4./3.*pi*ellipse_height))**(1./2.)
                height = 0.
            else:
                volume_cylinder = volume_tank-volume_ellipse
                height = volume_cylinder/(pi*max_radius**2.)
                radius = max_radius
        else:
            radius = (volume_tank/(4./3.*pi))**(1./3.)
            height = 0. 
   
    return radius, height


################### Start of variables declaration #######################################################
# it is possible that not all the variables have been declared, however this is not a critical problem for the code, which will run normally
#variables declaration is mainly done to speed up the computations

# Constants
cdef double R_G = 8.314469848  # J/K/mol (universal gas constant)
cdef double grav = pf.gravity

# liquid properties
cdef density_liquid
cdef double L_liquid # latent heat of evaporation, J/g
cdef double L_liquid_mol # latent heat of evaporation, J/mol
cdef Cp_liquid
cdef k_liquid
cdef mu_liquid
cdef nu_liquid
cdef alpha_liquid
cdef beta_liquid
cdef double V_total
cdef fill_level
cdef double V_liquid
cdef m_liquid_initial

cdef Pr

# vapour properties 
cdef density_vapour
cdef Cp_vapour
cdef double Cp_prop_vapour
cdef Cp_mix
cdef k_vapour
cdef mu_He
cdef nu_vapour
cdef alpha_vapour
cdef beta_vapour
cdef double V_vapour
cdef m_vapour_initial
cdef double p_vapour_initial = pf.initial_pressure  #Pa

# Fluid properties
cdef Ra, h_liquid, L_c_liquid, L_c_vap, h_vapour, Ra_vap, Nu_vap
cdef mixture = pf.propellant_mixture
cdef m_molecular_propell, m_molecular_vapour
cdef fluid, pressurant
cdef double p_start = pf.initial_pressure
cdef mode = pf.press_mode
cdef m_mix
cdef double fill_level_initial = pf.fill_level_initial

# MLI properties
cdef layers
cdef list MLI_layers = pf.MLI_layers

cdef double density_MLI 
cdef double k_MLI = pf.k_MLI
cdef double Cp_MLI = pf.Cp_MLI
cdef t_MLI
cdef double emissivity_MLI = pf.emissivity_MLI
cdef list layer_densities = pf.layer_densities
cdef layer_density
cdef rho_layer = pf.layer_specific_weight
cdef rho_spacer = pf.spacer_specific_weight
cdef t_nom_layer = pf.layer_nom_thickness
cdef t_nom_spacer = pf.spacer_nom_thickness

cdef double emissivity = pf.emissivity_coating
cdef double absorptivity = pf.absorptivity_coating

# Cooler properties
cdef list cooler_power = pf.cooler_power
cdef double Q_cooler = 0.
cdef eta

# Variables for time loop
cdef double runtime = pf.t_run
cdef timestep = pf.dt
cdef t_boil = 0.
cdef t_vent = 0.
cdef dt, t, i
cdef int nr_nodes = pf.total_number_nodes 
cdef int n_node
cdef int n_subnode
cdef int nn
cdef cyl_sect 
cdef int subnode

# Masses
cdef m_shell, m_cooler, m_MLI
cdef m_vapour
cdef double m_liquid
cdef double m_boil_off
cdef double dm_boil_off

# Temperatures
cdef double T_start = pf.initial_temperature
cdef double T_space = pf.T_space
cdef double T_sat
cdef double T_liquid
cdef double T_vapour

cdef double dT_rest
cdef double dT_liquid
cdef double dT_vapour

cdef double T_sl_avg
cdef double T_sv_avg
cdef T_reject = pf.T_rejection

cdef double T_previous
cdef double T_previous_vap

# Tank geometry
cdef radius_max = pf.maximum_tank_radius
cdef tank_endcap = pf.tank_heads
cdef cap_height = pf.cap_height
cdef length_increm
cdef list subnodes_vect = pf.insulation_structure
cdef list subnodes_structure = pf.subnodes_structure
cdef int nr_subnodes
cdef tank_node = pf.include_tank_node
cdef double height_calc
cdef double radius_calc
cdef frac

# wall properties
cdef double k_wall = pf.k_wall 
cdef double Cp_wall = pf.Cp_wall
cdef double density_wall = pf.density_wall
cdef double t_wall = pf.t_wall

# Heat loads
cdef dP_sA_R, dP_As_R
cdef dP_sA_L, dP_As_L

cdef dP_A_next_sub_R, dP_A_opp_R, dP_i_prev_sub_R, dP_i_next_sub_R, dP_A_prev_n_R, dP_A_next_n_R
cdef dP_A_next_sub_L, dP_A_opp_L, dP_i_prev_sub_L, dP_i_next_sub_L, dP_A_prev_n_L, dP_A_next_n_L

cdef dP_Dl_R, dP_Dv_R
cdef dP_Dl_L, dP_Dv_L

cdef dP_lv, dP_cooler
cdef double dP_liquid_net
cdef P_rest

cdef dP_vl
cdef double dP_vapour_net

cdef double dP_system
cdef double dP_in
cdef double dP_out
cdef double dP_tot

# Other heat load sources inside tank
cdef Q_mixer
cdef dP_penetration
cdef dP_mixer

# Venting process
cdef venting_rate = 0 # just for initialization, it will be overwritten in the code with the calculated value
cdef m_rate = 1 # just for initialization, it will be overwritten in the code with the calculated value
cdef p_vent_target = pf.desired_venting_pressure
cdef p_vent 
cdef venting_started='FALSE'
cdef gamma, m_vapour_new
cdef vent_number
cdef p_max = pf.p_max
cdef double m_vent
cdef double dm_vent

# Vapor Cooled Shield system
cdef d_VCS = pf.diameter_VCS_tube
cdef mu_VCS, Cp_VCS, k_VCS, Re_VCS, Pr_VCS, Nu_VCS, h_VCS, A_VCS
cdef dP_VCS
cdef length_VCS = pf.length_VCS_tube
cdef VCS_on = pf.VCS_activity
cdef M_shield 
cdef t_shield = pf.shield_thickness
cdef Cp_shield = pf.shield_Cp
cdef density_shield = pf.shield_density
cdef double T_VCS = 0 # just for initialization, it will be calculated later in the code
cdef dt_VCS = pf.VCS_timestep
cdef VCS_position = pf.VCS_location

# Variable Density MLI (VDMLI)
cdef VD_MLI_layers = pf.VD_MLI_number_of_layers
cdef VD_MLI = pf.VD_MLI_activity
cdef layer
cdef single_venting = pf.venting_in_one_step
cdef Ns3, Ns2, Ns1

# Thermal environment data imported from ESATAN-TMS
cdef period, step
cdef radiative_case = pf.heat_fluxes_ESATAN
cdef Radiative_model = pf.Radiative_model_import

# Define vectors and matrices
cdef cnp.ndarray T_VD_MLI_L, T_VD_MLI_R, dP_VD_MLI_R, dP_VD_MLI_L
cdef cnp.ndarray layers_vector
cdef cnp.ndarray layer_density_vector
cdef cnp.ndarray density_MLI_vector

cdef cnp.ndarray T_matrix_R
cdef cnp.ndarray T_matrix_L
cdef cnp.ndarray dT_matrix_R
cdef cnp.ndarray dT_matrix_L
cdef cnp.ndarray Cp_vect

cdef cnp.ndarray Mass_sect_cyl
cdef cnp.ndarray Mass_sect_sph

cdef cnp.ndarray dP_net_matrix_R
cdef cnp.ndarray dP_net_matrix_L
cdef cnp.ndarray dP_liquid_matrix
cdef cnp.ndarray dP_vapour_matrix

# Vectors for plots
cdef save_vectors = pf.save_vectors_for_plots
cdef cnp.ndarray dm_vent_vector
cdef cnp.ndarray m_boiloff_vector
cdef cnp.ndarray p_vapour_vector
cdef cnp.ndarray m_vapour_vector
cdef cnp.ndarray time_vect
cdef cnp.ndarray T_vapour_vect
cdef cnp.ndarray T_VCS_vect
cdef cnp.ndarray T_liquid_vect
cdef cnp.ndarray dP_vapour_vector
cdef cnp.ndarray dP_liquid_vector
cdef cnp.ndarray dP_radiative_net, dP_emitted_back 
cdef cnp.ndarray dP_space

# Parameters for data storage and output generation
cdef double Boiloff_from_start_of_month
cdef m_boil_off_previous_month, m_liquid_initial_previous_month
cdef double Boil_off_monthly

cdef option_nr = 0
cdef filename
cdef BASE_PATH = Path.cwd()/ 'Output_data'
cdef option

cdef bint one_month = False
cdef bint two_months = False
cdef bint three_months = False
cdef bint four_months = False
cdef bint five_months = False
cdef bint six_months = False

################################# End of variables declaration #############################################

################## INITIATE PROGRAM ###################################
print("\n initiate program...")

cdef list options = list(it.product(layer_densities, MLI_layers, cooler_power))
options = list(it.chain(options, pf.custom_designs))

print(options)

##################### Define radiative case ################
with open(Radiative_model,'rb') as f:

    RadData = np.load(f,allow_pickle = 'True')
    abs_flux_tot_L = RadData['abs_flux_tot_L'] # total absorbed fluxes R side
    abs_flux_tot_R = RadData['abs_flux_tot_R'] # total absorbed fluxes L side
    inc_flux_ISIA_R = RadData['inc_flux_ISIA_R'] # solar + abedo incoming fluxes R side
    inc_flux_ISIA_L = RadData['inc_flux_ISIA_L'] # solar + abedo incoming fluxes L side
    inc_flux_IP_R = RadData['inc_flux_IP_R'] # planetary incoming fluxes R side
    inc_flux_IP_L = RadData['inc_flux_IP_L'] # planetary incoming fluxes L side
    time_steps = RadData['time_steps']

time_steps_original = time_steps #save original time steps vector
period = time_steps_original[-1]
#print('Total absorbed fluxes right side')
#print(abs_flux_tot_R) 

#print('Total absorbed fluxes left side')
#print(abs_flux_tot_L)

############################################################

for option in tqdm(options):
    print("\n New combination")

    option_nr += 1

    t_boil = 0.
    t_vent = 0.

    one_month = False
    two_months = False
    three_months = False
    four_months = False
    five_months = False
    six_months = False

    layer_density, layers, Q_cooler = option

    if mixture == 'LH2':
        ### Hydrogen properties ###
        print('Propellant is Liquid Hydrogen')
        fluid = 'Hydrogen'
    elif mixture == 'CH4':
        ### Methane properties ###
        print('Propellant is Liquid Methane')
        fluid = 'Methane'

    if mode=='Autogenous':
        T_start =  CP.PropsSI('T','P',p_start,'Q',0,fluid) #K
        T_sat = T_start
        density_liquid = CP.PropsSI('D','P',p_start,'Q',0,fluid) #kg/m3
        Cp_liquid = CP.PropsSI('CP0MASS','P',p_start,'Q',0,fluid) #J/kg*K
        k_liquid = CP.PropsSI('CONDUCTIVITY','P',p_start,'Q',0,fluid) #W/m*K
        mu_liquid = CP.PropsSI('VISCOSITY','P',p_start,'Q',0,fluid) #Pa*s
        nu_liquid = mu_liquid/density_liquid 
        alpha_liquid = k_liquid/density_liquid/Cp_liquid #m2/s
        beta_liquid = CP.PropsSI('ISOBARIC_EXPANSION_COEFFICIENT','P',p_start,'Q',0,fluid) #1/K
        m_molecular_propell = CP.PropsSI('M','P',p_start,'Q',0,fluid)*1e3 #g/mol
        L_liquid = (CP.PropsSI('H','P',p_start,'Q',1,fluid) - CP.PropsSI('H','P',p_start,'Q',0,fluid))*1e-3  # Latent heat of vaporization J/g
        L_liquid_mol = L_liquid*m_molecular_propell #J/mol

        ### Autogenous pressurant properties ###
        pressurant = 'Hydrogen'
        density_vapour = CP.PropsSI('D','P',p_start,'Q',1,pressurant) #kg/m3
        Cp_vapour = CP.PropsSI('CP0MASS','P',p_start,'Q',1,pressurant) #J/kg*K
        k_vapour = CP.PropsSI('CONDUCTIVITY','P',p_start,'Q',1,pressurant) #W/m*K
        mu_He = CP.PropsSI('VISCOSITY','P',p_start,'Q',1,pressurant) #Pa*s
        nu_vapour = mu_He/density_vapour
        alpha_vapour = k_vapour/density_vapour/Cp_vapour #m2/s
        beta_vapour = CP.PropsSI('ISOBARIC_EXPANSION_COEFFICIENT','P',p_start,'Q',1,pressurant) #1/K
        m_molecular_vapour= CP.PropsSI('M','P',p_start,'Q',1,pressurant)*1e3 #g/mol

    elif mode == 'HeliumPress':

        density_liquid = CP.PropsSI('D','P',p_start,'T',T_start,fluid) #kg/m3
        Cp_liquid = CP.PropsSI('CP0MASS','P',p_start,'T',T_start,fluid) #J/kg*K
        k_liquid = CP.PropsSI('CONDUCTIVITY','P',p_start,'T',T_start,fluid) #W/m*K
        mu_liquid = CP.PropsSI('VISCOSITY','P',p_start,'T',T_start,fluid) #Pa*s
        nu_liquid = mu_liquid/density_liquid 
        alpha_liquid = k_liquid/density_liquid/Cp_liquid #m2/s
        beta_liquid = CP.PropsSI('ISOBARIC_EXPANSION_COEFFICIENT','P',p_start,'T',T_start,fluid) #1/K
        m_molecular_propell = CP.PropsSI('M','P',p_start,'T',T_start,fluid)*1e3 #g/mol
        L_liquid = (CP.PropsSI('H','P',p_start,'Q',1,fluid) - CP.PropsSI('H','P',p_start,'Q',0,fluid))*1e-3  # Latent heat of vaporization J/g
        L_liquid_mol = L_liquid*m_molecular_propell #J/mol

        ### Helium properties ###
        pressurant = 'Helium'
        density_vapour = CP.PropsSI('D','P',p_start,'T',T_start,pressurant) #kg/m3
        Cp_vapour = CP.PropsSI('CP0MASS','P',p_start,'T',T_start,pressurant) #J/kg*K
        k_vapour = CP.PropsSI('CONDUCTIVITY','P',p_start,'T',T_start,pressurant) #W/m*K
        mu_He = CP.PropsSI('VISCOSITY','P',p_start,'T',T_start,pressurant) #Pa*s
        nu_vapour = mu_He/density_vapour
        alpha_vapour = k_vapour/density_vapour/Cp_vapour #m2/s
        beta_vapour = CP.PropsSI('ISOBARIC_EXPANSION_COEFFICIENT','P',p_start,'T',T_start,pressurant) #1/K
        m_molecular_vapour= CP.PropsSI('M','P',p_start,'T',T_start,pressurant)*1e3 #g/mol

    ### Structure and Masses (MLI, Shell, Cooler)###

    radius_calc, height_calc = volume_calc(mass = pf.propellant_mass, fill_ratio = fill_level_initial, max_radius = radius_max, propellant_type = pf.propellant_mixture, propellant_density = density_liquid, tank_head_shape = tank_endcap, ellipse_height = cap_height)
    print('The propellant tank radius is %s [m] and the height is %s [m]' %(radius_calc,round(height_calc, 1)))

    # compute MLI specific weight (kg/m2/layer) using relation from Chai et al.
    density_MLI = rho_layer + rho_spacer/t_nom_spacer*(1./(100.*layer_density)-t_nom_layer)

    if tank_endcap == 'sphere' or cap_height == radius_max:

        prop = structure_sphere(height = height_calc, radius= radius_calc, thickness=t_wall, sphere_parts=pf.sphere_parts_tank, cylinder_parts=pf.cylinder_parts_tank, wall_density=density_wall)
        m_MLI_sph = (4.*pi*(prop.radius+prop.thickness)**2.) * density_MLI*layers

    elif tank_endcap == 'ellipse':

        prop = structure_ellipse(height = height_calc, radius= radius_calc, thickness=t_wall, sphere_parts=pf.sphere_parts_tank, cylinder_parts=pf.cylinder_parts_tank, wall_density=density_wall, height_cap = cap_height)
        m_MLI_sph = (2.*pi*(prop.radius+prop.thickness)**2. + pi*(cap_height+prop.thickness)**2./prop.ecc*log((1+prop.ecc)/(1-prop.ecc)))*density_MLI*layers
  
    # cylinder calculations are the same for both spherical and ellpisoidal endcaps, thus they can go outside the loop
    m_MLI_cyl = (2.*pi*(prop.radius+prop.thickness)*prop.height)*density_MLI*layers
    m_MLI = m_MLI_cyl + m_MLI_sph

    # Masses of sectors (nodes) to be analysed 
    m_sect_MLI_sph = m_MLI_sph/prop.sphere_parts/2. #divide by 2 because nodes for MLI in code are 2 (left L and right R)
    m_sect_MLI_cyl = m_MLI_cyl/prop.cylinder_parts/2.

    # VCS shield mass
    if VCS_on == 'TRUE':
        if tank_endcap == 'sphere' or cap_height == radius_max:
            M_shield = ((pi * ((prop.radius + prop.thickness + t_shield) ** 2. - (prop.radius + prop.thickness ) ** 2.) * prop.height) \
                + 4 / 3. * pi * ((prop.radius + prop.thickness + t_shield) ** 3. - (prop.radius + prop.thickness ) ** 3.)) * density_shield

        elif tank_endcap == 'ellipse':
            M_shield = ((pi * ((prop.radius + prop.thickness + t_shield) ** 2. - (prop.radius + prop.thickness) ** 2.) * prop.height) \
                + 4 / 3. * pi * ((cap_height+prop.thickness + t_shield)*(prop.radius + prop.thickness + t_shield) ** 2. - (cap_height+prop.thickness)*(prop.radius + prop.thickness) ** 2.)) * density_shield
        
        print('VCS shield mass is:', M_shield)
   
    # Cooler and shell mass
    if Q_cooler == 0.:
        eta = 0.
        m_cooler = 0.
    else:
        eta = efficiency(Q_cooler)
        m_cooler = Cryocooler_mass(Q_cooler, 20., T_reject, eta) #T_H = 263, which should be the spacecraft temperature according to Chai 2014.

    m_shell = prop.wall_mass

    # Vapour and liquid volume and initial mass

    V_total = prop.internalvolume+prop.endcap_internalvolume
    V_liquid = fill_level_initial*V_total
    V_vapour = V_total-V_liquid
    m_liquid_initial = density_liquid*V_liquid
    m_vapour_initial = p_start*V_vapour*m_molecular_vapour*1e-3/(R_G*T_start)  
    print('Initial vapour mass is:',m_vapour_initial)

    # heat leak from propellant mixer 
    Q_mixer = V_total/16.8*1.
    
    ###### VDMLI calculations (when applicable) #########
    if VD_MLI == 'TRUE':
        
        # set vectors for layers density and number of layers in MLI sector
        layer_density_vector = np.array([layer_density, 2.*layer_density/3., layer_density/3. ])
        Ns3 = ceil(layers/(100.*layer_density_vector[2])/(1./(100.*layer_density_vector[0]) + 1./(100.*layer_density_vector[2]) + layer_density_vector[1]/(100.*layer_density_vector[0]*layer_density_vector[2])))
        Ns2 = ceil(layer_density_vector[1]*Ns3/layer_density_vector[0])
        Ns1 = layers - Ns3 - Ns2
        layers_vector = np.array([Ns3, Ns2, Ns1])

        print('VDMLI # of layers:', layers_vector, 'and densities vector:', layer_density_vector)

        t_MLI = np.sum(layers_vector*1./(layer_density_vector*100.))
        t_MLI_vector = layers_vector*1./(layer_density_vector*100.)
        density_MLI_vector = rho_layer + rho_spacer/t_nom_spacer*(1./(100.*layer_density_vector)-t_nom_layer)

        Mass_sect_cyl_VDMLI = (2.*pi*(prop.radius+prop.thickness)*prop.height/prop.cylinder_parts)*layers_vector*density_MLI_vector 

        if tank_endcap == 'sphere' or cap_height == radius_max:

            Mass_sect_sph_VDMLI = (4.*pi * (prop.radius + prop.thickness) ** 2./prop.sphere_parts) *layers_vector * density_MLI_vector 
            m_MLI_sph = np.sum( (4.*pi*(prop.radius+prop.thickness)**2.)*density_MLI_vector*layers_vector)

        elif tank_endcap == 'ellipse':
            
            Mass_sect_sph_VDMLI = (2.*pi*(prop.radius+prop.thickness)**2. + pi*(cap_height+prop.thickness)**2./prop.ecc*log((1+prop.ecc)/(1-prop.ecc)))/prop.sphere_parts*density_MLI_vector*layers_vector
            m_MLI_sph = np.sum( (2.*pi*(prop.radius+prop.thickness)**2. + pi*(cap_height+prop.thickness)**2./prop.ecc*log((1+prop.ecc)/(1-prop.ecc)))*density_MLI_vector*layers_vector)

        m_MLI = np.sum(density_MLI_vector*layers_vector*(2.*pi*(prop.radius+prop.thickness)*prop.height))+ m_MLI_sph #recompute MLI mass   

    else:
        t_MLI = layers*1./(layer_density*100.)

    print('MLI mass is:', m_MLI)

    ###### TIME ######
    dt = 1. # 0.1
    t = 0.
    vent_number = 0 # to keep track of number of venting processes

    fill_level = fill_level_initial
    p_vapour = p_vapour_initial
    p_vent = p_vent_target

    ### Get number of subnodes for matrices initialization
    nr_subnodes = int(len(subnodes_structure))
    print('Number of subnodes is:',nr_subnodes)

    ###### INITIAL TEMPERATURES ######
    T_matrix_R = np.ones((nr_nodes//2,nr_subnodes))*T_start
    T_matrix_L = np.ones((nr_nodes//2,nr_subnodes))*T_start
    T_VD_MLI_R = np.ones((nr_nodes//2,len(VD_MLI_layers)))*T_start 
    T_VD_MLI_L = np.ones((nr_nodes//2,len(VD_MLI_layers)))*T_start 
    T_liquid = T_start
    T_previous = T_start
    T_vapour = T_start
    T_previous_vap = T_start

    if abs(m_liquid_initial/pf.propellant_mass-1.) > 0.001:
        print('initial propellant mass differs %s percent from target' %((m_liquid_initial/pf.propellant_mass-1.)*100.))

    # Masses initialization
    m_vapour = m_vapour_initial
    m_liquid = m_liquid_initial
    m_vent = 0.
    dm_vent = 0.
    m_boil_off = 0.
    dm_boil_off = 0.

    ####### Vectors (and matrices) initialization
    dT_matrix_R = np.ones((nr_nodes//2,nr_subnodes)) 
    dT_matrix_L = np.ones((nr_nodes//2,nr_subnodes))

    if tank_node == 'True':
        Cp_vect = np.array([Cp_MLI, Cp_MLI, Cp_wall])
        m_sect_wall_sph = prop.endcap_increment_mass/2.
        m_sect_wall_cyl = prop.increment_mass/2.
        Mass_sect_cyl = np.array([m_sect_MLI_cyl, m_sect_MLI_cyl, m_sect_wall_cyl])
        Mass_sect_sph = np.array([m_sect_MLI_sph, m_sect_MLI_sph, m_sect_wall_sph])
        if VD_MLI == 'True':
            Mass_sect_cyl_VDMLI = np.append(Mass_sect_cyl_VDMLI,m_sect_wall_cyl)
            Mass_sect_cyl_VDMLI = np.append(Mass_sect_cyl_VDMLI,m_sect_wall_sph)
    else:
        Cp_vect = np.array([Cp_MLI, Cp_MLI])
        Mass_sect_cyl = np.array([m_sect_MLI_cyl, m_sect_MLI_cyl])
        Mass_sect_sph = np.array([m_sect_MLI_sph, m_sect_MLI_sph])

    dP_net_matrix_R = np.zeros((nr_nodes//2,nr_subnodes))
    dP_net_matrix_L = np.zeros((nr_nodes//2,nr_subnodes))
    dP_liquid_matrix = np.zeros((nr_nodes//2,2))
    dP_vapour_matrix = np.zeros((nr_nodes//2,2))

    dP_VD_MLI_R = np.zeros((nr_nodes//2,len(VD_MLI_layers)))
    dP_VD_MLI_L = np.zeros((nr_nodes//2,len(VD_MLI_layers)))
    
    # initialise vectors for plots (nox exact allocations but based on runtime value*10), they can be set based on the user's needs
    if save_vectors == 'TRUE':
        time_vect = np.zeros(int(10*runtime))

        dm_vent_vector = np.zeros(int(10*runtime))
        m_boiloff_vector = np.zeros(int(10*runtime))

        p_vapour_vector = np.ones(int(10*runtime))*p_vapour
        m_vapour_vector = np.zeros(int(10*runtime))
    
        T_vapour_vect = np.ones(int(10*runtime))*T_start
        T_liquid_vect = np.ones(int(10*runtime))*T_start
        T_VCS_vect = np.zeros(int(10*runtime))

        dP_vapour_vector = np.zeros(int(10*runtime))
        dP_liquid_vector = np.zeros(int(10*runtime))

        dP_radiative_net = np.zeros(int(10*runtime))
        dP_emitted_back = np.zeros(int(10*runtime))
    
    #dP_node_net = np.zeros((nr_nodes//2,2))
    dP_space = np.zeros((nr_nodes//2,2))
    dP_sA = np.zeros((nr_nodes//2,2))
    dP_As = np.zeros((nr_nodes//2,2))
    ##########################################################################


    ############## Start time loop analysis ###################################
    nn = 1
    step = 0
    time_steps = time_steps_original #reset time steps vector for variable fluxes along orbit

    while t <= runtime:

        if t > 60.*60.*24.*14.:  #to change (increase) timestep once 2 weeks of storage passed
            dt = timestep

        t = t+dt

        # variable fluxes based on orbit data from ESATAN
        if step == (len(time_steps)-1) and t > time_steps[step]:
            step = 0
            time_steps = time_steps + period
        elif t >= time_steps[step+1]:
            step +=1

        if (T_liquid-T_previous) >= 0.5:
            print('\n change liquid properties')
            
            if mode == 'HeliumPress':
                density_liquid = CP.PropsSI('D','P',p_start,'T',T_liquid,fluid) #kg/m3
                Cp_liquid = CP.PropsSI('CP0MASS','P',p_start,'T',T_liquid,fluid) #J/kg*K
                k_liquid = CP.PropsSI('CONDUCTIVITY','P',p_start,'T',T_liquid,fluid) #W/m*K
                mu_liquid = CP.PropsSI('VISCOSITY','P',p_start,'T',T_liquid,fluid) #Pa*s
                nu_liquid = mu_liquid/density_liquid 
                alpha_liquid = k_liquid/density_liquid/Cp_liquid
                beta_liquid = CP.PropsSI('ISOBARIC_EXPANSION_COEFFICIENT','P',p_start,'T',T_liquid,fluid) #1/K
                T_previous = T_liquid
            elif mode == 'Autogenous':
                density_liquid = CP.PropsSI('D','T',T_liquid,'Q',0,fluid) #kg/m3
                Cp_liquid = CP.PropsSI('CP0MASS','T',T_liquid,'Q',0,fluid) #J/kg*K
                k_liquid = CP.PropsSI('CONDUCTIVITY','T',T_liquid,'Q',0,fluid) #W/m*K
                mu_liquid = CP.PropsSI('VISCOSITY','T',T_liquid,'Q',0,fluid) #Pa*s
                nu_liquid = mu_liquid/density_liquid 
                alpha_liquid = k_liquid/density_liquid/Cp_liquid
                beta_liquid = CP.PropsSI('ISOBARIC_EXPANSION_COEFFICIENT','T',T_liquid,'Q',0,fluid) #1/K
                T_previous = T_liquid


        if (T_vapour-T_previous_vap) >= 1.0:
            print('change vapour properties')

            if mode == 'HeliumPress':
                density_vapour = CP.PropsSI('D','P',p_start,'T',T_vapour,pressurant) #kg/m3
                Cp_vapour = CP.PropsSI('CP0MASS','P',p_start,'T',T_vapour,pressurant) #J/kg*K
                k_vapour = CP.PropsSI('CONDUCTIVITY','P',p_start,'T',T_vapour,pressurant) #W/m*K
                mu_He = CP.PropsSI('VISCOSITY','P',p_start,'T',T_vapour,pressurant) #Pa*s
                nu_vapour = mu_He/density_vapour
                alpha_vapour = k_vapour/density_vapour/Cp_vapour
                beta_vapour = CP.PropsSI('ISOBARIC_EXPANSION_COEFFICIENT','P',p_start,'T',T_vapour,pressurant) #1/K
                T_previous_vap = T_vapour
            elif mode == 'Autogenous':
                density_vapour = CP.PropsSI('D','P',p_vapour,'T',T_vapour,pressurant) #kg/m3
                Cp_vapour = CP.PropsSI('CP0MASS','P',p_vapour,'T',T_vapour,pressurant) #J/kg*K
                k_vapour = CP.PropsSI('CONDUCTIVITY','P',p_vapour,'T',T_vapour,pressurant) #W/m*K
                mu_He = CP.PropsSI('VISCOSITY','P',p_vapour,'T',T_vapour,pressurant) #Pa*s
                nu_vapour = mu_He/density_vapour
                alpha_vapour = k_vapour/density_vapour/Cp_vapour
                beta_vapour = CP.PropsSI('ISOBARIC_EXPANSION_COEFFICIENT','P',p_vapour,'T',T_vapour,pressurant) #1/K
                T_previous_vap = T_vapour

        L_c_liquid = fill_level*(2.*prop.radius+prop.height)
        L_c_vap = (1-fill_level)*(2.*prop.radius+prop.height)

        if nr_nodes == 4: # for spherical tank
            T_sl_avg = (np.sum(T_matrix_R[:,-1]) + np.sum(T_matrix_L[:,-1]))
            T_sv_avg = (np.sum(T_matrix_R[:,-1]) + np.sum(T_matrix_L[:,-1]))
        else:
            if fill_level > (0.5*prop.endcap_internalvolume+prop.internalvolume)/V_total:
                T_sl_avg = (np.sum(T_matrix_R[:,-1]) + np.sum(T_matrix_L[:,-1]))/12.
                T_sv_avg = (T_matrix_R[0,-1] + T_matrix_L[0,-1])/2.
            elif fill_level > (0.5*prop.endcap_internalvolume+0.75*prop.internalvolume)/V_total:
                T_sl_avg = (np.sum(T_matrix_R[1:,-1]) + np.sum(T_matrix_L[1:,-1]))/10.
                T_sv_avg = (np.sum(T_matrix_R[0:2,-1]) + np.sum(T_matrix_L[0:2,-1]))/4.
            elif fill_level > (0.5*prop.endcap_internalvolume+0.5*prop.internalvolume)/V_total:
                T_sl_avg = (np.sum(T_matrix_R[2:,-1]) + np.sum(T_matrix_L[2:,-1]))/8.
                T_sv_avg = (np.sum(T_matrix_R[0:3,-1]) + np.sum(T_matrix_L[0:3,-1]))/6.
            elif fill_level > (0.5*prop.endcap_internalvolume+0.25*prop.internalvolume)/V_total:
                T_sl_avg = (np.sum(T_matrix_R[3:,-1]) + np.sum(T_matrix_L[3:,-1]))/6.
                T_sv_avg = (np.sum(T_matrix_R[0:4,-1]) + np.sum(T_matrix_L[0:4,-1]))/8.
            elif fill_level > (0.5*prop.endcap_internalvolume)/V_total:
                T_sl_avg = (np.sum(T_matrix_R[4:,-1]) + np.sum(T_matrix_L[4:,-1]))/4.
                T_sv_avg = (np.sum(T_matrix_R[0:5,-1]) + np.sum(T_matrix_L[0:5,-1]))/8.
            else:
                T_sl_avg = (T_matrix_R[5,-1] + T_matrix_L[5,-1])/2.
                T_sv_avg = (np.sum(T_matrix_R[:,-1]) + np.sum(T_matrix_L[:,-1]))/12.

        Ra = Rayleigh(L_c_liquid, grav, beta_liquid, alpha_liquid, nu_liquid, T_sl_avg, T_liquid)
        h_liquid = h(k_liquid, L_c_liquid, Ra)

        Ra_vap = Rayleigh(L_c_vap, grav, beta_vapour, alpha_vapour, nu_vapour, T_sv_avg, T_vapour)
        h_vapour = h_vap(k_vapour, L_c_vap, Ra_vap)

        ########### Iteration for nodal analysis starts here
        cyl_sect = pf.cylinder_parts_tank/2. # Variable for counting the current cylindrical section that is being analysed, can also be written as (nr_nodes - nr_spherical nodes )/2

        n_node = 1
        while n_node <= nr_nodes//2:
            n_subnode = 0
            subnode = 0
            if subnodes_vect [subnode] == "MLI":  
                # NODE A
                if n_node == 1 or n_node == nr_nodes//2: # this is the top (or bottom) of the tank, there is no previous node (or next) 
                # NODE A top or bottom of tank (spherical part)
                    # Previous SUBnode (space), when using ESATAN model, true area (and not frontal) is used

                    # If absorbed fluxes from ESATAN are used
                    if radiative_case == 'Absorbed_Fluxes':
                        dP_sA_R = abs_flux_tot_R[step,n_node-1]*prop.endcap_increment_area #take first time step for now, then implement in code for all steps
                        dP_sA_L = abs_flux_tot_L[step,n_node-1]*prop.endcap_increment_area
 
                    # If direct fluxes from ESATAN are used
                    elif radiative_case == 'Direct_Fluxes':
                        dP_sA_R = absorptivity*inc_flux_ISIA_R[step,n_node-1]*prop.endcap_increment_area \
                                + emissivity*inc_flux_IP_R[step,n_node-1]*prop.endcap_increment_area

                        dP_sA_L = absorptivity*inc_flux_ISIA_L[step,n_node-1]*prop.endcap_increment_area \
                                + emissivity*inc_flux_IP_L[step,n_node-1]*prop.endcap_increment_area

                    # heat emitted back by spacecraft 
                    dP_As_R = radiation(emissivity,T_matrix_R[n_node-1,n_subnode],T_space)*prop.endcap_increment_area
                    dP_As_L = radiation(emissivity,T_matrix_L[n_node-1,n_subnode],T_space)*prop.endcap_increment_area
                    
                    # Next SUBnode
                    if VD_MLI == 'TRUE':
                        
                        for layer in VD_MLI_layers:

                            if layer == 0: 
                                dP_A_next_sub_R = Lockheed(T_VD_MLI_R[n_node-1,layer], T_VD_MLI_R[n_node-1,layer+1], emissivity_MLI, layer_density_vector[layer], layers_vector[layer])[0]*prop.endcap_increment_area
                                dP_A_next_sub_L = Lockheed(T_VD_MLI_L[n_node-1,layer], T_VD_MLI_L[n_node-1,layer+1], emissivity_MLI, layer_density_vector[layer], layers_vector[layer])[0]*prop.endcap_increment_area
                                dP_from_previous_sub_R = dP_A_next_sub_R
                                dP_from_previous_sub_L = dP_A_next_sub_L

                            elif layer == VD_MLI_layers[-1]: #if last sector reached (end of VDMLI) go to i node
                                dP_i_prev_sub_R = - dP_from_previous_sub_R
                                dP_i_prev_sub_L = - dP_from_previous_sub_L
                            else:
                                dP_i_VDMLI_prev_sub_R = - dP_from_previous_sub_R 
                                dP_i_VDMLI_next_sub_R = Lockheed(T_VD_MLI_R[n_node-1,layer], T_VD_MLI_R[n_node-1,layer+1], emissivity_MLI, layer_density_vector[layer], layers_vector[layer])[0]*prop.endcap_increment_area

                                dP_i_VDMLI_prev_sub_L = - dP_from_previous_sub_L 
                                dP_i_VDMLI_next_sub_L = Lockheed(T_VD_MLI_L[n_node-1,layer], T_VD_MLI_L[n_node-1,layer+1], emissivity_MLI, layer_density_vector[layer], layers_vector[layer])[0]*prop.endcap_increment_area

                                dP_VD_MLI_R[n_node-1,layer] = dP_i_VDMLI_prev_sub_R + dP_i_VDMLI_next_sub_R
                                dP_VD_MLI_L[n_node-1,layer] = dP_i_VDMLI_prev_sub_L + dP_i_VDMLI_next_sub_L

                                dP_from_previous_sub_R = dP_i_VDMLI_next_sub_R
                                dP_from_previous_sub_L = dP_i_VDMLI_next_sub_L

                    else:
                        dP_A_next_sub_R = Lockheed(T_matrix_R[n_node-1,n_subnode], T_matrix_R[n_node-1,n_subnode+1], emissivity_MLI, layer_density, layers)[0]*prop.endcap_increment_area
                        dP_A_next_sub_L = Lockheed(T_matrix_L[n_node-1,n_subnode], T_matrix_L[n_node-1,n_subnode+1], emissivity_MLI, layer_density, layers)[0]*prop.endcap_increment_area
                        dP_i_prev_sub_R = - dP_A_next_sub_R
                        dP_i_prev_sub_L = - dP_A_next_sub_L 

                    # Opposite node
                    dP_A_opp_R = conduction(k_MLI, T_matrix_R[n_node-1,n_subnode]-T_matrix_L[n_node-1,n_subnode], prop.endcap_increment_length)*pi*(prop.radius+prop.thickness)*t_MLI
                    dP_A_opp_L = - dP_A_opp_R
 
                # NODE A (cylindrical part)
                else : 
                    # Previous SUBnode (space)
                    
                    # If absorbed fluxes from ESATAN are used
                    if radiative_case == 'Absorbed_Fluxes':
                        dP_sA_R = abs_flux_tot_R[step,n_node-1]*prop.increment_area
                        dP_sA_L = abs_flux_tot_L[step,n_node-1]*prop.increment_area

                    # If direct fluxes from ESATAN are used
                    elif radiative_case == 'Direct_Fluxes':
                        dP_sA_R = absorptivity*inc_flux_ISIA_R[step,n_node-1]*prop.increment_area \
                                + emissivity*inc_flux_IP_R[step,n_node-1]*prop.increment_area
                            
                        dP_sA_L = absorptivity*inc_flux_ISIA_L[step,n_node-1]*prop.increment_area \
                                + emissivity*inc_flux_IP_L[step,n_node-1]*prop.increment_area
                                
                    # heat emitted back by spacecraft
                    dP_As_R = radiation(emissivity,T_matrix_R[n_node-1,n_subnode],T_space)*prop.increment_area
                    dP_As_L = radiation(emissivity,T_matrix_L[n_node-1,n_subnode],T_space)*prop.increment_area

                    # Next SUBnode 
                    if VD_MLI == 'TRUE':
                        
                        for layer in VD_MLI_layers:

                            if layer == 0: 
                                dP_A_next_sub_R = Lockheed(T_VD_MLI_R[n_node-1,layer], T_VD_MLI_R[n_node-1,layer+1], emissivity_MLI, layer_density_vector[layer], layers_vector[layer])[0]*prop.increment_area
                                dP_A_next_sub_L = Lockheed(T_VD_MLI_L[n_node-1,layer], T_VD_MLI_L[n_node-1,layer+1], emissivity_MLI, layer_density_vector[layer], layers_vector[layer])[0]*prop.increment_area
                                dP_from_previous_sub_R = dP_A_next_sub_R
                                dP_from_previous_sub_L = dP_A_next_sub_L

                            elif layer == VD_MLI_layers[-1]: #if last sector reached (end of VDMLI) go to i node
                                dP_i_prev_sub_R = - dP_from_previous_sub_R
                                dP_i_prev_sub_L = - dP_from_previous_sub_L
                            else:
                                dP_i_VDMLI_prev_sub_R = - dP_from_previous_sub_R 
                                dP_i_VDMLI_next_sub_R = Lockheed(T_VD_MLI_R[n_node-1,layer], T_VD_MLI_R[n_node-1,layer+1], emissivity_MLI, layer_density_vector[layer], layers_vector[layer])[0]*prop.increment_area

                                dP_i_VDMLI_prev_sub_L = - dP_from_previous_sub_L 
                                dP_i_VDMLI_next_sub_L = Lockheed(T_VD_MLI_L[n_node-1,layer], T_VD_MLI_L[n_node-1,layer+1], emissivity_MLI, layer_density_vector[layer], layers_vector[layer])[0]*prop.increment_area

                                dP_VD_MLI_R[n_node-1,layer] = dP_i_VDMLI_prev_sub_R + dP_i_VDMLI_next_sub_R
                                dP_VD_MLI_L[n_node-1,layer] = dP_i_VDMLI_prev_sub_L + dP_i_VDMLI_next_sub_L

                                dP_from_previous_sub_R = dP_i_VDMLI_next_sub_R
                                dP_from_previous_sub_L = dP_i_VDMLI_next_sub_L

                    else:
                        dP_A_next_sub_R = Lockheed(T_matrix_R[n_node-1,n_subnode], T_matrix_R[n_node-1,n_subnode+1], emissivity_MLI, layer_density, layers)[0]*prop.increment_area
                        dP_A_next_sub_L = Lockheed(T_matrix_L[n_node-1,n_subnode], T_matrix_L[n_node-1,n_subnode+1], emissivity_MLI, layer_density, layers)[0]*prop.increment_area
                        dP_i_prev_sub_R = - dP_A_next_sub_R
                        dP_i_prev_sub_L = - dP_A_next_sub_L 

                    # Opposite node
                    dP_A_opp_R = conduction(k_MLI, T_matrix_R[n_node-1,n_subnode]-T_matrix_L[n_node-1,n_subnode], pi*(prop.radius+prop.thickness+t_MLI))*t_MLI*prop.increment_length
                    dP_A_opp_L = -dP_A_opp_R

                # Matrix for space fluxes
                dP_space [n_node-1,0] = dP_sA_R + dP_As_R
                dP_space [n_node-1,1] = dP_sA_L + dP_As_L
                dP_sA [n_node-1,0] = dP_sA_R
                dP_sA [n_node-1,1] = dP_sA_L
                dP_As [n_node-1,0] = dP_As_R
                dP_As [n_node-1,1] = dP_As_L
                
                if n_node != 1: # exclude case of top of the tank
                    # Previous node exists 
                    if n_node == 2:
                        length_increm = prop.endcap_increment_length
                    else:
                        length_increm = prop.increment_length

                    dP_A_prev_n_R = conduction(k_MLI, T_matrix_R[n_node-1,n_subnode]-T_matrix_R[n_node-2,n_subnode], length_increm)*pi*(prop.radius+prop.thickness)*t_MLI
                    dP_A_prev_n_L = conduction(k_MLI, T_matrix_L[n_node-1,n_subnode]-T_matrix_L[n_node-2,n_subnode], length_increm)*pi*(prop.radius+prop.thickness)*t_MLI

                if n_node != nr_nodes//2: # exclude case of bottom of the tank
                    # Next node exists
                    if n_node == nr_nodes//2-1:
                        length_increm = prop.endcap_increment_length
                    else:
                        length_increm = prop.increment_length

                    dP_A_next_n_R = conduction(k_MLI, T_matrix_R[n_node-1,n_subnode]-T_matrix_R[n_node,n_subnode], length_increm)*pi*(prop.radius+prop.thickness)*t_MLI
                    dP_A_next_n_L = conduction(k_MLI, T_matrix_L[n_node-1,n_subnode]-T_matrix_L[n_node,n_subnode], length_increm)*pi*(prop.radius+prop.thickness)*t_MLI

                # TOTAL 
                # Node A
                if n_node == 1:  # top of the tank, term of previous node missing
                    dP_net_matrix_R [n_node-1,n_subnode] = dP_sA_R + dP_As_R + dP_A_next_n_R + dP_A_next_sub_R + dP_A_opp_R
                    dP_net_matrix_L [n_node-1,n_subnode] = dP_sA_L + dP_As_L + dP_A_next_n_L + dP_A_next_sub_L  + dP_A_opp_L
                elif n_node == nr_nodes//2: # bottom of the tank, term of next node missing
                    dP_net_matrix_R [n_node-1,n_subnode] = dP_sA_R + dP_As_R + dP_A_prev_n_R + dP_A_next_sub_R + dP_A_opp_R
                    dP_net_matrix_L [n_node-1,n_subnode] = dP_sA_L + dP_As_L + dP_A_prev_n_L + dP_A_next_sub_L  + dP_A_opp_L
                else: 
                    dP_net_matrix_R [n_node-1,n_subnode] = dP_sA_R + dP_As_R + dP_A_next_n_R + dP_A_prev_n_R + dP_A_next_sub_R + dP_A_opp_R
                    dP_net_matrix_L [n_node-1,n_subnode] = dP_sA_L + dP_As_L + dP_A_next_n_L + dP_A_prev_n_L + dP_A_next_sub_L  + dP_A_opp_L
                
                if VD_MLI == 'TRUE': 
                    # first value of dP VDMLI is net power of node A
                    dP_VD_MLI_R[n_node-1,0] = dP_net_matrix_R [n_node-1,n_subnode]
                    dP_VD_MLI_L[n_node-1,0] = dP_net_matrix_L [n_node-1,n_subnode]


                n_subnode +=1 #update index for subnode (go to node i)
                subnode +=1 #update index for subnodes_vect
            ######################################################################################à

            if subnodes_vect [subnode] == 'WALL':

                if nr_nodes == 4: # for spherical tank

                    dP_Dl_R = h_liquid*-(T_matrix_R[n_node-1,n_subnode]-T_liquid)*prop.endcap_increment_area
                    dP_Dv_R = 0.

                    dP_Dl_L = h_liquid*-(T_matrix_L[n_node-1,n_subnode]-T_liquid)*prop.endcap_increment_area
                    dP_Dv_L = 0.

                else:

                    if n_node == 1 or n_node == nr_nodes//2: # this is the top of the tank, there is no previous node
                        
                        if tank_node == 'True':
                            # Next SUBnode (tank wall)
                            dP_C_next_sub_R = conduction(k_wall, T_matrix_R[n_node-1,n_subnode] - T_matrix_R[n_node-1,n_subnode+1], prop.thickness)*prop.endcap_increment_area
                            dP_C_next_sub_L = conduction(k_wall, T_matrix_L[n_node-1,n_subnode] - T_matrix_L[n_node-1,n_subnode+1], prop.thickness)*prop.endcap_increment_area
                            # Opposite node
                            dP_C_opp_R = conduction(k_wall, T_matrix_R[n_node-1,n_subnode] - T_matrix_L[n_node-1,n_subnode], prop.endcap_increment_length)*pi*(prop.radius)*prop.thickness
                            dP_C_opp_L = - dP_C_opp_R
                            
                            n_subnode +=1

                        # Next SUBnode (can be liquid or vapour based on fill level)
                        if n_node == 1:
                            if fill_level > (0.5*prop.endcap_internalvolume+prop.internalvolume)/V_total:
                                frac = max((fill_level-(0.5*prop.endcap_internalvolume+prop.internalvolume)/V_total)/(0.5*prop.endcap_internalvolume/V_total),1)
                                
                                dP_Dl_R = frac*(h_liquid*-(T_matrix_R[n_node-1,n_subnode]-T_liquid)*prop.endcap_increment_area)
                                dP_Dv_R = (1-frac)*(h_vapour*-(T_matrix_R[n_node-1,n_subnode]-T_vapour)*prop.endcap_increment_area)

                                dP_Dl_L = frac*(h_liquid*-(T_matrix_L[n_node-1,n_subnode]-T_liquid)*prop.endcap_increment_area)
                                dP_Dv_L = (1-frac)*(h_vapour*-(T_matrix_L[n_node-1,n_subnode]-T_vapour)*prop.endcap_increment_area)

                            else:
                                dP_Dv_R = h_vapour*-(T_matrix_R[n_node-1,n_subnode]-T_vapour)*prop.endcap_increment_area
                                dP_Dl_R = 0.

                                dP_Dv_L = h_vapour*-(T_matrix_L[n_node-1,n_subnode]-T_vapour)*prop.endcap_increment_area
                                dP_Dl_L = 0.

                        elif n_node == nr_nodes//2:
                            dP_Dl_R = h_liquid*-(T_matrix_R[n_node-1,n_subnode]-T_liquid)*prop.endcap_increment_area
                            dP_Dv_R = 0.

                            dP_Dl_L = h_liquid*-(T_matrix_L[n_node-1,n_subnode]-T_liquid)*prop.endcap_increment_area
                            dP_Dv_L = 0.

                    else:
                        if tank_node == 'True':                       
                            # Next SUBnode(tank wall)
                            dP_C_next_sub_R = conduction(k_wall, T_matrix_R[n_node-1,n_subnode] - T_matrix_R[n_node-1,n_subnode+1], prop.thickness)*prop.increment_area
                            dP_C_next_sub_L = conduction(k_wall, T_matrix_L[n_node-1,n_subnode] - T_matrix_L[n_node-1,n_subnode+1], prop.thickness)*prop.increment_area
                            # Opposite node
                            dP_C_opp_R = conduction(k_wall, T_matrix_R[n_node-1,n_subnode] - T_matrix_L[n_node-1,n_subnode], pi*(prop.radius))*prop.thickness*prop.increment_length
                            dP_C_opp_L = - dP_C_opp_R
    
                            n_subnode +=1

                        # Next SUBnode (can be liquid or vapour based on fill level)
                        if fill_level > (0.5*prop.endcap_internalvolume+cyl_sect*prop.internalvolume/(12/2-2))/V_total:
                            dP_Dl_R = h_liquid*-(T_matrix_R[n_node-1,n_subnode]-T_liquid)*prop.increment_area
                            dP_Dv_R = 0.

                            dP_Dl_L = h_liquid*-(T_matrix_L[n_node-1,n_subnode]-T_liquid)*prop.increment_area
                            dP_Dv_L = 0.

                        else:
                            dP_Dv_R = h_vapour*-(T_matrix_R[n_node-1,n_subnode]-T_vapour)*prop.increment_area
                            dP_Dl_R = 0.

                            dP_Dv_L = h_vapour*-(T_matrix_L[n_node-1,n_subnode]-T_vapour)*prop.increment_area
                            dP_Dl_L = 0.

                        cyl_sect -= 1; #prepare cylinder part for next node

                    if tank_node == 'True':

                        n_subnode -=1 #go back to node C (tank wall)
                        if n_node != 1: # previous node exists
                        # NODE C 
                        # Previous node 
                            if n_node == 2:
                                length_increm = prop.endcap_increment_length
                            else:
                                length_increm = prop.increment_length

                            dP_C_prev_n_R = conduction(k_wall, T_matrix_R[n_node-1,n_subnode] - T_matrix_R[n_node-2,n_subnode], length_increm)*pi*(prop.radius)*prop.thickness
                            dP_C_prev_n_L = conduction(k_wall, T_matrix_L[n_node-1,n_subnode] - T_matrix_L[n_node-2,n_subnode], length_increm)*pi*(prop.radius)*prop.thickness

                        if n_node != nr_nodes//2: # Next node exists
                            # NODE C
                            # Next node
                            if n_node == nr_nodes//2-1:
                                length_increm = prop.endcap_increment_length
                            else:
                                length_increm = prop.increment_length

                            dP_C_next_n_R = conduction(k_wall, T_matrix_R[n_node-1,n_subnode] - T_matrix_R[n_node,n_subnode], length_increm)*pi*(prop.radius)*prop.thickness
                            dP_C_next_n_L = conduction(k_wall, T_matrix_L[n_node-1,n_subnode] - T_matrix_L[n_node,n_subnode], length_increm)*pi*(prop.radius)*prop.thickness

                # Update dp matrix for node i temperature calculation
                if tank_node == 'True':

                    if n_node == 1:  # top of the tank, term of previous node missing
                        dP_net_matrix_R[n_node-1, n_subnode] = dP_i_prev_sub_R + dP_C_next_n_R + dP_C_next_sub_R + dP_C_opp_R 
                        dP_net_matrix_L[n_node-1, n_subnode] = dP_i_prev_sub_L + dP_C_next_n_L + dP_C_next_sub_L + dP_C_opp_L 

                    elif n_node == nr_nodes//2: # bottom of the tank, term of next node missing
                        dP_net_matrix_R[n_node-1, n_subnode] = dP_i_prev_sub_R + dP_C_prev_n_R + dP_C_next_sub_R + dP_C_opp_R
                        dP_net_matrix_L[n_node-1, n_subnode] = dP_i_prev_sub_L + dP_C_prev_n_L + dP_C_next_sub_L + dP_C_opp_L

                    else: 
                        dP_net_matrix_R[n_node-1, n_subnode] = dP_i_prev_sub_R + dP_C_prev_n_R + dP_C_next_n_R + dP_C_next_sub_R + dP_C_opp_R
                        dP_net_matrix_L[n_node-1, n_subnode] = dP_i_prev_sub_L + dP_C_prev_n_L + dP_C_next_n_L + dP_C_next_sub_L + dP_C_opp_L

                    
                    #balance for tank node
                    dP_net_matrix_R [n_node-1,n_subnode+1] = - dP_C_next_sub_R + dP_Dl_R + dP_Dv_R 
                    dP_net_matrix_L [n_node-1,n_subnode+1] = - dP_C_next_sub_L + dP_Dl_L + dP_Dv_L 

                    if VD_MLI == 'TRUE':
                        # last value of dP VDMLI is net power of node i 
                        dP_VD_MLI_R[n_node-1,-1] = dP_net_matrix_R[n_node-1, n_subnode]
                        dP_VD_MLI_L[n_node-1,-1] = dP_net_matrix_L[n_node-1, n_subnode]


                else:  # no tank node is present  
                    dP_net_matrix_R [n_node-1,n_subnode] = dP_i_prev_sub_R + dP_Dl_R + dP_Dv_R 
                    dP_net_matrix_L [n_node-1,n_subnode] = dP_i_prev_sub_L + dP_Dl_L + dP_Dv_L 

                    if VD_MLI == 'TRUE':
                        # last value of dP VDMLI is net power of node i 
                        dP_VD_MLI_R[n_node-1,-1] = dP_i_prev_sub_R + dP_Dl_R + dP_Dv_R 
                        dP_VD_MLI_L[n_node-1,-1] = dP_i_prev_sub_L + dP_Dl_L + dP_Dv_L 

                n_subnode +=1 #update index for subnode (actually not needed but added for completeness) 
                subnode +=1 #update index for subnodes_vect (actually not needed but added for completeness)

        ############################  END of calculation FOR SUBNODES #########################################################

            # Net values per node (sum of each row of the matrices), to calculate once each subnode has been analysed
            #dP_node_net [n_node-1,0] = np.sum(dP_net_matrix_R[n_node-1,:])  # nodes on the right side of tank (R)
            #dP_node_net [n_node-1,1] = np.sum(dP_net_matrix_L[n_node-1,:]) # nodes on the left side of tank (L)
            
            # Once insulation sequence is complete, calculate heat fluxes to LIQUID and VAPOUR nodes 
            # (these need to be treated as two nodes, they go in the for loop for nodes)
                            
            ######## LIQUID NODE ########
            # NODE l(iquid)
            dP_liquid_matrix[n_node-1,0] = - dP_Dl_R 
            dP_liquid_matrix[n_node-1,1] = - dP_Dl_L    
            ########################

            ######## VAPOUR NODE ########
            # NODE v(apour)
            dP_vapour_matrix[n_node-1,0] = - dP_Dv_R 
            dP_vapour_matrix[n_node-1,1] = - dP_Dv_L
            ########################
            
            n_node +=1 
        ############################  END of FOR LOOP FOR NODES #########################################################

        # Other (fixed) power inputs for liquid and vapour nodes
        # Liquid-vapour interface
        dP_lv = h_liquid* -(T_liquid - T_vapour) * pi * prop.radius ** 2.
        dP_vl = -dP_lv
        # Cooler power
        dP_cooler = -Q_cooler
        # Mixer power (not considered here)
        dP_mixer = Q_mixer
        # Penetration leaks (not considered here)
        dP_penetration = 0.0025*-(T_liquid-T_space)*(V_total)**(0.5)
        
        dP_liquid_net = np.sum(dP_liquid_matrix) + dP_lv + dP_cooler #+ dP_penetration #+dP_mixer
        dP_vapour_net = np.sum(dP_vapour_matrix) + dP_vl

        P_rest = 0. #initialised for later

        # # ALL CONDUCTION IN THE SYSTEM PLUS RADIATION IN AND OUT
        #dP_system = np.sum(dP_node_net) + dP_liquid_net+dP_vapour_net
        # ALL RADIATION FROM SPACE TO THE TANK DUE TO THE SUN
        #dP_in = np.sum(dP_sA)
        # ALL RADIATION FROM THE TANK BACK TO SPACE DUE TO IR
        dP_out = np.sum(dP_As)

        #dP_tot = dP_system-dP_in-dP_out

        ######### Temperature increase ###########################
        if nr_nodes > 4:
            dT_matrix_R = np.divide(dP_net_matrix_R, np.multiply(Cp_vect,Mass_sect_cyl))
            dT_matrix_L = np.divide(dP_net_matrix_L, np.multiply(Cp_vect,Mass_sect_cyl))

        # Redo and overwrite for endcaps
        #Top of tank
        dT_matrix_R[0,:] = np.divide(dP_net_matrix_R[0,:], np.multiply(Cp_vect,Mass_sect_sph))
        dT_matrix_L[0,:] = np.divide(dP_net_matrix_L[0,:], np.multiply(Cp_vect,Mass_sect_sph))

        #Bottom of tank
        dT_matrix_R[-1,:] = np.divide(dP_net_matrix_R[-1,:], np.multiply(Cp_vect,Mass_sect_sph))
        dT_matrix_L[-1,:] = np.divide(dP_net_matrix_L[-1,:], np.multiply(Cp_vect,Mass_sect_sph))

        # calculate new node and subnodes temperatures
        T_matrix_R = T_matrix_R + dt*dT_matrix_R
        T_matrix_L = T_matrix_L + dt*dT_matrix_L

        if VD_MLI == 'TRUE':
            dT_matrix_VDMLI_R = np.divide(dP_VD_MLI_R, np.multiply(Cp_MLI,Mass_sect_cyl_VDMLI)) 
            dT_matrix_VDMLI_L = np.divide(dP_VD_MLI_L, np.multiply(Cp_MLI,Mass_sect_cyl_VDMLI))

            # Redo and overwrite for endcaps
            #Top of tank
            dT_matrix_VDMLI_R[0,:] = np.divide(dP_VD_MLI_R[0,:], np.multiply(Cp_MLI,Mass_sect_sph_VDMLI))
            dT_matrix_VDMLI_L[0,:] = np.divide(dP_VD_MLI_L[0,:], np.multiply(Cp_MLI,Mass_sect_sph_VDMLI))

            #Bottom of tank
            dT_matrix_VDMLI_R[5,:] = np.divide(dP_VD_MLI_R[5,:], np.multiply(Cp_MLI,Mass_sect_sph_VDMLI))
            dT_matrix_VDMLI_L[5,:] = np.divide(dP_VD_MLI_L[5,:], np.multiply(Cp_MLI,Mass_sect_sph_VDMLI))

            T_VD_MLI_R = T_VD_MLI_R + dt*dT_matrix_VDMLI_R 
            T_VD_MLI_L = T_VD_MLI_L + dt*dT_matrix_VDMLI_L 

            # in case of VDMLI, node i Temperature is the last column of T VD MLI and node A temperature is the first coloumn of T VD MLI
            T_matrix_R[:,1] = T_VD_MLI_R[:,-1]
            T_matrix_L[:,1] = T_VD_MLI_L[:,-1]

            T_matrix_R[:,0] = T_VD_MLI_R[:,0]
            T_matrix_L[:,0] = T_VD_MLI_L[:,0]  
        ###################################################################  

        ########## Pressure control and venting calculations ################
        if mode == 'HeliumPress':
            T_sat = CP.PropsSI('T','P',max(p_vapour_initial, min(p_max, p_vapour)),'Q',0,fluid)
            Cp_prop_vapour = CP.PropsSI('CP0MASS','P',p_vapour,'T',T_vapour,fluid)

            Cp_mix = (m_boil_off*Cp_prop_vapour+m_vapour_initial*Cp_vapour)/m_vapour
            m_mix = (m_boil_off*m_molecular_propell+m_vapour_initial*m_molecular_vapour)/m_vapour*1e-3 # m molecular is in g/mol

        elif mode == 'Autogenous':
            Cp_mix = Cp_vapour
            m_mix = m_molecular_vapour*1e-3 #kg/mol

        if T_liquid < T_sat:

            dT_liquid = dP_liquid_net*dt/(Cp_liquid*m_liquid)
            T_liquid = T_liquid+dT_liquid
            dm_boil_off = 0

            if T_liquid > T_sat:
                if t_boil == 0.:
                    t_boil = t/60./60./24.

                dT_rest = T_liquid-T_sat
                P_rest = dT_rest*Cp_liquid*m_liquid
                dm_boil_off = P_rest/L_liquid/1000.
                m_boil_off = m_boil_off+dm_boil_off
                m_liquid = m_liquid-dm_boil_off
                T_liquid = T_sat

            V_liquid = m_liquid/density_liquid
            fill_level = V_liquid/V_total

            m_vapour = m_vapour + dm_boil_off
            V_vapour = V_total-V_liquid
            dT_vapour = dP_vapour_net*dt/(Cp_mix*m_vapour)
            T_vapour = T_vapour+dT_vapour
            dp_vapour = pressure_increase(T_sat,L_liquid,m_molecular_propell,V_vapour,P_rest) + (R_G*dP_vapour_net*dt/(m_mix*V_vapour*Cp_mix))
            p_vapour = p_vapour + dp_vapour

        elif T_liquid >= T_sat:
            if t_boil == 0.:
                t_boil = t/60./60./24.

            dm_boil_off = dP_liquid_net*dt/L_liquid/1000.
            m_boil_off = m_boil_off+dm_boil_off
            m_liquid = m_liquid - dm_boil_off

            V_liquid = m_liquid/density_liquid
            fill_level = V_liquid/V_total
            
            m_vapour = m_vapour + dm_boil_off

            dT_vapour = dP_vapour_net*dt/(Cp_mix*m_vapour)
            T_vapour = T_vapour+dT_vapour
            V_vapour = V_total-V_liquid
            
            #Calculate pressure increase in ullage due to boil-off vapours and thermal expansion of existing gas
            dp_vapour = pressure_increase(T_sat,L_liquid,m_molecular_propell,V_vapour,dP_liquid_net*dt) + (R_G*dP_vapour_net*dt/(m_mix*V_vapour*Cp_mix))
            p_vapour = p_vapour + dp_vapour
            
        ##### Now check if overpressure is reached and if yes, start venting #################
        if p_vapour > p_max and venting_started == 'FALSE':
            vent_number +=1
            if t_vent == 0:
                t_vent = t # value is returned in seconds

            venting_started = 'TRUE'
            gamma = CP.PropsSI('ISENTROPIC_EXPANSION_COEFFICIENT','P',p_vapour,'Q',1,fluid)

            T_vapour_after_vent = T_vapour*(p_vapour/p_vent_target)**((1.-gamma)/gamma) # New vapour temperature after expansion in ullage due to venting
            m_vapour_new = p_vent_target*V_vapour*m_mix/(R_G*T_vapour_after_vent) # calculate new vapour mass in ullage after venting (already includes the effect of evaporative cooling liquid)
            dm_vent = m_vapour - m_vapour_new #and vented mass

            # evaporation cooling of liquid
            if mode == 'Autogenous':
                T_liquid_new = CP.PropsSI('T','P',p_vent_target,'Q',0,fluid) #new liquid temperature will be saturation T at new pressure p_vent
             
            venting_rate = dm_vent/dt
 
            if single_venting == 'TRUE':
                T_vapour = T_vapour_after_vent
                m_vent = m_vent + dm_vent #update vent mass
                m_vapour = m_vapour_new #update vapour mass (equals to m vapour new)  
                p_vapour = p_vent_target # set by user
                if mode == 'Autogenous':
                    T_liquid = T_liquid_new
                    T_sat = CP.PropsSI('T','P',p_vapour,'Q',0,fluid)
                venting_started = 'FALSE' #venting process is over

        elif p_vapour <= p_max and venting_started == 'FALSE':
                dm_vent = 0 # if maximum pressure is not reached, set dm_vent = 0
            
        if m_liquid <= 0. or V_liquid <= 0.:
            print('V_liquid = ', V_liquid, 'm_liquid =', m_liquid)
            break
        elif m_vapour < 0. or V_vapour < 0.:
            print('V_vapour =', V_vapour, 'm_vapour =', m_vapour)
            break
        ##################################################################################

        ######## VCS COMPUTATIONS
        # VCS has uniform temperature, this meaning that conductive heat transfer to surrounding nodes does not happen!
        # the only heat transfer to be calculated is the one between fluid and shield
        if dm_vent > 0 and VCS_on == 'TRUE':

            if VD_MLI == 'TRUE':
                if VCS_position == 'middle':
                    dP_in_VCS = np.sum(dP_VD_MLI_R[:,1] + dP_VD_MLI_L[:,1]) #VCS is placed where middle MLI is
                elif VCS_position == 'inner':
                    dP_in_VCS = np.sum(dP_VD_MLI_R[:,2] + dP_VD_MLI_L[:,2]) #VCS is placed where inner MLI is

            else: 
                # TOTAL Power input from last MLI layer to VCS shield for all nodes 
                dP_in_VCS = np.sum(dP_net_matrix_R[:,1] + dP_net_matrix_L[:,1])

            m_rate = venting_rate

            # Values initialization
            mu_VCS = CP.PropsSI('VISCOSITY','P',p_vapour,'Q',1,fluid) #Pa*s
            Cp_VCS = CP.PropsSI('CP0MASS','P',p_vapour,'Q',1,fluid) 
            k_VCS = CP.PropsSI('CONDUCTIVITY','P',p_vapour,'Q',1,fluid) 

            if VD_MLI == 'TRUE':
                if VCS_position == 'middle':
                    T_VCS_old = np.average(np.maximum(T_VD_MLI_R[:,1],T_VD_MLI_L[:,1])) # K average of temperatures of maximum between R and L MLI middle nodes
                elif VCS_position == 'inner': 
                    T_VCS_old = np.average(np.maximum(T_VD_MLI_R[:,2],T_VD_MLI_L[:,2])) # K average of temperatures of maximum between R and L MLI inner nodes

            else: 
                T_VCS_old = np.average(np.append(T_matrix_R[:,1],T_matrix_L[:,1])) # K average of temperatures of all i nodes

            print('T VCS old is', T_VCS_old)
            T_fluid_out = T_sat

            # VCS geometry and fluid properties
            A_cross = pi*d_VCS**2./4.
            A_VCS = pi*d_VCS*length_VCS

            Re_VCS = 4*m_rate/(mu_VCS*pi*d_VCS)
            Pr_VCS = Cp_VCS*mu_VCS/k_VCS
            Nu_VCS = NusseltfromRe(Re_VCS,Pr_VCS)
            h_VCS = Nu_VCS*k_VCS/d_VCS

            # Initialization
            tt = 0
            T_fluid_in = T_sat

            while tt < dt:
                tt = tt + dt_VCS

                # forward Euler
                dP_VCS = h_VCS*A_VCS*(T_VCS_old - (T_fluid_in + T_fluid_out)/2.)
                T_VCS_new = T_VCS_old + dt_VCS*(dP_in_VCS - dP_VCS)/(M_shield*Cp_shield)
                T_fluid_out = T_fluid_in + dt_VCS*dP_VCS/(m_rate*Cp_VCS) #

                T_VCS_old = T_VCS_new

            T_VCS = T_VCS_new # while loop returns "stabilized" VCS temperature

            # print('VCS temperature is', T_VCS,'dP VCS is', dP_VCS, 'dP in VCS is', dP_in_VCS,'T in fluid', T_fluid_in, 'T_out_fluid', T_fluid_out,'venting rate is',m_rate)
            
            if T_VCS < 0:
                print('T VCS', T_VCS)
                break

            # now update temperatures of insulation layers to be used in next iteration 
            if VD_MLI == 'TRUE': # middle MLI layer is VCS
                if VCS_position == 'middle':
                    T_VD_MLI_R [:,1] = T_VCS
                    T_VD_MLI_L [:,1] = T_VCS
                elif VCS_position == 'inner':
                    T_VD_MLI_R [:,2] = T_VCS
                    T_VD_MLI_L [:,2] = T_VCS
                    
            elif VD_MLI == 'FALSE':#node i temperature is T VCS
                T_matrix_R[:,1] = T_VCS
                T_matrix_L[:,1] = T_VCS

        # save values for plots
        if save_vectors == 'TRUE':
            time_vect[nn] = t
            dm_vent_vector [nn] = dm_vent 
            m_boiloff_vector [nn] = dm_boil_off

            p_vapour_vector [nn] = p_vapour
            m_vapour_vector [nn] = m_vapour

            T_vapour_vect [nn] = round(T_vapour,1)
            T_liquid_vect [nn] = round(T_liquid,1)
            T_VCS_vect [nn] = T_VD_MLI_L [2,1]

            dP_vapour_vector [nn] = dP_vapour_net
            dP_liquid_vector [nn] = dP_liquid_net

            dP_radiative_net [nn] = np.sum(dP_space) #net heat incoming from space (absorbed + emitted)
            dP_emitted_back [nn] = dP_out

        nn +=1

        ###################### Write on Excel if 1,2,3,4,5,6 months are reached #####################################################
   
        if t >= (60.*60.*24.*30.*1.) and one_month == False:
            # initialise values
            m_boil_off_previous_month = 0.
            m_liquid_initial_previous_month = m_liquid_initial  

            Boil_off_monthly = (m_boil_off/m_liquid_initial/t*60.*60.*24.*30.)*100.
            Boiloff_from_start_of_month = ((m_boil_off - m_boil_off_previous_month)/m_liquid_initial_previous_month)*100.
            with open('mass_and_boil-off_v15-1m.csv', 'a') as filewriter:
                writer = csv.writer(filewriter)
                writer.writerow([layer_density, layers, m_MLI, Q_cooler, m_cooler, prop.thickness, m_shell, m_boil_off, m_vent, m_boil_off+m_MLI+m_shell+m_cooler, Boil_off_monthly, t_boil, (t/60./60./24.), T_liquid, T_vapour,  Boiloff_from_start_of_month])
            filewriter.close()
            #now save final temperatures and heat leaks
            filename = 'T_and_dP_Data-{}m-option{}'.format(months,option_nr) 
            np.savez(os.path.join(BASE_PATH,filename), T_matrix_R = T_matrix_R, T_matrix_L = T_matrix_L,  dP_net_matrix_R = dP_net_matrix_R, dP_net_matrix_L = dP_net_matrix_L)
            #save data for next month boiloff calculation
            m_boil_off_previous_month = m_boil_off
            m_liquid_initial_previous_month = m_liquid

            one_month = True

        elif t >= (60.*60.*24.*30.*2.) and two_months == False:
            Boil_off_monthly = (m_boil_off/m_liquid_initial/t*60.*60.*24.*30.)*100.
            Boiloff_from_start_of_month = ((m_boil_off - m_boil_off_previous_month)/m_liquid_initial_previous_month)*100.
            with open('mass_and_boil-off_v15-2m.csv', 'a') as filewriter:
                writer = csv.writer(filewriter)
                writer.writerow([layer_density, layers, m_MLI, Q_cooler, m_cooler, prop.thickness, m_shell, m_boil_off, m_vent, m_boil_off+m_MLI+m_shell+m_cooler, Boil_off_monthly, t_boil, (t/60./60./24.), T_liquid, T_vapour, Boiloff_from_start_of_month])
            filewriter.close()
            #now save final temperatures and heat leaks
            filename = 'T_and_dP_Data-{}m-option{}'.format(months,option_nr) 
            np.savez(os.path.join(BASE_PATH,filename), T_matrix_R = T_matrix_R, T_matrix_L = T_matrix_L,  dP_net_matrix_R = dP_net_matrix_R, dP_net_matrix_L = dP_net_matrix_L)
            #save data for next month boiloff calculation
            m_boil_off_previous_month = m_boil_off
            m_liquid_initial_previous_month = m_liquid

            two_months = True

        elif t >= (60.*60.*24.*30.*3.) and three_months == False:
            Boil_off_monthly = (m_boil_off/m_liquid_initial/t*60.*60.*24.*30.)*100.
            Boiloff_from_start_of_month = ((m_boil_off - m_boil_off_previous_month)/m_liquid_initial_previous_month)*100.
            with open('mass_and_boil-off_v15-3m.csv', 'a') as filewriter:
                writer = csv.writer(filewriter)
                writer.writerow([layer_density, layers, m_MLI, Q_cooler, m_cooler, prop.thickness, m_shell, m_boil_off, m_vent, m_boil_off+m_MLI+m_shell+m_cooler, Boil_off_monthly, t_boil, (t/60./60./24.), T_liquid, T_vapour, Boiloff_from_start_of_month])            
            filewriter.close()
            #now save final temperatures and heat leaks
            filename = 'T_and_dP_Data-{}m-option{}'.format(months,option_nr) 
            np.savez(os.path.join(BASE_PATH,filename), T_matrix_R = T_matrix_R, T_matrix_L = T_matrix_L,  dP_net_matrix_R = dP_net_matrix_R, dP_net_matrix_L = dP_net_matrix_L)
            #save data for next month boiloff calculation
            m_boil_off_previous_month = m_boil_off
            m_liquid_initial_previous_month = m_liquid

            three_months = True

        elif t >= (60.*60.*24.*30.*4.) and four_months == False:
            Boil_off_monthly = (m_boil_off/m_liquid_initial/t*60.*60.*24.*30.)*100.
            Boiloff_from_start_of_month = ((m_boil_off - m_boil_off_previous_month)/m_liquid_initial_previous_month)*100.
            with open('mass_and_boil-off_v15-4m.csv', 'a') as filewriter:
                writer = csv.writer(filewriter)
                writer.writerow([layer_density, layers, m_MLI, Q_cooler, m_cooler, prop.thickness, m_shell, m_boil_off, m_vent, m_boil_off+m_MLI+m_shell+m_cooler, Boil_off_monthly, t_boil, (t/60./60./24.), T_liquid, T_vapour, Boiloff_from_start_of_month])            
            filewriter.close()
            #now save final temperatures and heat leaks
            filename = 'T_and_dP_Data-{}m-option{}'.format(months,option_nr) 
            np.savez(os.path.join(BASE_PATH,filename), T_matrix_R = T_matrix_R, T_matrix_L = T_matrix_L,  dP_net_matrix_R = dP_net_matrix_R, dP_net_matrix_L = dP_net_matrix_L)
            #save data for next month boiloff calculation
            m_boil_off_previous_month = m_boil_off
            m_liquid_initial_previous_month = m_liquid
            four_months = True

        elif t >= (60.*60.*24.*30.*5.) and five_months == False:
            Boil_off_monthly = (m_boil_off/m_liquid_initial/t*60.*60.*24.*30.)*100.
            Boiloff_from_start_of_month = ((m_boil_off - m_boil_off_previous_month)/m_liquid_initial_previous_month)*100.
            with open('mass_and_boil-off_v15-5m.csv', 'a') as filewriter:
                writer = csv.writer(filewriter)
                writer.writerow([layer_density, layers, m_MLI, Q_cooler, m_cooler, prop.thickness, m_shell, m_boil_off, m_vent, m_boil_off+m_MLI+m_shell+m_cooler, Boil_off_monthly, t_boil, (t/60./60./24.), T_liquid, T_vapour, Boiloff_from_start_of_month])            
            filewriter.close()
            #now save final temperatures and heat leaks
            filename = 'T_and_dP_Data-{}m-option{}'.format(months,option_nr) 
            np.savez(os.path.join(BASE_PATH,filename), T_matrix_R = T_matrix_R, T_matrix_L = T_matrix_L,  dP_net_matrix_R = dP_net_matrix_R, dP_net_matrix_L = dP_net_matrix_L)
            #save data for next month boiloff calculation
            m_boil_off_previous_month = m_boil_off
            m_liquid_initial_previous_month = m_liquid
            five_months = True

        elif t >= (60.*60.*24.*30.*6.) and six_months == False:
            Boil_off_monthly = (m_boil_off/m_liquid_initial/t*60.*60.*24.*30.)*100.
            Boiloff_from_start_of_month = ((m_boil_off - m_boil_off_previous_month)/m_liquid_initial_previous_month)*100.
            with open('mass_and_boil-off_v15-6m.csv', 'a') as filewriter:
                writer = csv.writer(filewriter)
                writer.writerow([layer_density, layers, m_MLI, Q_cooler, m_cooler, prop.thickness, m_shell, m_boil_off, m_vent, m_boil_off+m_MLI+m_shell+m_cooler, Boil_off_monthly, t_boil, (t/60./60./24.), T_liquid, T_vapour, Boiloff_from_start_of_month])            
            filewriter.close()
            #now save final temperatures and heat leaks
            filename = 'T_and_dP_Data-{}m-option{}'.format(months,option_nr) 
            np.savez(os.path.join(BASE_PATH,filename), T_matrix_R = T_matrix_R, T_matrix_L = T_matrix_L,  dP_net_matrix_R = dP_net_matrix_R, dP_net_matrix_L = dP_net_matrix_L)
            #save data for next month boiloff calculation
            m_boil_off_previous_month = m_boil_off
            m_liquid_initial_previous_month = m_liquid
            six_months = True

    ############## Write final results on Excel ##################
    
    Boil_off_monthly = (m_boil_off/m_liquid_initial/t*60.*60.*24.*30.)*100.
    Boiloff_from_start_of_month = 0 #this parameter was only saved for consecutive months of storage

    if months != 1 and months < 2:
        with open('mass_and_boil-off_v15-%sm.csv' %months, 'a') as filewriter:
            writer = csv.writer(filewriter, delimiter=',')
            print('you entered the right line')
            writer.writerow([layer_density, layers, m_MLI, Q_cooler, m_cooler, prop.thickness, m_shell, m_boil_off, m_vent, m_boil_off+m_MLI+m_shell+m_cooler, Boil_off_monthly, t_boil, (t/60./60./24.), T_liquid, T_vapour, Boiloff_from_start_of_month])
        filewriter.close()
        #now save final temperatures and heat leaks
        filename = 'T_and_dP_Data-{}m-option{}'.format(months,option_nr) 
        np.savez(os.path.join(BASE_PATH,filename), T_matrix_R = T_matrix_R, T_matrix_L = T_matrix_L,  dP_net_matrix_R = dP_net_matrix_R, dP_net_matrix_L = dP_net_matrix_L)

    elif months > 2 and months < 3:
        with open('mass_and_boil-off_v15-%sm.csv' %months, 'a') as filewriter:
            writer = csv.writer(filewriter)
            writer.writerow([layer_density, layers, m_MLI, Q_cooler, m_cooler, prop.thickness, m_shell, m_boil_off, m_vent, m_boil_off+m_MLI+m_shell+m_cooler, Boil_off_monthly, t_boil, (t/60./60./24.), T_liquid, T_vapour, Boiloff_from_start_of_month])        
        filewriter.close()
        #now save final temperatures and heat leaks
        filename = 'T_and_dP_Data-{}m-option{}'.format(months,option_nr) 
        np.savez(os.path.join(BASE_PATH,filename), T_matrix_R = T_matrix_R, T_matrix_L = T_matrix_L,  dP_net_matrix_R = dP_net_matrix_R, dP_net_matrix_L = dP_net_matrix_L)

    elif months > 3 and months < 4:
        with open('mass_and_boil-off_v15-%sm.csv' %months, 'a') as filewriter:
            writer = csv.writer(filewriter)
            writer.writerow([layer_density, layers, m_MLI, Q_cooler, m_cooler, prop.thickness, m_shell, m_boil_off, m_vent, m_boil_off+m_MLI+m_shell+m_cooler, Boil_off_monthly, t_boil, (t/60./60./24.), T_liquid, T_vapour, Boiloff_from_start_of_month])        
        filewriter.close()
        #now save final temperatures and heat leaks
        filename = 'T_and_dP_Data-{}m-option{}'.format(months,option_nr) 
        np.savez(os.path.join(BASE_PATH,filename), T_matrix_R = T_matrix_R, T_matrix_L = T_matrix_L,  dP_net_matrix_R = dP_net_matrix_R, dP_net_matrix_L = dP_net_matrix_L)

    elif months > 4 and months < 5:
        with open('mass_and_boil-off_v15-%sm.csv' %months, 'a') as filewriter:
            writer = csv.writer(filewriter)
            writer.writerow([layer_density, layers, m_MLI, Q_cooler, m_cooler, prop.thickness, m_shell, m_boil_off, m_vent, m_boil_off+m_MLI+m_shell+m_cooler, Boil_off_monthly, t_boil, (t/60./60./24.), T_liquid, T_vapour, Boiloff_from_start_of_month])        
        filewriter.close()
        #now save final temperatures and heat leaks
        filename = 'T_and_dP_Data-{}m-option{}'.format(months,option_nr) 
        np.savez(os.path.join(BASE_PATH,filename), T_matrix_R = T_matrix_R, T_matrix_L = T_matrix_L,  dP_net_matrix_R = dP_net_matrix_R, dP_net_matrix_L = dP_net_matrix_L)

    elif months > 5:
        with open('mass_and_boil-off_v15-%sm.csv' %months, 'a') as filewriter:
            writer = csv.writer(filewriter)
            writer.writerow([layer_density, layers, m_MLI, Q_cooler, m_cooler, prop.thickness, m_shell, m_boil_off, m_vent, m_boil_off+m_MLI+m_shell+m_cooler, Boil_off_monthly, t_boil, (t/60./60./24.), T_liquid, T_vapour, Boiloff_from_start_of_month])        
        filewriter.close()
        #now save final temperatures and heat leaks
        filename = 'T_and_dP_Data-{}m-option{}'.format(months,option_nr) 
        np.savez(os.path.join(BASE_PATH,filename), T_matrix_R = T_matrix_R, T_matrix_L = T_matrix_L,  dP_net_matrix_R = dP_net_matrix_R, dP_net_matrix_L = dP_net_matrix_L)


    ####### Reduce plot vectors to actual size
    if save_vectors == 'TRUE':
        time_vect = time_vect [0:nn]

        dm_vent_vector = dm_vent_vector [0:nn]
        m_boiloff_vector = m_boiloff_vector [0:nn]

        p_vapour_vector = p_vapour_vector [0:nn]
        m_vapour_vector = m_vapour_vector [0:nn]

        T_vapour_vect = T_vapour_vect [0:nn]
        T_liquid_vect = T_liquid_vect [0:nn]
        T_VCS_vect = T_VCS_vect[0:nn]
        
        dP_vapour_vector = dP_vapour_vector [0:nn]
        dP_liquid_vector = dP_liquid_vector [0:nn]

        dP_radiative_net = dP_radiative_net[0:nn]
        dP_emitted_back = dP_emitted_back[0:nn]

   ##############################################

    print('Venting starts at day : ', t_vent/(60*60*24))
    print('Boiloff starts at day: ', t_boil/(60*60*24))
    print('total heat leak into liquid and vapour is:', dP_liquid_net, dP_vapour_net,'corresponding to:', ((dP_liquid_net+dP_vapour_net)*(60.*60.*24.*30.)/((L_liquid*1e3)*m_liquid_initial))*100.,'and L_liquid and m initial',L_liquid, m_liquid_initial)
    
################# Save data in .npz files for plots ##################################

if VD_MLI == 'TRUE':
    thickness_vector=np.append([t_wall],t_MLI_vector)
else:
    thickness_vector = np.array([t_wall,t_MLI])

##################################

np.savez(os.path.join(BASE_PATH,'Tank_geometry'), tank_endcap = tank_endcap, sphere_parts = prop.sphere_parts, cylinder_parts = prop.cylinder_parts, radius = prop.radius, height = prop.height, surf_area_cap_node = prop.endcap_increment_area, surf_area_cyl_node = prop.increment_area)
if save_vectors == 'TRUE':
    np.savez(os.path.join(BASE_PATH,'PlotDataR'), T_VD_MLI_R = T_VD_MLI_R, T_matrix_R = T_matrix_R, time_vect = time_vect, thickness_vector = thickness_vector)
    np.savez(os.path.join(BASE_PATH,'PlotDataL'), T_VD_MLI_L = T_VD_MLI_L, T_matrix_L = T_matrix_L)
    np.savez(os.path.join(BASE_PATH,'dPData'), dP_vapour_vector = dP_vapour_vector, dP_liquid_vector = dP_liquid_vector, dP_radiative_net = dP_radiative_net, dP_VD_MLI_R = dP_VD_MLI_R, dP_VD_MLI_L = dP_VD_MLI_L, dP_emitted_back = dP_emitted_back) 
    np.savez(os.path.join(BASE_PATH,'VentData'), p_vapour_vector = p_vapour_vector, m_vapour_vector = m_vapour_vector, T_vapour_vect = T_vapour_vect, T_liquid_vect = T_liquid_vect, dm_vent_vector = dm_vent_vector, m_boiloff_vector = m_boiloff_vector, T_VCS_vect = T_VCS_vect)
    np.savez(os.path.join(BASE_PATH,'MLIcheck'), T_matrix_L = T_matrix_L, T_liquid_vect = T_liquid_vect, T_VCS_vect = T_VCS_vect)
