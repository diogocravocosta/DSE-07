
import pandas as pd
import csv
import itertools as it
from math import *
import pathlib
import numpy as np
from tqdm.auto import tqdm
import BO_tool_main_body

################## MAIN DESIGN PARAMETERS ##############################
# parameters that can be vectors are indicated with []
cooler_power = [0.] 

dt = 4.
months = 0.25

MLI_layers = [ 45.] #[ 10., 20., 25., 35., 45., 65.]
layer_densities = [ 16.] #[16., 40.]

heat_fluxes_ESATAN = 'Absorbed_Fluxes'# 'Absorbed_Fluxes' or 'Direct_Fluxes'
Radiative_model_import = 'Radiative_model_GEO_single.npz' #put name of Radiative model from ESATAN to be used

press_mode = 'Autogenous' # 'HeliumPress' or 'Autogenous'

VD_MLI_activity = 'TRUE' # 'TRUE' or 'FALSE'. If TRUE, set VDMLI data below

VCS_activity = 'TRUE' # 'TRUE' or 'FALSE'
VCS_location = 'middle' # 'middle' or 'inner'

include_tank_node = 'False' # 'TRUE' or 'FALSE'. If false, remember to also remove "WALL" from subnodes_structure
subnodes_structure = ["MLI","MLI"] # ["MLI","MLI"] or ["MLI","MLI","WALL","WALL"]

save_vectors_for_plots = 'TRUE' #set as false if t run is more than 3 months or if design options are multiple

############################################################################################

##### Other design parameters #######################

# convert storage duration in seconds
t_run = months*2592000.

# Propellant properties
propellant_mixture = 'LH2' # Boil-off tool has been written for'LH2', however in theory it can work with any fluid in the CoolProp library
propellant_mass = 34600  # kg
fill_level_initial = 0.9
initial_pressure = 1.3*1e5 # Pa 
p_max = 3*1e5 # Pa
maximum_tank_radius = 2.7 # m

initial_temperature = 20. #K, not relevant for autogenous pressurization case but set it for initialization of values

# Tank geometry and nodal structure
tank_heads = 'sphere' #alternative is 'ellipse'
cap_height = 0.5*maximum_tank_radius # only applicable for 'ellipse' tank head, cap_height should be always less than tank radius!!! (oblate spheroid)
sphere_parts_tank = 4.
cylinder_parts_tank = 8.
total_number_nodes = 12 #sphere_parts_tank+cylinder_parts_tank

VD_MLI_number_of_layers = [0,1,2] # outer (A), middle and inner(i)
insulation_structure = ["MLI", "WALL"]

# MLI properties
layer_specific_weight = 0.00881 #kg/m2 Double aluminized Mylar
spacer_specific_weight = 0.00635 #kg/m2 Dacron net spacers
layer_nom_thickness = 0.0064*1e-3 #m
spacer_nom_thickness = 0.16*1e-3 #m 

k_MLI = 0.24 #W/mK
Cp_MLI = 1170. #J/kgK. Mylar (DAM Cp) 

emissivity_MLI = 0.03 #inner, from Thermal Control handbook
emissivity_coating = 0.66 # outer coating Teflon Aluminium Backing 2 mil
absorptivity_coating = 0.10 # outer coating Teflon Aluminium Backing 2 mil

# Tank wall properties (when applicable)
Cp_wall =  526.3 # Titanium: 526.3 Aluminum 2195-T8: 900. Stainless Steel:500
k_wall =  6.7 # Titanium: 6.7 Aluminum 2195-T8: 130.Stainless Steel:16.3
density_wall =  4430. # Titanium: 4430. Aluminum 2195-T8: 3000. Stainless Steel:8000
t_wall =  1.27*1e-2 # Titanium: 0.005 # Validation: tank_thickness Aluminum&Stainless Steel: 1.27 cm 

# Vapor Cooled Shield parameters
diameter_VCS_tube = 11.7 * 1e-3 
length_VCS_tube = 2*52.19 #m --> set this according to tank dimensions
shield_density = 2660 #kg/m3. aluminum 2660, copper: 8960 
shield_thickness = 0.1*1e-3  #m
shield_Cp = 900 #J/kgK

VCS_timestep = 1e-4 #s

# Tank venting
desired_venting_pressure = 1.3*1e5 #Pa
venting_in_one_step = 'TRUE' #DEACTIVATE if you want venting to last longer thank 1 timestep. 
##### WARNING: the part corresponding to venting_in_one_step ='FALSE' needs to be written in the code!!

# Environment
gravity = 1e-7* 9.81 #m/s2. 1e-6* 9.81 for LEO, 1e-7* 9.81 for GEO
T_space = 3. #K
T_rejection = 273. #K

#Custom designs in format custom_design = [(layer_dens, MLI layers, Q_cooler), (option 2)], etc
# The customed designs will be added to the combinations that have already been chosen above
custom_designs = []