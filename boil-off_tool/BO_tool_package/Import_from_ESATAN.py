
import time
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import os
from pathlib import Path
from math import *

##### Some notes about the use of this code:############################################################################à
# the code takes the data in *.txt format from the radiative analysis on ESATAN
# this code is MODEL-specific, meaning that the loop has been set according to the order in which the data related to each model face is reported from ESATAN
# in particular, in this case, we start from node 4 (bottom cylinder) on the right side, and go to the top taking alternatively nodes from R and L 
# in the end, endcaps are done, first bottom part (index 5), then top part (index 0)
# data are stored in two matrices: total fluxes on right side and total fluxes on left side. Each matrix has n rows= time steps from ESATAN and n columns=nodes on L or R
# before importing data, all the NON RELEVANT faces (like solar panels and inner surfaces) have been deleted from the *.txt file
####################################################################################################################################


############# PARAMETERS TO CHANGE EVERY TIME A NEW GEOMETRY OR DIFFERENT TIME STEPS ARE IMPORTED ##########################à
REPORT_NAME = 'GEO_rad_report_fluxes_single.txt' 
BASE_PATH = Path.cwd()/ 'Radiative_reports_from_ESATAN' #check if path is correct
SAVE_WITH_NAME =  'Radiative_model_GEO_single' #indicate the name of the .npz radiative output file you want to generate
positions =  20 # change this according to orbit positions set in ESATAN 
store_data = 'YES SOLAR' # when ready to store data on .npz file, set this to 'YES ALL' or 'YES SOLAR' (only for solar fluxes)
#############################################################################################################

rows_number = positions + 4 #(positions + 4)
n_skip = rows_number + 3 # (rows_number + 3)

node_index = 4
right_side = 'False' 
endcaps = 'False'

stay_in_loop = 'True'
first_loop = 'True'

skipped_rows = 0
counter = 0

while stay_in_loop == 'True':
   
    X = pd.read_csv(os.path.join(BASE_PATH,REPORT_NAME), delim_whitespace = True, header=None,skiprows=skipped_rows,nrows=rows_number)
    

    Data_Array = X.to_numpy(dtype=None, copy=False) #convert dataframe to numpy array
    time_steps = Data_Array[2:-1,1]
    abs_solar = Data_Array[2:-1,5]
    abs_albedo = Data_Array[2:-1,6]
    abs_planetary = Data_Array[2:-1,7]
    inc_solar = Data_Array[2:-1,2]
    inc_albedo = Data_Array[2:-1,3]
    inc_planetary = Data_Array[2:-1,4]

    if first_loop == 'True':
        abs_flux_tot_R = np.zeros((len(time_steps),6))
        abs_flux_tot_L = np.zeros((len(time_steps),6))
        abs_flux_S_R = np.zeros((len(time_steps),6))
        abs_flux_S_L = np.zeros((len(time_steps),6))
        inc_flux_ISIA_R = np.zeros((len(time_steps),6))
        inc_flux_ISIA_L = np.zeros((len(time_steps),6))
        inc_flux_IP_R = np.zeros((len(time_steps),6))
        inc_flux_IP_L = np.zeros((len(time_steps),6))
        inc_flux_IA_R = np.zeros((len(time_steps),6))
        inc_flux_IA_L = np.zeros((len(time_steps),6))
        inc_flux_IS_R = np.zeros((len(time_steps),6))
        inc_flux_IS_L = np.zeros((len(time_steps),6))
    
    for i in range(0, len(time_steps)):
        time_steps[i] = float(time_steps[i])
        abs_solar[i] = float(abs_solar[i])
        abs_albedo[i] = float(abs_albedo[i]) 
        abs_planetary[i] = float(abs_planetary[i])
        inc_solar [i] = float(inc_solar[i])
        inc_albedo [i] = float(inc_albedo[i])
        inc_planetary [i]= float(inc_planetary[i])

    abs_tot = abs_solar + abs_albedo + abs_planetary
    inc_ISIA = inc_solar + inc_albedo
    
    if right_side == 'False':
        abs_flux_tot_R[:,node_index] = abs_tot
        abs_flux_S_R[:,node_index] = abs_solar
        inc_flux_ISIA_R[:,node_index] = inc_ISIA
        inc_flux_IP_R[:,node_index] = inc_planetary
        inc_flux_IA_R[:,node_index] = inc_albedo
        inc_flux_IS_R [:,node_index] = inc_solar
        right_side = 'True'

    elif right_side == 'True':
        abs_flux_tot_L[:,node_index] = abs_tot
        abs_flux_S_L[:,node_index] = abs_solar
        inc_flux_ISIA_L[:,node_index] = inc_ISIA
        inc_flux_IP_L[:,node_index] = inc_planetary
        inc_flux_IA_L[:,node_index] = inc_albedo
        inc_flux_IS_L [:,node_index] = inc_solar
        right_side = 'False'
        node_index +=-1
    


    if node_index == 0 and endcaps == 'False': 
        node_index = 5
        endcaps = 'True'
    elif node_index == 4 and endcaps == 'True':
        node_index = 0
    elif node_index == 0 and endcaps == 'True':
        counter += 1
    elif node_index == -1 and endcaps == 'True':
        counter += 1

    if counter == 2:
        stay_in_loop = 'False'
    
    if first_loop == 'True':
        
        skipped_rows = n_skip #79 #43 #17 (rows_number + 3)
        first_loop = 'False'    
    else:
        skipped_rows = skipped_rows + n_skip + 1 #80 #18 #44 (rows_number + 4)


############# STORE data in npz file to export them to main code
if store_data == 'YES ALL':
    np.savez(SAVE_WITH_NAME, time_steps = time_steps, abs_flux_tot_L = abs_flux_tot_L, abs_flux_tot_R = abs_flux_tot_R, inc_flux_ISIA_R = inc_flux_ISIA_R, inc_flux_ISIA_L = inc_flux_ISIA_L, inc_flux_IP_R = inc_flux_IP_R, inc_flux_IP_L = inc_flux_IP_L)
elif store_data == 'YES SOLAR':
    np.savez(SAVE_WITH_NAME, time_steps = time_steps, abs_flux_tot_L = abs_flux_S_L, abs_flux_tot_R = abs_flux_S_R, inc_flux_ISIA_R = inc_flux_IS_R, inc_flux_ISIA_L = inc_flux_IS_L, inc_flux_IP_R = 0*inc_flux_IP_R, inc_flux_IP_L = 0*inc_flux_IP_L)


    
# print('Total abs fluxes right side')
# print(abs_flux_tot_R) 
# print('abs flux tot rows are',len(abs_flux_tot_R))

#print('Total abs fluxes left side')
#print(abs_flux_tot_L)

#print('Total inc fluxes right side')
#print(inc_flux_ISIA_R) 

# print('Total inc fluxes left side')
# print(inc_flux_ISIA_L)

# print('Total inc planet fluxes right side')
# print(inc_flux_IP_R) 

# print('Total inc planet fluxes left side')
# print(inc_flux_IP_L)


### try here to implement time variable fluxes
t = 0
nn = 0 
step = 0 
dt = 10
runtime = 60*60*24*2 


dP_radiative_matrix = np.zeros((3,int(10*runtime)))
time_vect = np.zeros(int(10*runtime))
TIME_VECT = time_steps
one_period = time_steps

period = time_steps[-1]

while t <= runtime:
    t = t + dt
    time_vect[nn] = t
    if step == (len(time_steps)-1) and t > time_steps[step]:
        step = 0
        time_steps = time_steps + period
        TIME_VECT = np.append(TIME_VECT,time_steps)
    elif t >= time_steps[step+1]:
        step +=1

    dP_radiative_matrix [0,nn] = inc_flux_ISIA_L[step,0] + inc_flux_IP_R[step,0] #just to try for node 1
    dP_radiative_matrix [1,nn] = inc_flux_ISIA_L[step,2] + inc_flux_IP_R[step,2] # central node of cylinder
    dP_radiative_matrix [2,nn] = inc_flux_ISIA_L[step,5] + inc_flux_IP_R[step,5]             
    nn +=1

dP_radiative_matrix = dP_radiative_matrix[:,0:nn]
time_vect = time_vect[0:nn]


inc_flux_IS_R_cyl = inc_flux_IS_R [:,2]
inc_flux_IS_L_cyl = inc_flux_IS_L [:,2]

##########################################################

CB_color_cycle = ['#377eb8', '#ff7f00', '#4daf4a',
                  '#f781bf', '#a65628', '#984ea3',
                  '#999999', '#e41a1c', '#dede00']

##################################### Radiative model ###############################################################
cyl_sections = 8
sphere_parts = 4
radius = 2.7
height = 20.

increment_frontalarea = 2.*radius*height/(cyl_sections/2.)
increment_area = pi*radius*height/(cyl_sections/2.)

endcap_increment_area = 4.*pi*radius**2./sphere_parts
endcap_increment_frontalarea = pi*(radius)**2./(0.5*sphere_parts) #every spherical part has an area of half circle
 
absorptivity = 0.10
emissivity = 0.66 

# Solar heat loads, left side of tank

plot1 = plt.figure(1)

plt.plot(one_period/(60*60),abs_flux_S_L[:,0]*endcap_increment_area, linestyle='solid', color = CB_color_cycle[0])
plt.plot(one_period/(60*60),abs_flux_S_L[:,1]*increment_area, linestyle='solid', color = CB_color_cycle[1])
plt.plot(one_period/(60*60),abs_flux_S_L[:,2]*increment_area, linestyle='solid', color = CB_color_cycle[2])
plt.plot(one_period/(60*60),abs_flux_S_L[:,3]*increment_area, linestyle='solid', color = CB_color_cycle[3])
plt.plot(one_period/(60*60),abs_flux_S_L[:,4]*increment_area, linestyle='solid', color = CB_color_cycle[4])
plt.plot(one_period/(60*60),abs_flux_S_L[:,5]*endcap_increment_area, linestyle='solid', color = CB_color_cycle[5])

plt.xlabel('Time [hours]')
plt.ylabel('Heat load [W]')
plt.legend(["Endcap top","Node 2", "Node 3", "Node 4", "Node 5","Endcap bottom"])
plt.title('Solar heat loads on left side of the tank')
plt.grid()

# Solar heat loads, right side of tank

plot2 = plt.figure(2)
plt.plot(one_period/(60*60),abs_flux_S_R[:,0]*endcap_increment_area, linestyle='solid', color = CB_color_cycle[0])
plt.plot(one_period/(60*60),abs_flux_S_R[:,1]*increment_area, linestyle='solid', color = CB_color_cycle[1])
plt.plot(one_period/(60*60),abs_flux_S_R[:,2]*increment_area, linestyle='solid', color = CB_color_cycle[2])
plt.plot(one_period/(60*60),abs_flux_S_R[:,3]*increment_area, linestyle='solid', color = CB_color_cycle[3])
plt.plot(one_period/(60*60),abs_flux_S_R[:,4]*increment_area, linestyle='solid', color = CB_color_cycle[4])
plt.plot(one_period/(60*60),abs_flux_S_R[:,5]*endcap_increment_area, linestyle='solid', color = CB_color_cycle[5])

plt.xlabel('Time [hours]')
plt.ylabel('Heat load [W]')
plt.legend(["Endcap top","Node 2", "Node 3", "Node 4", "Node 5","Endcap bottom"])
plt.title('Solar heat loads on right side of the tank')
plt.grid()

# ALL heat loads, left side of tank

plot3 = plt.figure(3)
plt.plot(one_period/(60*60),abs_flux_tot_L[:,0]*endcap_increment_area, linestyle='solid', color = CB_color_cycle[0])
plt.plot(one_period/(60*60),abs_flux_tot_L[:,1]*increment_area, linestyle='solid', color = CB_color_cycle[1])
plt.plot(one_period/(60*60),abs_flux_tot_L[:,2]*increment_area, linestyle='solid', color = CB_color_cycle[2])
plt.plot(one_period/(60*60),abs_flux_tot_L[:,3]*increment_area, linestyle='solid', color = CB_color_cycle[3])
plt.plot(one_period/(60*60),abs_flux_tot_L[:,4]*increment_area, linestyle='solid', color = CB_color_cycle[4])
plt.plot(one_period/(60*60),abs_flux_tot_L[:,5]*endcap_increment_area, linestyle='solid', color = CB_color_cycle[5])

plt.xlabel('Time [hours]')
plt.ylabel('Heat load [W]')
plt.legend(["Endcap top","Node 2", "Node 3", "Node 4", "Node 5","Endcap bottom"])
plt.title('All heat loads on left side of the tank')
plt.grid()

# ALL heat loads, right side of tank
plot4 = plt.figure(4)
plt.plot(one_period/(60*60),abs_flux_tot_R[:,0]*endcap_increment_area, linestyle='solid', color = CB_color_cycle[0])
plt.plot(one_period/(60*60),abs_flux_tot_R[:,1]*increment_area, linestyle='solid', color = CB_color_cycle[1])
plt.plot(one_period/(60*60),abs_flux_tot_R[:,2]*increment_area, linestyle='solid', color = CB_color_cycle[2])
plt.plot(one_period/(60*60),abs_flux_tot_R[:,3]*increment_area, linestyle='solid', color = CB_color_cycle[3])
plt.plot(one_period/(60*60),abs_flux_tot_R[:,4]*increment_area, linestyle='solid', color = CB_color_cycle[4])
plt.plot(one_period/(60*60),abs_flux_tot_R[:,5]*endcap_increment_area, linestyle='solid', color = CB_color_cycle[5])

plt.xlabel('Time [hours]')
plt.ylabel('Heat load [W]')
plt.legend(["Endcap top","Node 2", "Node 3", "Node 4", "Node 5","Endcap bottom"])
plt.title('All heat loads on right side of the tank')
plt.grid()

plt.show()