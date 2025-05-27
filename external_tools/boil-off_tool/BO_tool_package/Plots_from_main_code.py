
from turtle import left
import numpy as np
import matplotlib.pyplot as plt
from pathlib import Path
import os

########## Please manually set this value #####
VD_MLI_activity = 'TRUE'
CASE = 4
#####################################################

BASE_PATH = Path.cwd()/ 'Output_data' 

CB_color_cycle = ['#377eb8', '#ff7f00', '#4daf4a',
                  '#f781bf', '#a65628', '#984ea3',
                  '#999999', '#e41a1c', '#dede00']

# import tank geometry from code for heat flux calculation
with open(os.path.join(BASE_PATH,'Tank_geometry.npz'),'rb') as f:
   TankGeom = np.load(f, allow_pickle = 'True')
   tank_endcap = TankGeom['tank_endcap']
   endcap_parts = TankGeom['sphere_parts']
   cylinder_parts = TankGeom['cylinder_parts']
   radius = TankGeom['radius']
   height = TankGeom['height']
   surf_area_cap_node = TankGeom['surf_area_cap_node']
   surf_area_cyl_node = TankGeom['surf_area_cyl_node']

surf_area_caps = surf_area_cap_node*endcap_parts
surf_area_cyl = surf_area_cyl_node*cylinder_parts
total_surf_area = surf_area_caps + surf_area_cyl
print('tot surface area is:',total_surf_area)


#read data from simulation
if CASE == 1:
    with open(os.path.join(BASE_PATH,'PlotDataR.npz'),'rb') as f:

        PlotDataR = np.load(f)
        T_VD_MLI_R = PlotDataR['T_VD_MLI_R']
        T_matrix_R = PlotDataR['T_matrix_R']
        thick_vector = PlotDataR['thickness_vector']
        time_vect = PlotDataR['time_vect']


    with open(os.path.join(BASE_PATH,'PlotDataL.npz'),'rb') as f:

        PlotDataL = np.load(f)
        T_VD_MLI_L = PlotDataL['T_VD_MLI_L']
        T_matrix_L = PlotDataL['T_matrix_L']

    with open(os.path.join(BASE_PATH,'dPData.npz'),'rb') as f:

        dPData = np.load(f)
        dP_vapour_vector = dPData['dP_vapour_vector']
        dP_liquid_vector = dPData['dP_liquid_vector']
        dP_radiative_net = dPData['dP_radiative_net']
        dP_emitted_back = dPData['dP_emitted_back'] 

    with open(os.path.join(BASE_PATH,'VentData.npz'),'rb') as f:

        VentData = np.load(f)
        p_vapour_vector = VentData ['p_vapour_vector']
        m_vapour_vector = VentData ['m_vapour_vector'] 
        T_vapour_vect = VentData ['T_vapour_vect'] 
        T_liquid_vect = VentData ['T_liquid_vect'] 
        dm_vent_vector = VentData ['dm_vent_vector'] 
        m_boiloff_vector = VentData ['m_boiloff_vector']
        T_VCS_vect = VentData['T_VCS_vect']

    if VD_MLI_activity == 'TRUE':
        # Define vector for thickness (cumulative)
        t_wall = thick_vector[0]
        print('wall thickness',t_wall)
        t_wall = 0.01

        t_MLI_inner = thick_vector[1]
        t_MLI_middle = thick_vector[2]
        t_MLI_outer = thick_vector[3]
        thickness = 100*np.array([ t_wall + t_MLI_inner/2, t_wall + t_MLI_inner+t_MLI_middle/2, t_wall + t_MLI_inner+t_MLI_middle+t_MLI_outer/2]) 
        thickness_wall = 100*np.array([0,t_wall + t_MLI_inner/2, t_wall + t_MLI_inner+t_MLI_middle/2, t_wall + t_MLI_inner+t_MLI_middle+t_MLI_outer/2]) 
        MLI_layers_thickness = 100*np.array([ t_wall, t_wall + t_MLI_inner, t_wall + t_MLI_inner+t_MLI_middle, t_wall + t_MLI_inner+t_MLI_middle+t_MLI_outer])
        # in cm

        plot1 = plt.figure(1)
        
        plt.plot(thickness_wall, np.append(T_matrix_R[0,1],np.flip(T_VD_MLI_R[0,:])),'d-', color = CB_color_cycle[0])
        plt.plot(thickness_wall, np.append(T_matrix_R[1,1],np.flip(T_VD_MLI_R[1,:])),'d-', color = CB_color_cycle[1])
        plt.plot(thickness_wall, np.append(T_matrix_R[2,1],np.flip(T_VD_MLI_R[2,:])),'d-', color = CB_color_cycle[2])
        plt.plot(thickness_wall, np.append(T_matrix_R[3,1],np.flip(T_VD_MLI_R[3,:])),'d-', color = CB_color_cycle[3])
        plt.plot(thickness_wall, np.append(T_matrix_R[4,1],np.flip(T_VD_MLI_R[4,:])),'d-', color = CB_color_cycle[4])
        plt.plot(thickness_wall, np.append(T_matrix_R[5,1],np.flip(T_VD_MLI_R[5,:])),'d-', color = CB_color_cycle[5])
        
        plt.legend(['Node 1 (endcap top)','Node 2','Node 3','Node 4','Node 5', 'Node 6 (endcap bottom)'], loc="best", bbox_transform=plot1.transFigure, ncol=2)
        plt.title('Tank shadow side nodes temperature distribution')
        plt.xlabel('Thickness [cm]')
        plt.ylabel('Temperature [K]')

        plot2 = plt.figure(2)
        
        plt.plot(thickness_wall, np.append(T_matrix_L[0,1], np.flip(T_VD_MLI_L[0,:])),'d-',color = CB_color_cycle[0])
        plt.plot(thickness_wall, np.append(T_matrix_L[1,1],np.flip(T_VD_MLI_L[1,:])),'d-', color = CB_color_cycle[1])
        plt.plot(thickness_wall, np.append(T_matrix_L[2,1],np.flip(T_VD_MLI_L[2,:])),'d-', color = CB_color_cycle[2])
        plt.plot(thickness_wall, np.append(T_matrix_L[3,1],np.flip(T_VD_MLI_L[3,:])),'d-', color = CB_color_cycle[3])
        plt.plot(thickness_wall, np.append(T_matrix_L[4,1],np.flip(T_VD_MLI_L[4,:])),'d-', color = CB_color_cycle[4])
        plt.plot(thickness_wall, np.append(T_matrix_L[5,1],np.flip(T_VD_MLI_L[5,:])),'d-', color = CB_color_cycle[5])
        plt.vlines(MLI_layers_thickness, 0, 225, colors='k', linestyles='dotted')
        plt.vlines(0, 0, 225, colors='k', linestyles='solid')
        plt.vlines(t_wall*100, 0, 225, colors='k', linestyles='solid')

        plt.title('Tank sunlit side nodes temperature distribution')
        plt.xlabel('Thickness [cm]')
        plt.ylabel('Temperature [K]')
        plt.grid()



    fig, ax1 = plt.subplots() 
    
    ax1.set_xlabel('Time [days]') 
    ax1.set_ylabel('Vapors mass [kg]', color = 'red') 
    ax1.plot(time_vect/(60*60*24), dm_vent_vector, color = 'red')
    ax1.tick_params(axis ='y', labelcolor = 'red')  
    # Adding Twin Axes
    ax2 = ax1.twinx() 
    ax2.set_ylabel('Vapors mass [g]', color = 'black') 
    ax2.plot(time_vect/(60*60*24), m_boiloff_vector*1e3, color = 'black') #marker="o",
    ax2.tick_params(axis ='y', labelcolor = 'black') 
    ax1.grid()

    fig, (ax1, ax2) = plt.subplots(2)

    ax1.plot(time_vect/(60*60*24), p_vapour_vector*1e-5,color = CB_color_cycle[0])
    ax1.set_xlabel('Time [days]',fontsize=14)
    ax1.set_ylabel('Ullage pressure [bar]',fontsize=14)
    ax1.tick_params(axis='both', which='major', labelsize=14)
    ax1.grid()

    ax2.plot(time_vect/(60*60*24), T_liquid_vect,color = CB_color_cycle[5])
    ax2.plot(time_vect/(60*60*24), T_vapour_vect,color = CB_color_cycle[1])
    plt.plot(time_vect/(60*60*24), T_VCS_vect,color = CB_color_cycle[2]) 
    #plt.xlim(left = 20)
    ax2.set_xlabel('Time [days]',fontsize=14)
    ax2.set_ylabel('Temperature [K]',fontsize=14)
    ax2.tick_params(axis='both', which='major', labelsize=14)
    ax2.legend(["Liquid","Vapor","VCS"],loc="best",fontsize=14)#,"T VCS"]) #
    ax2.grid()

    plot6 = plt.figure(6)
    plt.plot(time_vect/(60*60*24), p_vapour_vector*1e-5,color = CB_color_cycle[0])
    plt.xlabel('Time [days]',fontsize=14)
    plt.ylabel('Ullage pressure [bar]',fontsize=14)
    plt.tick_params(axis='both', which='major', labelsize=14)
    #plt.ylim(bottom = 0)
    plt.grid()

    plot5 = plt.figure(5)
    plt.plot(time_vect/(60*60*24), dP_vapour_vector,color = CB_color_cycle[0])
    plt.plot(time_vect/(60*60*24), dP_liquid_vector,color = CB_color_cycle[1])
    plt.legend(["dP vapour net","dP liquid net"])
    plt.xlabel('Time [days]')
    plt.ylabel('dP net (all nodes)[W]')
    plt.grid()

    plot8 = plt.figure(8)
    plt.plot(time_vect/(60*60*24), m_vapour_vector)
    plt.xlabel('Time [days]', fontsize = 14)
    plt.ylabel('Vapor mass in ullage [kg]', fontsize = 14)
    plt.tick_params(axis='both', which='major', labelsize=14)
    plt.grid()


elif CASE == 2:
    with open(os.path.join(BASE_PATH,'T_and_dP_Data-1m-option1.npz'),'rb') as f:

        TandPData = np.load(f)
        T_matrix_R = TandPData['T_matrix_R']
        T_matrix_L = TandPData['T_matrix_L']
        dP_net_matrix_R = TandPData['dP_net_matrix_R']
        dP_net_matrix_L = TandPData['dP_net_matrix_L']
 
    plot20 = plt.figure(20)
    plt.plot(['A','i'],T_matrix_R[0,:], 'o-') #['A','i','B','C','D']
    plt.plot(['A','i'],T_matrix_R[1,:], 'o-')
    # plt.plot(['A','i'],T_matrix_R[2,:], 'o-')
    # plt.plot(['A','i'],T_matrix_R[3,:], 'o-')
    # plt.plot(['A','i'],T_matrix_R[4,:], 'o-')
    # plt.plot(['A','i'],T_matrix_R[5,:], 'o-')
    # plt.plot(['A','A','A','A','A','A'], T_hot[:,0],'x')
    plt.xlabel('Subnode [-]')
    plt.ylabel('Temperature [K]')
    plt.legend(["1","2","3","4","5","6"])
    plt.title("Final temperatures right side of tank")
    plt.grid()


    plot21 = plt.figure(21)
    # plt.plot(['A','i'],T_matrix_L[0,:], 'o')
    plt.plot(['A','i'],T_matrix_L[4,:], 'o')
    # plt.plot(['A','i'],T_matrix_L[2,:], 'o')
    # plt.plot(['A','i'],T_matrix_L[3,:], 'o')
    # plt.plot(['A','i'],T_matrix_L[4,:], 'o')
    # plt.plot(['A','i'],T_matrix_L[5,:], 'o')
    # plt.plot(['A','A','A','A','A','A'], T_hot[:,1],'x')
    plt.xlabel('Subnode [-]')
    plt.ylabel('Temperature [K]')
    # plt.legend(["1","2","3","4","5","6"])
    plt.title("Final temperatures left side of tank")
    plt.grid()

    plot22 = plt.figure(22)
    # plt.plot(['Out','Mid', 'i'],dP_VD_MLI_R[0,:], 'o-')
    plt.plot(['A','i'],dP_net_matrix_L[1,:], 'o-')
    # plt.plot(['A','i'],dP_net_matrix_L[1,:], 'o-')
    # # plt.plot(['A','i'],dP_net_matrix_R[3,:], 'o-')
    # # plt.plot(['A','i'],dP_net_matrix_R[4,:], 'o-')
    # # plt.plot(['A','i'],dP_net_matrix_R[5,:], 'o-')
    plt.xlabel('Subnode [-]')
    plt.ylabel('Heat load [W]')
    # plt.legend(["1","2","3","4","5","6"])
    plt.title("Final heat loads left side of tank")
    plt.grid()

    # print(dP_in_tank_chai)

elif CASE == 3:
    with open(os.path.join(BASE_PATH,'PlotDataR.npz'),'rb') as f:

        PlotDataR = np.load(f)
        T_VD_MLI_R = PlotDataR['T_VD_MLI_R']
        T_matrix_R = PlotDataR['T_matrix_R']
        thick_vector = PlotDataR['thickness_vector']
        time_vect = PlotDataR['time_vect']

    with open(os.path.join(BASE_PATH,'VentData.npz'),'rb') as f:

        VentData = np.load(f)
        p_vapour_vector = VentData ['p_vapour_vector']
        T_vapour_vect = VentData ['T_vapour_vect'] 
        T_liquid_vect = VentData ['T_liquid_vect']
        dm_vent_vector = VentData ['dm_vent_vector'] 
        m_boiloff_vector = VentData ['m_boiloff_vector']


    fig, ax1 = plt.subplots() 
    
    ax1.set_xlabel('Time [days]') 
    ax1.set_ylabel('Vapors mass [kg]', color = 'red') 
    ax1.plot(time_vect/(60*60*24), dm_vent_vector, color = 'red')
    ax1.tick_params(axis ='y', labelcolor = 'red')  
    # Adding Twin Axes
    ax2 = ax1.twinx() 
    ax2.set_ylabel('Vapors mass [g]', color = 'black') 
    ax2.plot(time_vect/(60*60*24), m_boiloff_vector*1e3,) #marker="o", color = 'black'
    ax2.tick_params(axis ='y', labelcolor = 'black') 
    ax1.grid()

    plot40 = plt.figure(40)
    plt.plot(time_vect/(60*60*24), p_vapour_vector*1e-5)
    plt.xlabel('Time [days]')
    plt.ylabel('P vapour [bar]')
    plt.grid()
   
    plot60 = plt.figure(60)
    plt.plot(time_vect/(60*60*24), T_liquid_vect,color = CB_color_cycle[0])
    plt.plot(time_vect/(60*60*24), T_vapour_vect,color = CB_color_cycle[1])
    plt.xlabel('Time [days]')
    plt.ylabel('Temperature [K]')
    plt.legend(["T liquid","T_vapour"]) #
    plt.grid()


elif CASE == 4:

    with open(os.path.join(BASE_PATH,'VentData.npz'),'rb') as f:

        VentData = np.load(f)
        p_vapour_vector = VentData ['p_vapour_vector']
        m_vapour_vector = VentData ['m_vapour_vector'] 
        T_vapour_vect = VentData ['T_vapour_vect'] 
        T_liquid_vect = VentData ['T_liquid_vect'] 
        dm_vent_vector = VentData ['dm_vent_vector'] 
        m_boiloff_vector = VentData ['m_boiloff_vector']
        
    with open(os.path.join(BASE_PATH,'PlotDataR.npz'),'rb') as f:

        PlotDataR = np.load(f)
        thick_vector = PlotDataR['thickness_vector']
        time_vect = PlotDataR['time_vect']


    with open(os.path.join(BASE_PATH,'MLIcheck.npz'),'rb') as f:

        MLIcheck = np.load(f)
        T_matrix_L = MLIcheck['T_matrix_L']
        T_liquid_vect = MLIcheck['T_liquid_vect']
        T_VCS_vect = MLIcheck['T_VCS_vect'] 
        


    plot2 = plt.figure(2)
    plt.plot(time_vect/(60*60*24), T_VCS_vect,color = CB_color_cycle[0])
    plt.xlabel('Time [days]')
    plt.ylabel('Temperature [K]')
    plt.legend(["VCS","fluid out", "fluid in"]) #
    plt.grid()



    fig, (ax1, ax2) = plt.subplots(2)
    #fig.suptitle('Vertically stacked subplots')
    ax1.plot(time_vect/(60*60*24), p_vapour_vector*1e-5,color = CB_color_cycle[0])
    ax1.set_xlabel('Time [days]',fontsize=14)
    ax1.set_ylabel('Ullage pressure [bar]',fontsize=14)
    ax1.tick_params(axis='both', which='major', labelsize=14)
    ax1.grid()
    


    #plot6 = plt.figure(6)
    ax2.plot(time_vect/(60*60*24), T_liquid_vect,color = CB_color_cycle[5])
    ax2.plot(time_vect/(60*60*24), T_vapour_vect,color = CB_color_cycle[1])
    ax2.plot(time_vect/(60*60*24), T_VCS_vect,color = CB_color_cycle[2])
    ax2.set_xlabel('Time [days]',fontsize=14)
    ax2.set_ylabel('Temperature [K]',fontsize=14)
    ax2.tick_params(axis='both', which='major', labelsize=14)
    ax2.legend(["Liquid","Vapor","VCS"],loc="best",fontsize=14)#,"T VCS"])
    # ax2.legend(["VCS shield"],loc="best",fontsize=14)
    ax2.grid()
    
plt.show() 
