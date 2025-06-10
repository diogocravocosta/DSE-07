#imports 
import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
from tank_sizing import TankSizer
import csv
import pandas as pd
import os

class BluntBody:
    '''
    Class: Defining geometry and calculating aerodynamic characteristics of blunt bodies
    Author: Soham Katewale

    Initialised Variables:
        Defined geometry:
        cone_length = length of the cone
        cone_max_radius = maximum radius of the cone (base of the cone)
        cone_min_radius = cone minimum radius (tip of the cone)
        base_arc_height = arc height of the spherical base
        mass = mass of the vehicle
    
        Calculated geometry:
        sphere_geometry = radius of the sphere used to create the curved base4
        mu_b = half arc angle of the curved base
        alpha_max = maximum angle of attack before numerical integration is required (i.e becomes too complex)
        cone_half_angle = half angle of the cone
    '''
    def __init__(self, 
                 cone_length: float, #Cone length
                 cone_max_radius: float, #Cone maximum radius
                 cone_min_radius: float, #Cone minimumr radius
                 base_arc_height: float, #Arc height of spherical base
                 mass: float #Mass of vehicle
                 ):
        #initialisations
        self.cone_length = cone_length
        self.cone_max_radius = cone_max_radius
        self.cone_min_radius = cone_min_radius
        self.base_arc_height = base_arc_height
        self.mass = mass

        #Calculated geometry 
        self.sphere_radius = (2*self.cone_max_radius)**2 / (8*self.base_arc_height) + self.base_arc_height/2 #Radius of the hemisphere (base spherical base of the blunt body)
        assert self.sphere_radius >= self.cone_max_radius, "Hemisphere radius is too small, sphereRadius must be >= coneMaxRadius"
        self.mu_b = np.acos((self.sphere_radius - self.base_arc_height)/self.sphere_radius) #half arc angle of hemisphere
        self.alpha_max = np.pi/2 - self.mu_b #maximum angle of attack before requiring numerical integration
        self.cone_half_angle = np.arctan(self.cone_max_radius / ((self.cone_max_radius*self.cone_length) / (self.cone_max_radius-self.cone_min_radius))) #cone half angle
        self.cap_centroid = self.sphere_radius
        self.cone_centroid = self.sphere_radius*(1 - np.cos(self.mu_b)) + (self.cone_length - (2*self.cone_length/3*(np.cos(self.cone_half_angle)**2)))
        
    def compute_aerodynamics(self):
        '''
        Function: Computes the aerodynamic coefficients

        Variables:
            c_x_cap = axial force coefficient acting on the cap (spherical base of the blunt body)
            c_y_cap = normal force coefficient acting on the cap
            c_x_cone = axial force coefficient acting on the cone
            c_y_cone = normal force coefficient acting on the cone
            aoa_all = whole range of angle of attacks
            rho_str = angle for tangential incidence at the cone
            c_l = lift coefficients
            c_d = drag coefficients
            lift_over_drag = L/D force ratio
            normal_over_axial = N/A force ratio
        '''
        '''
        calculating axial and normal coefficients
        '''
        self.aoa_all = np.linspace(0.1 * np.pi/180, self.alpha_max , 100)
        self.c_x_cap = 0.5 * (np.sin(self.aoa_all)**2) * (np.sin(self.mu_b)**2) + \
        (1 + np.cos(self.mu_b)**2) * (np.cos(self.aoa_all)**2)
        self.c_y_cap = np.sin(self.aoa_all) * np.cos(self.aoa_all) * (np.sin(self.mu_b)**2)

        self.c_x_cone = []
        self.c_y_cone = []
        for angle in self.aoa_all:
            rho_str = np.arccos(-1 * np.tan(self.cone_half_angle) / np.tan(angle))

            if np.abs(angle) >= self.cone_half_angle:
                c_x_cone_1 = -1 * ((np.sin(angle))**2) * ((np.cos(self.cone_half_angle))**2) / np.pi
                c_x_cone_2 = (1 + 2*(np.cos(rho_str))**2) * (np.pi - rho_str)
                c_x_cone_3 = 3 * np.sin(rho_str) * np.cos(rho_str)

                self.c_x_cone.append(c_x_cone_1 * (c_x_cone_2 + c_x_cone_3))

                c_y_cone_1 = 2 * (((np.sin(angle))**2) * ((np.cos(self.cone_half_angle))**2) / np.pi) * (1 / np.tan(self.cone_half_angle))
                c_y_cone_2 = (np.pi - rho_str) * np.cos(rho_str)
                c_y_cone_3 = (np.sin(rho_str) / 3) * (2 + (np.cos(rho_str))**2)
                self.c_y_cone.append(c_y_cone_1 * (c_y_cone_2 + c_y_cone_3))
                                        
            else: 
                self.c_x_cone.append(0)
                self.c_y_cone.append(0)
        
        self.c_axial = (self.c_x_cap + self.c_x_cone)
        
        '''
        calculating lift and drag coefficients
        '''
        self.c_l = (((self.c_y_cap + self.c_y_cone)*np.cos(self.aoa_all)) - ((self.c_x_cap + self.c_x_cone)*np.sin(self.aoa_all))) * -1
        self.c_d = ((self.c_y_cap + self.c_y_cone)*np.sin(self.aoa_all)) + ((self.c_x_cap + self.c_x_cone)*np.cos(self.aoa_all))
        self.lift_over_drag = self.c_l / self.c_d
        self.normal_over_axial = (self.c_y_cap + self.c_y_cone) / (self.c_x_cap + self.c_x_cone)

        '''
        Calculating cl_alpha
        '''
        cl_max = max(self.c_l)
        cl_min = min(self.c_l)
        cl_max_idx = np.where(self.c_l == cl_max)
        self.cl_alpha_calc = (cl_max - cl_min) / (self.aoa_all[cl_max_idx]) * -1
        

        

        
    def flight_profile(self, velocities, densities, altitudes, angle_of_attack):
        '''
        Function: Plots the lift and drag over the flight profile

        Variables:

        '''

            


    def draw_geometry(self):
        '''
        Function: draws the geometry of the vehicle
        '''

        '''
        Hemisphere mesh grid
        '''
        mu = np.linspace(0, self.mu_b, 100)
        rho = np.linspace(0, 2*np.pi, 100)
        mu, rho = np.meshgrid(mu, rho)

        '''
        Hemisphere parametric equations
        '''
        self.x_sphere = self.sphere_radius * (1 - np.cos(mu))
        self.y_sphere = self.sphere_radius * np.sin(mu) * np.cos(rho)
        z_sphere = self.sphere_radius * np.sin(mu) * np.sin(rho)

        '''
        Cone mesh grid
        '''
        self.x_cone = np.linspace(self.base_arc_height, self.base_arc_height + self.cone_length, 100)
        self.rho_cone = np.linspace(0, 2*np.pi, 100)
        x_cone, rho_cone = np.meshgrid(self.x_cone, self.rho_cone)

        '''
        Cone parametric equations
        '''
        self.r_cone = self.cone_max_radius - ((self.cone_max_radius - self.cone_min_radius)/self.cone_length) * (x_cone - self.base_arc_height)
        self.y_cone = self.r_cone * np.cos(self.rho_cone)
        z_cone = self.r_cone * np.sin(self.rho_cone)

        '''
        Plotting
        '''
        fig = plt.figure(figsize=(10, 7))
        ax = fig.add_subplot(111, projection='3d')

        # Hemisphere (front)
        ax.plot_surface(self.x_sphere, self.y_sphere, z_sphere, color='skyblue', alpha=0.9)

        # Cone (extending backward)
        ax.plot_surface(self.x_cone, self.y_cone, z_cone, color='lightcoral', alpha=0.9)

        # Labels and aspect
        ax.set_xlabel('x')
        ax.set_ylabel('y')
        ax.set_zlabel('z')
        ax.set_title('Blunt Body (Hemisphere + Cone)')
        ax.set_box_aspect([2, 1, 1])  # Wider x to show elongation

        plt.show()


    def RASAero(self, file_name, mode):
        '''
        Function: To assess stability of a geometry
        Input: Requires CSV file output from RASAero 

        Mode 1: Considers the default output from RASAero
        Mode 2: Considers the "Run Test" output from RASAero for a specific angle of attack
        
        Variables:
            cl_alpha = cl-alpha slope
            cmq_cmadot = stability coefficients
            moment_inertia = moment of inertia
            radius_gyration = radius of gyration
        '''
        print("\033[91mVerify that the output file from RASAero is updated for the current geometry\033[0m")
        self.moment_inertia = 231368.0069
        self.radius_gyration = np.sqrt(self.moment_inertia / self.mass)
                                       
        if mode == 1:
            current_dir = os.path.dirname(__file__)
            data_path = os.path.join(current_dir, '..', 'data', file_name)
            data_path = os.path.abspath(data_path)
            print("Resolved path:", data_path)
            data = pd.read_csv(data_path)

            data_alpha_0 = data[data["Alpha"] == 0]
            self.data_alpha_2 = data[data["Alpha"] == 2]
            self.data_alpha_4 = data[data["Alpha"] == 4]
            self.RASAero_mach = np.array(self.data_alpha_4["Mach"])

            self.cl_alpha_rasaero = np.array((self.data_alpha_4["CL"])/ 2)

            self.cmq_cmadot = (4*self.moment_inertia / (self.mass * (self.cone_max_radius*2)**2)) * self.cl_alpha_rasaero
            print(self.cmq_cmadot)

            data_drag_4 = self.data_alpha_4['CD']
            data_lift_4 = self.data_alpha_4["CL"]
            self.cd_4 = np.array(data_drag_4)
            self.cl_4 = np.array(data_lift_4)
        
        if mode == 2:
            self.cmq_cmadot = -1 * (4*self.moment_inertia / (self.mass * (self.cone_max_radius*2)**2)) * self.c_axial
            #self.cmq_cmadot = -1 * (4*self.moment_inertia / (self.mass * (self.sphere_radius*2)**2)) * self.c_d
            self.cmq = (self.cl_alpha_calc - (2*self.c_d)) / (self.mass * (self.sphere_radius*2)**2 / self.moment_inertia)


    
    def analysis(self):
        '''
        Function: Runs a bunch of stuff to visualise the impact of changing certain geometries
        '''

        '''
        Analysing effect of taper ratios 
        '''
        max_ld = []
        max_ld_aoa = []
        taper_ratios = np.arange(0, 0.98, 0.01)
        volume = 619
        for t in taper_ratios:
            ts = TankSizer(volume, taper = t)
            cone_length = ts.h
            cone_max_radius = ts.r_bottom
            cone_min_radius = ts.r_top
            base_arc_height = 2
            mass = 28000

            body = BluntBody(cone_length, cone_max_radius, cone_min_radius, base_arc_height, mass)
            body.compute_aerodynamics()
            max_ld_index = np.argmax(body.lift_over_drag)
            max_ld.append(body.lift_over_drag[max_ld_index])
            max_ld_aoa.append(body.aoa_all[max_ld_index])

        x = taper_ratios
        y = np.array(max_ld_aoa) * 180/np.pi
        z = np.array(max_ld)

        fig = plt.figure()
        ax = fig.add_subplot(projection='3d')
        ax.scatter(x, y, z, c=z, cmap="viridis")

        ax.set_xlabel('Taper Ratios')
        ax.set_ylabel('Angle of attack')
        ax.set_zlabel('L/D')

        plt.show()

        '''
        Analysing effect of sphericity of base
        '''
        max_ld = []
        max_ld_aoa = []
        sphericity_radius = np.linspace(0.1, cone_max_radius, 100)
        for s in sphericity_radius:
            body = BluntBody(cone_length = 15, cone_max_radius = 5, cone_min_radius = 1, base_arc_height = s, mass = 28000)
            body.compute_aerodynamics()
            max_ld_index = np.argmax(body.lift_over_drag)
            max_ld.append(body.lift_over_drag[max_ld_index])
            max_ld_aoa.append(body.aoa_all[max_ld_index])

        x = sphericity_radius
        y = np.array(max_ld_aoa) * 180/np.pi
        z = np.array(max_ld)

        fig = plt.figure()
        ax = fig.add_subplot(projection='3d')
        ax.scatter(x, y, z, c=z, cmap="viridis")

        ax.set_xlabel('Base arc height')
        ax.set_ylabel('Angle of attack')
        ax.set_zlabel('L/D')

        plt.show()

    def plots(self):
        '''
        Aerodynamic coefficient plots
        '''
        # plt.plot(self.aoa_all * 180/np.pi, self.c_y_cap + self.c_y_cone, label = "normal force coeff", color = "blue")
        # plt.plot(self.aoa_all * 180/np.pi, self.c_x_cap + self.c_x_cone, label = "axial force coeff", color = "red")
        # plt.plot(self.aoa_all * 180/np.pi, self.normal_over_axial, label = "normal / axial", color = "green")
        plt.figure()
        plt.plot(self.aoa_all * 180/np.pi, self.lift_over_drag, label="lift over drag", color = 'red')
        plt.plot(self.aoa_all * 180/np.pi, self.c_l, label="lift coefficient", color = 'green')
        plt.plot(self.aoa_all * 180/np.pi, self.c_d, label="drag coefficient", color = 'blue')
        plt.legend()

        # plt.figure()
        # plt.plot(self.RASAero_mach, self.cd_4, label='C_d vs Mach at a = 4 / RASAero')
        # plt.legend()

        # plt.figure()
        # plt.plot(self.RASAero_mach, self.cl_4, label='C_l vs Mach at a = 4 / RASAero')
        # plt.legend()

        '''
        Stability plots
        '''
        plt.figure()
        plt.plot(self.aoa_all * 180/np.pi, self.cmq_cmadot, label='stability')
        plt.legend()
        plt.show()

if __name__ == "__main__":
    #CURRENT DIMENSIONS
    cone_length = 13.95
    cone_max_radius = 4.95
    cone_min_radius = 2.46
    base_arc_height = 2
    mass = 28000


    # #APOLLO CAPSULE
    # # cone_length = 2.662
    # # cone_max_radius = 1.956
    # # cone_min_radius = 0.01
    # # base_arc_height = 0.5
    # # mass = 28000

    body = BluntBody(cone_length, cone_max_radius, cone_min_radius, base_arc_height, mass)

    #body.draw_geometry()

    body.compute_aerodynamics()

    body.RASAero(mode = 2, file_name = "HermesV1-RASAero.csv")

    #body.analysis()

    body.plots()





                
                

