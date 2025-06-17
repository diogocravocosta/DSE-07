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
        self.reference_area = np.pi * self.sphere_radius**2
        
    def hypersonic_aerodynamics(self):
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
        self.aoa_all = np.linspace(0.1 * np.pi / 180, self.alpha_max, 100)
        self.aoa_all = np.linspace(0.1 * np.pi / 180, 30 * np.pi / 180, 100)
        #print('alpha max', self.alpha_max)
        #print("angle", self.aoa_all)
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
        
        self.c_axial = self.c_x_cap + self.c_x_cone
        self.c_normal = self.c_y_cap + self.c_y_cone
        
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

    def subsonic_aerodynamics(self, file_name_base, file_name_nose):
        print("\033[91mVerify that the output file from RASAero is updated for the current geometry\033[0m")
        self.radius_gyration = np.sqrt(self.moment_inertia / self.mass)
                                       
        current_dir = os.path.dirname(__file__)
        data_path = os.path.join(current_dir, '..', 'data', file_name_base)
        data_path = os.path.abspath(data_path)
        print("Resolved path:", data_path)
        data_base = pd.read_csv(data_path)

        current_dir = os.path.dirname(__file__)
        data_path = os.path.join(current_dir, '..', 'data', file_name_nose)
        data_path = os.path.abspath(data_path)
        print("Resolved path:", data_path)
        data_nose = pd.read_csv(data_path)

        self.data_alpha_0 = data_base[data_base["Alpha"] == 0]
        self.data_alpha_2 = data_base[data_base["Alpha"] == 2]
        self.data_alpha_4 = data_base[data_base["Alpha"] == 4]

        self.cl_alpha_rasaero = np.array((self.data_alpha_4["CL"])/ 2)

        reentry_data = self.data_alpha_4[self.data_alpha_4["Mach"] < 0.9]
        self.mach_subsonic_reentry = np.array(reentry_data["Mach"])
        self.cl_subsonic_reentry = np.array(reentry_data["CL"])
        self.cd_subsonic_reentry = np.array(reentry_data["CD"])

        self.data_alpha_4_nose = data_nose[data_nose["Alpha"] == 4]
        liftoff_data = self.data_alpha_4_nose[self.data_alpha_4["Mach"] < 0.9]
        self.mach_subsonic_liftoff = np.array(liftoff_data["Mach"])
        self.cl_subsonic_liftoff = np.array(liftoff_data["CL"])
        self.cd_subsonic_liftoff = np.array(liftoff_data["CD"])

        

    def stability(self, velocity, density, deviation, moment_inertia, aoa, xcg):
        '''
        Function: To assess stability of a geometry
        Input: Requires CSV file output from RASAero and now a bunch more stuff
        
        Variables:
            cl_alpha = cl-alpha slope
            cmq_cmadot = stability coefficients
            moment_inertia = moment of inertia
            radius_gyration = radius of gyration

            more variable that should be added
        '''
        self.moment_inertia = moment_inertia
        self.cmq_cmadot = -1 * (4*self.moment_inertia / (self.mass * (self.cone_max_radius*2)**2)) * self.c_axial

    
        '''
        Force calculations Hypersonic
        '''
        self.dynamic_pressure = 0.5 * density * velocity**2
        self.cap_force = self.dynamic_pressure * self.reference_area * self.c_y_cap
        self.cone_force = self.dynamic_pressure * self.reference_area * np.array(self.c_y_cone)

        '''
        Force calculations subsonic
        '''
        self.subsonic_force = self.dynamic_pressure * self.reference_area * 0.35

        '''
        Calculating cm_alpha hypersonic
        '''
        self.cap_centroid = self.sphere_radius
        x_b = self.sphere_radius * (1 - np.cos(self.mu_b))
        self.cone_centroid = x_b + (self.cone_length - ((2 * self.cone_length)/(3 * (np.cos(self.cone_half_angle))**2)))
        self.total_moment = (self.cap_centroid * self.cap_force) + (self.cone_centroid * self.cone_force)
        self.total_moment_arm = ((self.sphere_radius * self.cap_force) + (self.cone_centroid * self.cone_force)) / (self.cap_force + self.cone_force)
        self.moment_xcg = (self.total_moment - xcg*self.cap_force) * -1
        self.c_moment = self.moment_xcg / (self.dynamic_pressure * self.reference_area * self.cone_max_radius)

        '''
        cm_alpha subsonic
        '''
        self.total_moment_subsonic = self.subsonic_force * self.total_moment_arm
        self.moment_xcg_subsonic = (self.total_moment_subsonic - xcg*self.subsonic_force) * -1
        self.c_moment_subsonic = self.moment_xcg_subsonic / (self.dynamic_pressure * self.reference_area * self.cone_max_radius)
        
        min_index = np.where(self.c_moment == min(self.c_moment))
        max_index = np.where(self.c_moment == max(self.c_moment))

        alpha_min_cm = self.aoa_all[min_index]
        alpha_max_cm = self.aoa_all[max_index]

        self.cm_alpha = (self.c_moment[min_index] / (alpha_min_cm))
        '''
        Oscillations constant velocity, constant density
        '''
        print("\033[91mCheck if the pitch-damping and axial-coeffs value have been changed for what you want\033[0m")

        self.pitch_damping = [-0.3775, -0.365, -0.35]
        #self.axial_coeff = [1,2,0.5]
        c_a = 1.75

        t1 = 2*self.mass / (density * velocity * self.reference_area * c_a)
        self.time = np.arange(t1, t1+0.5, 0.001)
        self.alpha_const_all = []
        self.alpha_decel_all = []
        for i in range(3):
            cmq_cmadot = self.pitch_damping[i]
            c_a = 2
            eta_1 = (density * velocity * self.reference_area * (self.cone_max_radius * 2)**2 ) * cmq_cmadot / (8 * self.moment_inertia)
            omega = np.sqrt(-1 * density * velocity**2 * self.reference_area * (self.cone_max_radius * 2) * self.cm_alpha / (2 * self.moment_inertia))

            alpha_const = (deviation * (np.pi / 180) * np.exp(eta_1 * self.time) * np.cos(omega * self.time))
            self.alpha_const_all.append(alpha_const)

        '''
        Oscillations decelleration, constant density
        '''
        for i in range(3):
            cmq_cmadot = self.pitch_damping[i]
            c_a = 2
            mu = ((self.mass * (self.cone_max_radius*2)**2 * cmq_cmadot) / (4 * self.moment_inertia * c_a)) + 1
            nu_sqrt = mu**2 - ((2*self.mass*self.cone_max_radius*2*self.cm_alpha)/(density*self.reference_area*self.moment_inertia*c_a**2))
            nu1 = (8*self.mass**2 * self.cm_alpha) / (np.pi * density * self.cone_max_radius*2 * self.moment_inertia * c_a**2) 
            nu_sqrt = (mu**2 - nu1)
            nu = np.sqrt(nu_sqrt)
            delta = (np.arctan2(mu,nu) - nu*np.log(t1)) 
            #delta = 0
            A = deviation / (t1**mu * np.cos(nu * np.log(t1) + delta))
            alpha_decel = A * ((self.time)**mu) * np.cos(nu * np.log(self.time) + delta)

            self.alpha_decel_all.append(alpha_decel)

        

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
            body.hypersonic_aerodynamics()
            max_ld_index = np.argmax(body.lift_over_drag)
            max_ld.append(body.lift_over_drag[max_ld_index])
            max_ld_aoa.append(body.aoa_all[max_ld_index])

        x1 = taper_ratios
        y = np.array(max_ld_aoa) * 180/np.pi
        z1 = np.array(max_ld)

        plt.show()

        '''
        Analysing effect of sphericity of base
        '''
        max_ld = []
        max_ld_aoa = []
        sphericity_radius = np.linspace(0.1, cone_max_radius, 100)
        for s in sphericity_radius:
            body = BluntBody(cone_length = 15, cone_max_radius = 5, cone_min_radius = 1, base_arc_height = s, mass = 28000)
            body.hypersonic_aerodynamics()
            max_ld_index = np.argmax(body.lift_over_drag)
            max_ld.append(body.lift_over_drag[max_ld_index])
            max_ld_aoa.append(body.aoa_all[max_ld_index])

        x2 = sphericity_radius
        y = np.array(max_ld_aoa) * 180/np.pi
        z2 = np.array(max_ld)
        '''
        Plots
        '''
        colors = {
        "cyan": "#00FFFF",
        "black": "#000000",
        "white": "#FFFFFF",
        "darkblue": "#0C2340",
        "turquoise": "#00B8C8",
        "royalblue": "#0076C2",
        "purple": "#6F1D77",
        "pink": "#EF60A3",
        "bordeaux": "#A50034",
        "red": "#E03C31",
        "orange": "#EC6842",
        "yellow": "#FFB81C",
        "green": "#6CC24A",
        "forestgreen": "#009B77",
        "darkgray": "#5C5C5C"
        }
        plt.figure()
        plt.scatter(x2, z2, label="Effect of Base Bluntness", color="purple")
        plt.xlabel("Arc height of base sphere[m]", fontsize=12)
        plt.ylabel("Maximum L/D", fontsize=12)
        plt.legend()
        plt.grid(True)
        plt.savefig("aero_base_bluntness_effect.pdf")

        plt.figure()
        plt.scatter(x1, z1, label="Effect of Cone Taper", color="royalblue")
        plt.xlabel("Taper ratio of Cone[min radius/max radius]", fontsize=12)
        plt.ylabel("Maximum L/D", fontsize=12)
        plt.legend()
        plt.grid(True)
        plt.savefig("aero_cone_taper_effect.pdf")

        plt.tight_layout()
        plt.show()

    def plots(self):
        # '''
        # Aerodynamic coefficient plots
        # '''
        # plt.plot(self.aoa_all * 180/np.pi, self.c_y_cap + self.c_y_cone, label = "normal force coeff", color = "blue")
        # plt.plot(self.aoa_all * 180/np.pi, self.c_x_cap + self.c_x_cone, label = "axial force coeff", color = "red")
        # plt.plot(self.aoa_all * 180/np.pi, self.normal_over_axial, label = "normal / axial", color = "green")
        # plt.figure()
        # plt.plot(self.aoa_all * 180/np.pi, self.lift_over_drag, label="lift over drag", color = 'red')
        # plt.plot(self.aoa_all * 180/np.pi, self.c_l, label="lift coefficient", color = 'green')
        # plt.plot(self.aoa_all * 180/np.pi, self.c_d, label="drag coefficient", color = 'blue')
        # plt.legend()

        # plt.figure()
        # plt.plot(self.aoa_all * 180/np.pi, self.c_axial, label="axial coefficient", color = 'green')
        # plt.plot(self.aoa_all * 180/np.pi, self.c_normal, label="normal coefficient", color = 'blue')
        # plt.legend()

        # # plt.figure()
        # # plt.plot(self.RASAero_mach, self.cd_4, label='C_d vs Mach at a = 4 / RASAero')
        # # plt.legend()

        # # plt.figure()
        # # plt.plot(self.RASAero_mach, self.cl_4, label='C_l vs Mach at a = 4 / RASAero')
        # # plt.legend()

        # '''
        # Stability plots
        # '''
        # plt.figure()
        # plt.plot(self.aoa_all * 180/np.pi, self.cmq_cmadot, label='pitch dampign coefficnet')
        # plt.legend()

        # plt.figure()
        # plt.plot(self.aoa_all * 180/np.pi, self.c_moment, label='moment coefficient')
        # plt.legend()

        # plt.figure()
        # plt.plot(self.aoa_all * 180/np.pi, self.c_moment, label='moment coefficient')
        # plt.legend()

        # plt.figure()
        # plt.plot(self.time, self.alpha_decel, label='stability for decelerating', color = "red")
        # plt.plot(self.time, self.alpha_const * 180/np.pi, label="stability for constant velocity", color="blue")
        # plt.legend()

        # plt.show()

        # Define TU Delft house colors
        colors = {
            "cyan": "#00FFFF",
            "black": "#000000",
            "white": "#FFFFFF",
            "darkblue": "#0C2340",
            "turquoise": "#00B8C8",
            "royalblue": "#0076C2",
            "purple": "#6F1D77",
            "pink": "#EF60A3",
            "bordeaux": "#A50034",
            "red": "#E03C31",
            "orange": "#EC6842",
            "yellow": "#FFB81C",
            "green": "#6CC24A",
            "forestgreen": "#009B77",
            "darkgray": "#5C5C5C"
        }

        deg = self.aoa_all * 180 / np.pi
        # --- Lift, Drag and Ratio Subsonic ---
        plt.figure(figsize=(10, 6))
        plt.plot(self.mach_subsonic_reentry, self.cl_subsonic_reentry / self.cd_subsonic_reentry, label="Subsonic L/D (Re-entry)", color=colors['bordeaux'], linewidth=2)
        plt.plot(self.mach_subsonic_reentry, self.cl_subsonic_reentry, label="Subsonic Lift Coefficient (Re-entry)", color=colors["green"], linewidth=2)
        plt.plot(self.mach_subsonic_reentry, self.cd_subsonic_reentry, label="Subsonic Drag Coefficient (Re-entry)", color=colors["royalblue"], linewidth=2)
        #plt.title("Lift and Drag Characteristics", fontsize=14)
        plt.xlabel("Mach Number", fontsize=12)
        plt.ylabel("Coefficient", fontsize=12)
        plt.legend()
        plt.grid(True)
        plt.savefig("aero_lift-drag-coeffs-subsonic-reentry.pdf")

        plt.figure(figsize=(10, 6))
        plt.plot(self.mach_subsonic_reentry, self.cl_subsonic_liftoff / self.cd_subsonic_liftoff, label="Subsonic L/D (Launch)", color=colors['bordeaux'], linewidth=2)
        plt.plot(self.mach_subsonic_liftoff, self.cl_subsonic_liftoff, label="Subsonic Lift Coefficient (Launch)", color=colors["green"], linewidth=2)
        plt.plot(self.mach_subsonic_liftoff, self.cd_subsonic_liftoff, label="Subsonic Drag Coefficient (Launch)", color=colors["royalblue"], linewidth=2)
        #plt.title("Lift and Drag Characteristics", fontsize=14)
        plt.xlabel("Mach Number", fontsize=12)
        plt.ylabel("Coefficient", fontsize=12)
        plt.legend()
        plt.grid(True)
        plt.savefig("aero_lift-drag-coeffs-subsonic-launch.pdf")

        # --- Lift, Drag and Ratio Hypersonic ---
        plt.figure(figsize=(10, 6))
        plt.plot(deg, self.lift_over_drag, label="Hypersonic L/D", color=colors["bordeaux"], linewidth=2)
        plt.plot(deg, self.c_l, label="Hypersonic Lift Coefficient", color=colors["green"], linewidth=2)
        plt.plot(deg, self.c_d, label="Hypesonic Drag Coefficient", color=colors["royalblue"], linewidth=2)
        #plt.title("Lift and Drag Characteristics", fontsize=14)
        plt.xlabel("Angle of Attack (deg)", fontsize=12)
        plt.ylabel("Coefficient", fontsize=12)
        plt.legend()
        plt.grid(True)
        plt.savefig("aero_lift-drag-coeffs.pdf")

        # --- Axial and Normal Hypersonic---
        plt.figure(figsize=(10, 6))
        plt.plot(deg, self.c_axial, label="Axial Coefficient", color=colors["royalblue"], linewidth=2)
        plt.plot(deg, self.c_normal, label="Normal Coefficient", color=colors["forestgreen"], linewidth=2)
        #plt.title("Axial vs Normal Coefficients", fontsize=14)
        plt.xlabel("Angle of Attack (deg)", fontsize=12)
        plt.ylabel("Coefficient", fontsize=12)
        plt.legend()
        plt.grid(True)
        plt.savefig("aero_normal-axial-coeffs.pdf")

        # --- Pitch Damping Coefficient ---
        plt.figure(figsize=(10, 6))
        plt.plot(deg, self.cmq_cmadot, label='Pitch Damping Coefficient', color=colors["purple"], linewidth=2)
        #plt.title("Pitch Damping Coefficient", fontsize=14)
        plt.xlabel("Angle of Attack (deg)", fontsize=12)
        plt.ylabel("Coefficient", fontsize=12)
        plt.legend()
        plt.grid(True)
        plt.savefig("aero_pitch-damping-coeff.pdf")

        # --- Moment Coefficient ---
        plt.figure(figsize=(10, 6))
        plt.plot(deg, self.c_moment, label='Moment Coefficient', color=colors["bordeaux"], linewidth=2)
        #plt.title("Moment Coefficient vs AoA", fontsize=14)
        plt.xlabel("Angle of Attack (deg)", fontsize=12)
        plt.ylabel("Moment Coefficient", fontsize=12)
        plt.legend()
        plt.grid(True)
        plt.savefig("aero_moment-coeffs.pdf")

        # --- Stability over Time ---
        plt.figure(figsize=(10, 6))
        plt.plot(self.time, self.alpha_const_all[0] * 180 / np.pi, label=f"cmq_cmadot = {self.pitch_damping[0]}", color=colors["darkblue"], linewidth=2,)
        plt.plot(self.time, self.alpha_const_all[1] * 180 / np.pi, label=f"cmq_cmadot = {self.pitch_damping[1]}", color=colors["bordeaux"], linewidth=2)
        plt.plot(self.time, self.alpha_const_all[2] * 180 / np.pi, label=f"cmq_cmadot = {self.pitch_damping[2]}", color=colors["forestgreen"], linewidth=2)
        #plt.title("Stability vs Time", fontsize=14)
        plt.xlabel("Time (s)", fontsize=12)
        plt.ylabel("Angle of Attack (deg)", fontsize=12)
        plt.legend()
        plt.grid(True)
        plt.savefig("aero_stability_constvelocity.pdf")

        plt.figure(figsize=(10, 6))
        plt.plot(self.time, self.alpha_decel_all[0], label=f"cmq_cmadot = {self.pitch_damping[0]}", color=colors["darkblue"], linewidth=2, linestyle = "dotted", alpha = 0.5)
        plt.plot(self.time, self.alpha_decel_all[1], label=f"cmq_cmadot = {self.pitch_damping[1]}", color=colors["bordeaux"], linewidth=2)
        plt.plot(self.time, self.alpha_decel_all[2], label=f"cmq_cmadot = {self.pitch_damping[2]}", color=colors["forestgreen"], linewidth=2, linestyle = "dotted", alpha = 0.5)
        plt.xlabel("Time (s)", fontsize=12)
        plt.ylabel("Angle of Attack (deg)", fontsize=12)
        plt.legend()
        plt.grid(True)
        plt.savefig("aero_stability_decelvelocity.pdf")

        plt.tight_layout()
        plt.show()

def run():
    # APOLLO CAPSULE
    # cone_length = 2.662
    # cone_max_radius = 1.956
    # cone_min_radius = 0.01
    # base_arc_height = 0.5
    # mass = 28000

    cone_length = 30
    cone_max_radius = 10.625/2
    cone_min_radius = 2.5
    base_arc_height = 0.5
    mass = 40000

    # cone_length = 13.95
    # cone_max_radius = 5
    # cone_min_radius = 2.5
    # base_arc_height = 2
    # mass = 28000

    body = BluntBody(cone_length, cone_max_radius, cone_min_radius, base_arc_height, mass)
    print("sphere radius: ", body.sphere_radius)
    #body.draw_geometry()
    body.hypersonic_aerodynamics()
    body.stability(velocity = 4000, density = 4.0084E-2, deviation = 3, moment_inertia=231368.0069, aoa = 16)
    #body.analysis()
    body.subsonic_aerodynamics(file_name_base= "HermesV1-RASAero.csv", file_name_nose="HermesV1-RASAero-Nose.csv")
    body.plots()


if __name__ == "__main__":
    run()
    






                
                

