#imports 
import numpy as np
import matplotlib.pyplot as plt

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
        self.sphere_radius = (2*self.cone_max_radius)**2 / (8*self.cone_length) + self.cone_length/2 #Radius of the hemisphere (base spherical base of the blunt body)
        assert self.sphere_radius >= self.cone_max_radius, "Hemisphere radiys is too small, sphereRadius must be >= coneMaxRadius"
        self.mu_b = np.acos((self.sphere_radius - self.base_arc_height)/self.sphere_radius) #half arc angle of hemisphere
        self.alpha_max = np.pi/2 - self.mu_b #maximum angle of attack before requiring numerical integration
        self.cone_half_angle = np.arctan(self.cone_max_radius / ((self.cone_max_radius*self.cone_length) / (self.cone_max_radius-self.cone_min_radius))) #cone half angle

    def compute_aerodynamics(self, mode):
        '''
        Function: Computes the aerodynamic coefficients
        Mode 1: Computes aerodynamic coefficients for a range of angle of attacks
        Mode 2: Computes axial and normal force and thus the drag and lift forces over a range of velocities and densities

        Variables:
            c_x_cap = axial force coefficient acting on the cap (spherical base of the blunt body)
            c_y_cap = normal force coefficient acting on the cap
            c_x_cone = axial force coefficient acting on the cone
            c_y_cone = normal force coefficient acting on the cone
            aoa_all = whole range of angle of attacks
            rho_str = angle for tangential incidence at the cone
            c_l = lift coefficients
            c_d = drag coefficients
        '''
        if mode == 1:
            '''
            calculating axial and normal coefficients
            '''
            self.aoa_all = np.linspace(0, self.alpha_max , 100)
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
            
            '''
            calculating lift and drag coefficients
            '''
            self.c_l = (((self.c_y_cap + self.c_y_cone)*np.cos(self.aoa_all)) - ((self.c_x_cap + self.c_x_cone)*np.sin(self.aoa_all))) * -1
            self.c_d = ((self.c_y_cap + self.c_y_cone)*np.sin(self.aoa_all)) + ((self.c_x_cap + self.c_x_cone)*np.cos(self.aoa_all))
            self.lift_over_drag = self.c_l / self.c_d
            self.normal_over_axial = (self.c_y_cap + self.c_y_cone) / (self.c_x_cap + self.c_x_cone)

            '''
            making the plots
            '''
            #plt.plot(self.aoa_all * 180/np.pi, self.c_y_cap + self.c_y_cone, label = "normal force coeff", color = "blue")
            #plt.plot(self.aoa_all * 180/np.pi, self.c_x_cap + self.c_x_cone, label = "axial force coeff", color = "red")
            #plt.plot(self.aoa_all * 180/np.pi, self.normal_over_axial, label = "normal / axial", color = "green")
            plt.plot(self.aoa_all * 180/np.pi, self.lift_over_drag, label="lift over drag", color = 'red')
            plt.plot(self.aoa_all * 180/np.pi, self.c_l, label="lift coefficient", color = 'green')
            plt.plot(self.aoa_all * 180/np.pi, self.c_d, label="drag coefficient", color = 'blue')
            plt.legend()
            plt.show()

    def draw_geometry(self):
        '''
        Function: draws the geometry of the vehicle
        '''

        '''
        Hemisphere mesh grid
        '''
        mu = np.linspace(0, self.mu_b, 100)
        rho = np.linspace(0, 2*np.pi)
        mu, rho = np.meshgrid(mu, rho)

        '''
        Hemisphere parametric equations
        '''
        x_sphere = self.sphere_radius * (1 - np.cos(mu))
        y_sphere = self.sphere_radius * np.sin(mu) * np.cos(rho)
        z_sphere = self.sphere_radius * np.sin(mu) * np.sin(rho)

        '''
        Cone mesh grid
        '''
        x_cone = np.linspace(self.base_arc_height, self.base_arc_height + self.cone_length, 100)
        rho_cone = np.linspace(0, 2*np.pi, 100)
        x_cone, rho_cone = np.meshgrid(x_cone, rho_cone)

        '''
        Cone parametric equations
        '''
        r_cone = self.cone_max_radius - ((self.cone_max_radius - self.cone_min_radius)/self.cone_length) * (x_cone - self.base_arc_height)
        y_cone = r_cone * np.cos(rho_cone)
        z_cone = r_cone * np.sin(rho_cone)

        '''
        Plotting
        '''
        fig = plt.figure(figsize=(10, 7))
        ax = fig.add_subplot(111, projection='3d')

        # Hemisphere (front)
        ax.plot_surface(x_sphere, y_sphere, z_sphere, color='skyblue', alpha=0.9)

        # Cone (extending backward)
        ax.plot_surface(x_cone, y_cone, z_cone, color='lightcoral', alpha=0.9)

        # Labels and aspect
        ax.set_xlabel('x')
        ax.set_ylabel('y')
        ax.set_zlabel('z')
        ax.set_title('Blunt Body (Hemisphere + Cone)')
        ax.set_box_aspect([2, 1, 1])  # Wider x to show elongation

        plt.show()

if __name__ == "__main__":
    #CURRENT DIMENSIONS
    # cone_length = 13.95
    # cone_max_radius = 4.95
    # cone_min_radius = 2.46
    # base_arc_height = 2
    # mass = 28000

    cone_length = 13.95
    cone_max_radius = 5.5
    cone_min_radius = 2.46
    base_arc_height = 2
    mass = 28000

    # cone_length = 24.5
    # cone_max_radius = 4.95
    # cone_min_radius = 4.95/100
    # base_arc_height = 2
    # mass = 28000

    #APOLLO CAPSULE
    # cone_length = 2.662
    # cone_max_radius = 1.956
    # cone_min_radius = 0.01
    # base_arc_height = 0.5
    # mass = 28000

    body = BluntBody(cone_length, cone_max_radius, cone_min_radius, base_arc_height, mass
    )

    body.draw_geometry()

    body.compute_aerodynamics(mode = 1)

    base_surface_area = np.pi * 4.95**2
    flap_area = (np.pi * cone_max_radius**2) - base_surface_area
    max_ld_index = np.argmax(body.lift_over_drag)
    print("maximum l/d = ", body.lift_over_drag[max_ld_index])
    print("at angle of attack (deg) =", body.aoa_all[max_ld_index] * 180 / np.pi)
    print("flap area required = ", flap_area)

                
                

