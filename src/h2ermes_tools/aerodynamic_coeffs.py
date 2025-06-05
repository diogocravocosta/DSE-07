#imports 
import numpy as np
import matplotlib.pyplot as plt

class bluntBody:
    '''
    Class: Defining geometry and calculating aerodynamic characteristics of blunt bodies
    Author: Soham Katewale
    '''
    def __init__(self, 
                 a: float, #Cone length
                 B: float, #Cone maximum radius
                 c: float, #Arc height of spherical base
                 m: float, #Mass of vehicle
                 )
        #initialisations
        self.a = a
        self.B = B
        self.c = c
        self.m = m

        #Calculated geometry 
        self.R = (2*B)**2 / (8*C) + c/2 #Radius of the hemisphere (base spherical base of the blunt body)
        assert self.R >= B, "Hemisphere radiys is too small, R must be >= B"
        self.mu_b = np.arcos((Self.R - c)/self.R) #half arc angle of hemisphere
        self.alpha_max = np.pi/2 - self.mu_b #maximum angle of attack before requiring numerical integration
        self.w = np.arctan(B / ((B*a) / (B-Rc))) #cone half angle

        def compute_aerodynamics(self, mode):
            '''
            Function: Computes the aerodynamic coefficients
            Mode 1: Computes aerodynamic coefficients for a range of angle of attacks
            Mode 2: Computes axial and normal force and thus the drag and lift forces over a range of velocities and densities
            '''
            if mode == 1:
                self.alpha = np.linspace(0, self.alpha_max , 100)
                '''
                C_x_cap = axial force coefficient acting on the cap (spherical base of the blunt body)
                C_y_cap = normal force coefficient acting on the cap
                C_x_cone = axial force coefficient acting on the cone
                C_y_cone = normal force coefficient acting on the cone
                '''
                self.C_x_cap = 0.5 * (np.sin(self.alpha)**2) * (np.sin(self.mu_b)**2) + \
                (1 + np.cos(self.mu_b)**2) * (np.cos(self.alpha)**2)
                self.C_y_cap = np.sin(self.alpha) * np.cos(self.alpha) * (np.sin(self.mu_b)**2)

                self.C_x_cone = []
                self.C_y_cone = []
                for a in self.alpha:
                    rho_str = np.arccos(-1 * np.tan(self.w) / np.tan(a))

