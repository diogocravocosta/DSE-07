# external
import numpy as np

# internal
import src.data.constants as ct

class Atmosphere:
    def __init__(self,
                 height: float):
        # VARIABLES AND CONSTANTS
        # gravitational acceleration at surface
        self.g_0 = ct.g_0           # m/s^2
        # radius of earth
        self.Re = ct.radius         # m
        # height
        self.h = height * 1000      # m

        # gas constants
        self.R_star = ct.R_star
        self.M = ct.molecular_mass
        self.R = self.R_star / self.M
        self.gamma = ct.gamma
        self.rho_0 = ct.density_sea_level

        # CALCULATIONS
        self.g = self.gravitational_acceleration()

    def gravitational_acceleration(self):
        g = self.g_0 / (1 + self.h / self.Re)**2

        return g

    def geopotential_altitude(self):
        z = self.h * (1 - self.h / self.Re)

        return z

    # FOR EXPONENTIAL MODEL
    def exponential_density(self, height, beta):
        rho = self.rho_0 * np.exp(-beta * height)
        return rho

    def exponential_atmosphere(self, height, constant_temp = 240):
        self.rho_0 = ct.density_sea_level

        if constant_temp is None:
            T_exp = 246
            scale_height = 7200
        else:
            T_exp = constant_temp
            scale_height = (self.R * T_exp) / self.g_0

        beta = 1 / scale_height

        rho_exp = self.exponential_density(height, beta)

        speed_of_sound = np.sqrt(self.gamma * self.R * T_exp)

        return rho_exp, speed_of_sound

