# external
import numpy as np

# internal
import src.data.constants as ct

class Atmosphere:
    def __init__(self,
                 height,
                 type = "exponential",
                 constant_temp = 240):
        # VARIABLES AND CONSTANTS
        # gravitational acceleration at surface
        self.g_0 = ct.g_0           # m/s^2
        # radius of earth
        self.Re = ct.radius         # m
        # height
        self.h = height             # m
        # type : whether exponential atmosphere used for simplified calculations
        self.type = type
        # gas constants
        self.R_star = ct.R_star
        self.M = ct.molecular_mass
        self.R = self.R_star / self.M


        # CALCULATIONS
        # gravitational acceleration
        self.g = self.gravitational_acceleration()      # m/s^2

        if self.type == "exponential":
            # assume isothermal and g = g0
            self.scale_height = 7200
            self.beta = 1/self.scale_height
            self.rho_0 = ct.density_sea_level

            self.T_exp = 246
            self.rho_exp = self.exponential_density()





    def gravitational_acceleration(self):
        g = self.g_0 / (1 + self.h / self.Re)**2

        return g

    def geopotential_altitude(self):
        z = self.h * (1 - self.h / self.Re)

    # FOR EXPONENTIAL MODEL
    def exponential_density(self):
        rho = self.rho_0 * np.exp(-self.beta * self.h)
        return rho


# class Exponential:
#     """Class for exponential atmosphere to be used for quick calculations."""
#     def __init__(self,
#                  height):

#
# class Standard(Atmosphere):
#     """Class for standard atmosphere to be used for more detailed calculations."""
#     def __init__(self,
#                  height):
#         super().__init__()