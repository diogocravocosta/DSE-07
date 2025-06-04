# external
import numpy as np
import math

# internal
import src.data.constants as ct

class Atmosphere:
    """
    Class for atmospheric calculations.

    Parameters:
        altitude (float): Altitude at which to perform calculations
    """
    def __init__(self,
                 altitude,
                 model) -> None:


        # VARIABLES AND CONSTANTS
        # gravitational acceleration at surface
        self.g_0 = ct.g_0           # m/s^2
        # radius of earth
        self.Re = ct.radius         # m
        # height
        self.h = altitude * 1000      # m

        # gas constants
        self.R_star = ct.R_star
        self.M = ct.molecular_mass
        self.R = self.R_star / self.M
        self.gamma = ct.gamma
        self.rho_0 = ct.density_sea_level

        # CALCULATIONS
        self.g = self.gravitational_acceleration()
        self.geop_altitude = self.geopotential_altitude()

        # EXPONENTIAL MODEL

        self.T_exp, self.scale_height, self.beta, self.speed_of_sound = self.exponential_atmosphere(scale_height=7050)

        self.rho_exp = self.exponential_density()



    def gravitational_acceleration(self) -> float:
        """
            Calculates the gravitational acceleration at a certain altitude.

            Returns:
                g: gravitational acceleration

        """
        g = self.g_0 / (1 + (self.h / self.Re))**2

        return g

    def geopotential_altitude(self):
        """
            Calculates the geopotential altitude equivalent at a certain altitude.

            Returns:
                gp_altitude : geopotential altitude

        """
        gp_altitude = self.h * (1 - self.h / self.Re)

        return gp_altitude

    # FOR EXPONENTIAL ATMOSPHERE MODEL
    def exponential_atmosphere(self, scale_height: float) -> tuple[int, float, float, float]:
        """
            Calculates the gravitational acceleration at a certain altitude.

            Returns:
                T_exp: the assumed constant temperature for the exponential atmosphere model
                scale_height: the scale height of the model
                beta: 1/scale_height
                speed_of_sound: the speed of sound at the given constant temperature

        """

        if scale_height is None:
            h_s = 7050
        else:
            h_s = scale_height

        T_exp = math.floor(h_s * self.g_0 / self.R)

        beta = 1 / scale_height

        speed_of_sound = round(math.sqrt(self.gamma * self.R * T_exp),1)

        return T_exp, scale_height, beta, speed_of_sound

    def exponential_density(self) -> float:
        """
            Calculates the exponential density at a certain altitude.
            Returns:
                rho: the density calculated using the exponential density

        """
        rho = self.rho_0 * np.exp(-self.beta * self.h)
        return rho


