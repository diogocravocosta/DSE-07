# external
import numpy as np
import math
import pandas as pd
from pyatmos import download_sw_nrlmsise00,read_sw_nrlmsise00
from pyatmos import nrlmsise00
from scipy.interpolate import CubicSpline

# internal
import data.constants as ct

# MAYBE MOVE TO RE-ENTRY FOLDER?

class ExponentialAtmosphere:
    """
    Class for exponential atmospheric calculations.

    Parameters:
        altitude (float): Altitude at which to perform calculations
    """
    def __init__(self,
                 altitude,
                 ) -> None:
        # VARIABLES AND CONSTANTS
        # gravitational acceleration at surface
        self.g_0 = ct.g_0           # m/s^2
        # radius of earth
        self.Re = ct.earth_radius         # m
        # height
        self.altitude = altitude * 1000      # m

        # gas constants
        self.R_star = ct.R_star
        self.M = ct.molecular_mass
        self.R = self.R_star / self.M
        self.gamma = ct.gamma
        self.rho_0 = ct.density_sea_level

        # CALCULATIONS

        self.altitude_gp = self.geopotential_altitude()

        # EXPONENTIAL MODEL
        self.g = self.gravitational_acceleration()

        self.T_exp, self.scale_height, self.beta, self.speed_of_sound = self.atmosphere(scale_height=7050)

        self.rho_exp = self.density()

    def gravitational_acceleration(self) -> float:
        """
            Calculates the gravitational acceleration at a certain altitude.

            Returns:
                g: gravitational acceleration

        """
        g = self.g_0 / (1 + (self.altitude / self.Re))**2

        return g

    def geopotential_altitude(self):
        """
            Calculates the geopotential altitude equivalent at a certain altitude.

            Returns:
                gp_altitude : geopotential altitude

        """
        gp_altitude = self.altitude * (1 - self.altitude / self.Re)

        return gp_altitude

    # FOR EXPONENTIAL ATMOSPHERE MODEL
    def atmosphere(self, scale_height: float) -> tuple[int, float, float, float]:
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

    def density(self) -> float:
        """
            Calculates the exponential density at a certain altitude.
            Returns:
                rho: the density calculated using the exponential density

        """
        rho = self.rho_0 * np.exp(-self.beta * self.altitude)
        return rho

class StandardAtmosphere:
    def __init__(self,
                 altitude: float,
                 ) -> None:
        # VARIABLES AND CONSTANTS
        # gravitational acceleration at surface
        self.g_0 = ct.g_0  # m/s^2
        # radius of earth
        self.Re = ct.earth_radius  # m
        # height
        self.altitude = altitude * 1000  # m

        # gas constants
        self.R_star = ct.R_star
        self.M = ct.molecular_mass
        self.R = self.R_star / self.M
        self.gamma = ct.gamma

        # sea-level
        self.rho_0 = ct.density_sea_level
        self.T_0 = ct.temperature_sea_level
        self.p_0 = ct.pressure_sea_level

        # infinity
        self.T_inf = 1000

        self.R_0 = 6356766  # m

        # CALCULATIONS
        self.g = self.gravitational_acceleration_std()
        self.altitude_gp = self.geopotential_altitude()
        self.T_std, self.p_std, self.rho_std = self.atmosphere()

    # FOR STANDARD ATMOSPHERE MODEL
    def gravitational_acceleration_std(self) -> float:
        """
            Calculates the gravitational acceleration used in standard atmosphere calculations.

            Returns:
                g: the gravitational acceleration

        """
        g = self.g_0 * (self.R_0 / (self.R_0 + self.altitude))**2

        return g

    def geopotential_altitude(self):
        """
            Calculates the geopotential altitude equivalent at a certain altitude.

            Returns:
                gp_altitude : geopotential altitude

        """
        gp_altitude = self.altitude * (1 - self.altitude / self.Re)

        return gp_altitude

    def atmosphere(self):
        atmosphere_df = pd.DataFrame({"geometric_altitude_min": [0, 11019, 20062, 32161, 47348, 51411, 71799, 86000],
                                      "geopotential_altitude": [0, 11000, 20000, 32000, 47000, 51000, 71000, 84852],
                                      "temperature_gradient": [-6.5, 0.0, 1.0, 2.8, 0.0, -2.8, -2.0, 0],
                                      })
        atmosphere_df.reset_index(inplace=True, names=["layer"])

        T_final = 0

        # layer 1 (0 - 11019)
        if self.altitude < 11019:
            scale = -6.5
            T_final = self.T_0 + scale * ((self.altitude_gp-0) / 1000)
            K = (self.g_0*self.M*1000/(self.R_star*scale))
            p_final = self.p_0 * (self.T_0/(self.T_0 + scale * ((self.altitude_gp-0)/1000)))**K
            rho_final = p_final * self.M / (self.R_star * T_final)
        else:
            scale = -6.5
            T_1 = self.T_0 + scale * (11000 / 1000)
            K = (self.g_0 * self.M * 1000 / (self.R_star * scale))
            p_1 = self.p_0 * (self.T_0/(self.T_0 + scale * ((11000-0) / 1000)))**K

            # layer 2 (11019 - 20062)
            if 11019 <= self.altitude < 20062:
                scale = 0.0
                T_final = T_1 + scale * ((self.altitude_gp-11000) / 1000)
                p_final = p_1 * np.exp(-self.g_0*self.M* (self.altitude_gp-11000)/(self.R_star*T_1))
                rho_final = p_final * self.M / (self.R_star * T_final)
            elif 20062 <= self.altitude:
                scale = 0.0
                T_2 = T_1 + scale * ((20000- 11000) / 1000)
                p_2 = p_1 * np.exp(-self.g_0*self.M* (20000-11000)/(self.R_star*T_1))

                # layer 3 (20062 - 32161)
                if 20062 <= self.altitude < 32161:
                    scale = 1.0
                    T_final = T_2 + scale * ((self.altitude_gp-20000) / 1000)
                    K = (self.g_0 * self.M * 1000 / (self.R_star * scale))
                    p_final = p_2 * (T_2 / (T_2 + scale * ((self.altitude_gp - 20000) / 1000))) ** K
                    rho_final = p_final * self.M / (self.R_star * T_final)

                elif 32161 <= self.altitude:
                    scale = 1.0
                    T_3 = T_2 + scale * ((32000-20000) / 1000)
                    K = (self.g_0 * self.M * 1000 / (self.R_star * scale))
                    p_3 = p_2 * (T_2 / (T_2 + scale * ((32000 - 20000) / 1000))) ** K

                    # layer 4 (32161 - 47348)
                    if 32161 <= self.altitude < 47348:
                        scale = 2.8
                        T_final = T_3 + scale * ((self.altitude_gp-32000) / 1000)
                        K = (self.g_0 * self.M * 1000 / (self.R_star * scale))
                        p_final = p_3 * (T_3 / (T_3 + scale * ((self.altitude_gp - 32000) / 1000))) ** K
                        rho_final = p_final * self.M / (self.R_star * T_final)
                    elif 47348 <= self.altitude:
                        scale = 2.8
                        T_4 = T_3 + scale * ((47000 - 32000) / 1000)
                        K = (self.g_0 * self.M * 1000 / (self.R_star * scale))
                        p_4 = p_3 * (T_3 / (T_3 + scale * ((47000 - 32000) / 1000))) ** K

                        # layer 5 (47348 - 51411)
                        if 47348 <= self.altitude < 51411:
                            scale = 0.0
                            T_final = T_4 + scale * ((self.altitude_gp-47000) / 1000)
                            p_final = p_4 * np.exp(-self.g_0*self.M* (self.altitude_gp-47000)/(self.R_star*T_4))
                            rho_final = p_final * self.M / (self.R_star * T_final)

                        elif 51411 <= self.altitude:
                            scale = 0.0
                            T_5 = T_4 + scale * ((51000 - 47000) / 1000)
                            p_5 = p_4 * np.exp(-self.g_0*self.M* (51000-47000)/(self.R_star*T_4))

                            # layer 6 (51411 - 71799)
                            if 51411 <= self.altitude < 71799:
                                scale = -2.8
                                T_final = T_5 + scale * ((self.altitude_gp-51000) / 1000)
                                K = (self.g_0 * self.M * 1000 / (self.R_star * scale))
                                p_final = p_5 * (T_5 / (T_5 + scale * ((self.altitude_gp - 51000) / 1000))) ** K
                                rho_final = p_final * self.M / (self.R_star * T_final)

                            elif 71799 <= self.altitude:
                                scale = -2.8
                                T_6 = T_5 + scale * ((71000 - 51000) / 1000)
                                K = (self.g_0 * self.M * 1000 / (self.R_star * scale))
                                p_final = p_5 * (T_5 / (T_5 + scale * ((71000 - 51000) / 1000))) ** K
                                p_6 = p_5 * (T_5 / (T_5 + scale * ((71000 - 51000) / 1000))) ** K

                                # layer 7 (71799 - 86000)
                                if 71799 <= self.altitude < 86000:
                                    scale = -2.0
                                    T_final = T_6 + scale * ((self.altitude_gp-71000) / 1000)
                                    K = (self.g_0 * self.M * 1000 / (self.R_star * scale))
                                    p_final = p_6 * (T_6 / (T_6 + scale * ((self.altitude_gp - 71000) / 1000))) ** K
                                    rho_final = p_final * self.M / (self.R_star * T_final)

                                elif 86000 <= self.altitude:
                                    scale = -2.0
                                    T_7_0 = T_6 + scale * ((84852 - 71000) / 1000)
                                    T_7 = T_7_0 * 0.999579
                                    K = (self.g_0 * self.M * 1000 / (self.R_star * scale))
                                    p_7 = p_6 * (T_6 / (T_6 + scale * ((84852 - 71000) / 1000))) ** K

                                    # layer 8 (86000 - 91000 89716)
                                    if 86000 <= self.altitude < 91000:
                                        scale = 0.0
                                        T_final = T_7 + scale * ((self.altitude-86000) / 1000)
                                        p_final = self.interpolate_pressure(self.altitude /1000)
                                        rho_final = self.interpolate_density(self.altitude / 1000)

                                    elif 91000 <= self.altitude:
                                        scale = 0.0
                                        T_8 = T_7 + scale * ((91000 - 86000) / 1000)

                                        # layer 9 (91000 - 110000)
                                        if 91000 <= self.altitude < 110000:
                                            T_c = 263.1905
                                            A = -76.3232
                                            a = -19.9429
                                            T_final = T_c + A * (1-((((self.altitude-91000)/1000)/a)**2)) ** 0.5
                                            p_final = self.interpolate_pressure(self.altitude /1000)
                                            rho_final = self.interpolate_density(self.altitude / 1000)

                                        elif 110000 <= self.altitude:
                                            T_c = 263.1905
                                            A = -76.3232
                                            a = -19.9429
                                            T_9 = T_c + A * (1 - ((((110000 - 91000)/1000) / a) ** 2)) ** 0.5

                                            # layer 10 (110000 - 120000)
                                            if 110000 <= self.altitude < 120000:
                                                scale = 12.0
                                                T_final = T_9 + scale * ((self.altitude-110000) / 1000)
                                                p_final = self.interpolate_pressure(self.altitude /1000)
                                                rho_final = self.interpolate_density(self.altitude / 1000)

                                            elif 120000 <= self.altitude:
                                                scale = 12.0
                                                T_10 = T_9 + scale * ((120000 - 110000) / 1000)

                                                # layer 11 (120000 - 1000000)
                                                if 120000 <= self.altitude < 1000000:
                                                    lbd = scale / (self.T_inf - T_10)
                                                    xi = ((self.altitude - 120000) * (self.R_0 + 120000) / (self.R_0 + self.altitude))/1000
                                                    T_final = self.T_inf - (self.T_inf - T_10) * np.exp(- lbd * xi)
                                                    p_final = self.interpolate_pressure(self.altitude /1000)
                                                    rho_final = self.interpolate_density(self.altitude /1000)

        return T_final, p_final, rho_final

    def interpolate_molecular_mass(self, height):
        # geometric altitude in km
        Z = [
            0, 5, 10, 15, 20, 25, 30, 35, 40, 45, 50, 55, 60, 65, 70, 75, 80, 85, 90, 95, 100,
            105, 110, 115, 120, 125, 130, 135, 140, 145, 150, 155, 160, 165, 170, 175, 180, 185, 190, 195,
            200, 205, 210, 215, 220, 225, 230, 235, 240, 245, 250, 255, 260, 265, 270, 275, 280, 285, 290, 295,
            300, 305, 310, 315, 320, 325, 330, 335, 340, 345, 350, 355, 360, 365, 370, 375, 380, 385, 390, 395,
            400, 405, 410, 415, 420, 425, 430, 435, 440, 445, 450, 455, 460, 465, 470, 475, 480, 485, 490, 495,
            500, 505, 510, 515, 520, 525, 530, 535, 540, 545, 550, 555, 560, 565, 570, 575, 580, 585, 590, 595,
            600, 605, 610, 615, 620, 625, 630, 635, 640, 645, 650, 655, 660, 665, 670, 675, 680, 685, 690, 695,
            700, 705, 710, 715, 720, 725, 730, 735, 740, 745, 750, 755, 760, 765, 770, 775, 780, 785, 790, 795,
            800, 805, 810, 815, 820, 825, 830, 835, 840, 845, 850, 855, 860, 865, 870, 875, 880, 885, 890, 895,
            900, 905, 910, 915, 920, 925, 930, 935, 940, 945, 950, 955, 960, 965, 970, 975, 980, 985, 990, 995, 1000
        ]

        # molecular mass
        M = [
            28.964, 28.964, 28.964, 28.964, 28.964, 28.964, 28.964, 28.964, 28.964, 28.964, 28.964, 28.964, 28.964,
            28.964, 28.964, 28.964, 28.964, 28.964, 28.964, 28.906, 28.734,
            28.4, 27.881, 27.267, 26.681, 26.205, 25.802, 25.437, 25.09, 24.753, 24.422, 24.1, 23.791, 23.49, 23.192,
            22.9, 22.616, 22.34, 22.072, 21.81, 21.552, 21.3, 21.056, 20.82, 20.591, 20.37, 20.156, 19.949, 19.749,
            19.556, 19.369, 19.19, 19.017, 18.85, 18.69, 18.536, 18.388, 18.246, 18.109, 17.978, 17.851, 17.73, 17.613,
            17.501, 17.393, 17.289, 17.188, 17.091, 16.998, 16.907, 16.82, 16.735, 16.652, 16.572, 16.493, 16.417,
            16.341, 16.267, 16.195, 16.122, 16.051, 15.98, 15.909, 15.838, 15.767, 15.695, 15.623, 15.55, 15.476,
            15.401, 15.324, 15.246, 15.166, 15.083, 14.999, 14.912, 14.823, 14.731, 14.636, 14.537, 14.435, 14.33,
            14.221, 14.108, 13.992, 13.871, 13.748, 13.621, 13.49, 13.356, 13.219, 13.079, 12.935, 12.788, 12.639,
            12.486, 12.33, 12.172, 12.01, 11.846, 11.679, 11.51, 11.338, 11.164, 10.988, 10.81, 10.631, 10.451, 10.27,
            10.089, 9.908, 9.727, 9.547, 9.367, 9.189, 9.012, 8.837, 8.664, 8.494, 8.326, 8.161, 8.0, 7.842, 7.688,
            7.538, 7.391, 7.248, 7.109, 6.973, 6.841, 6.712, 6.588, 6.466, 6.349, 6.235, 6.125, 6.018, 5.915, 5.816,
            5.72, 5.628, 5.54, 5.455, 5.374, 5.296, 5.222, 5.15, 5.082, 5.017, 4.955, 4.895, 4.839, 4.784, 4.733, 4.684,
            4.637, 4.593, 4.55, 4.51, 4.471, 4.435, 4.4, 4.367, 4.335, 4.305, 4.276, 4.249, 4.222, 4.197, 4.173, 4.15,
            4.128, 4.107, 4.086, 4.067, 4.047, 4.029, 4.01, 3.992, 3.975, 3.957, 3.94
        ]

        molecular_mass = np.interp(height, Z, M)
        return molecular_mass

    def interpolate_pressure(self, height):
        altitude = [86,100,115,130,150,175,200,250,300,400,500,600,700,800,900,1000]
        pressure_array = [3.7338e-1, 3.2011e-2, 4.0096e-3, 1.2505e-3, 4.5422e-4, 1.7936e-4, 8.4736e-5, 2.4767e-5, 8.7704e-6, 1.4518e-6, 3.0236e-7, 8.2130e-8, 3.1908e-8, 1.7036e-8, 1.0873e-8, 7.5138e-9]

        pressure_interp = CubicSpline(altitude, pressure_array)

        pressure_new = pressure_interp(height)

        return pressure_new

    def interpolate_density(self, height):
        altitude = [86, 100, 115, 130, 150, 175, 200, 250, 300, 400, 500, 600, 700, 800, 900, 1000]
        rho_array = [6.958e-6, 5.604e-7, 4.289e-8, 8.152e-9, 2.076e-9, 6.339e-10, 2.541e-10, 6.073e-11, 1.916e-11,
                     2.803e-12, 5.215e-13, 1.137e-13, 3.070e-14, 1.136e-14, 5.759e-15, 3.561e-15]

        rho_interp = CubicSpline(altitude, rho_array)
        rho_new = rho_interp(height)

        return rho_new







class NrlmsiseAtmosphere:
    def __init__(self,
                 altitude,
                 time,
                 latitude,
                 longitude,
                 ) -> None:
        # VARIABLES AND CONSTANTS
        # gravitational acceleration at surface
        self.g_0 = ct.g_0           # m/s^2
        # radius of earth
        self.Re = ct.earth_radius         # m
        # height
        self.altitude = altitude # km
        self.time = time # UTC
        self.latitude = latitude
        self.longitude = longitude

        # gas constants
        self.R_star = ct.R_star
        self.M = ct.molecular_mass
        self.R = self.R_star / self.M
        self.gamma = ct.gamma
        self.rho_0 = ct.density_sea_level

        # files
        self.sw_data = self.space_weather_data()

        # calculations
        self.T_nrl, self.rho_nrl = self.atmosphere()

    def space_weather_data(self):
        # download space weather file
        swfile = download_sw_nrlmsise00()

        # read space weather data
        swdata = read_sw_nrlmsise00(swfile)

        return swdata

    def atmosphere(self):
        nrl00 = nrlmsise00(self.time, (self.latitude, self.longitude, self.altitude), self.sw_data)

        T = nrl00.T
        rho = nrl00.rho

        return T, rho



