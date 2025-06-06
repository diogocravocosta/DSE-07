# external
import numpy as np
import math
import pandas as pd
import hvplot.pandas
from pyatmos import download_sw_nrlmsise00,read_sw_nrlmsise00
from pyatmos import nrlmsise00

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
        self.T_std = self.atmosphere()
        self.rho_std = 0  # TODO: add values

    # FOR STANDARD ATMOSPHERE MODEL
    def gravitational_acceleration_std(self) -> float:
        """
            Calculates the gravitational acceleration used in standard atmosphere calculations.

            Returns:
                g: the gravitational acceleration

        """
        g = self.g_0 * (self.R_0 / (self.R_0 + self.altitude))

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
            # TODO: add density
            scale = -6.5
            T_final = self.T_0 + scale * ((self.altitude_gp-0) / 1000)
            K = (self.g_0*self.M*1000/(self.R_star*scale))
            p_final = self.p_0 * (self.T_0/(self.T_0 + scale * ((self.altitude_gp-0)/1000)))**K
        else:
            # TODO: add density
            scale = -6.5
            T_1 = self.T_0 + scale * (11000 / 1000)
            K = (self.g_0 * self.M * 1000 / (self.R_star * scale))
            p_1 = self.p_0 * (self.T_0/(self.T_0 + scale * ((11000-0) / 1000)))**K

            # layer 2 (11019 - 20062)
            if 11019 <= self.altitude < 20062:
                # TODO: add density
                scale = 0.0
                T_final = T_1 + scale * ((self.altitude_gp-11000) / 1000)
                p_final = p_1 * np.exp(-self.g_0*self.M* (self.altitude_gp-11000)/(self.R_star*T_1))
            elif 20062 <= self.altitude:
                # TODO: add density
                scale = 0.0
                T_2 = T_1 + scale * ((20000- 11000) / 1000)
                p_2 = p_1 * np.exp(-self.g_0*self.M* (20000-11000)/(self.R_star*T_1))

                # layer 3 (20062 - 32161)
                if 20062 <= self.altitude < 32161:
                    # TODO: add density
                    scale = 1.0
                    T_final = T_2 + scale * ((self.altitude_gp-20000) / 1000)
                    K = (self.g_0 * self.M * 1000 / (self.R_star * scale))
                    p_final = p_2 * (T_2 / (T_2 + scale * ((self.altitude_gp - 20000) / 1000))) ** K

                elif 32161 <= self.altitude:
                    # TODO: add density
                    scale = 1.0
                    T_3 = T_2 + scale * ((32000-20000) / 1000)
                    K = (self.g_0 * self.M * 1000 / (self.R_star * scale))
                    p_3 = p_2 * (T_2 / (T_2 + scale * ((32000 - 20000) / 1000))) ** K

                    # layer 4 (32161 - 47348)
                    if 32161 <= self.altitude < 47348:
                        # TODO: add density
                        scale = 2.8
                        T_final = T_3 + scale * ((self.altitude_gp-32000) / 1000)
                        K = (self.g_0 * self.M * 1000 / (self.R_star * scale))
                        p_final = p_3 * (T_3 / (T_3 + scale * ((self.altitude_gp - 32000) / 1000))) ** K
                    elif 47348 <= self.altitude:
                        # TODO: add density
                        scale = 2.8
                        T_4 = T_3 + scale * ((47000 - 32000) / 1000)
                        K = (self.g_0 * self.M * 1000 / (self.R_star * scale))
                        p_4 = p_3 * (T_3 / (T_3 + scale * ((47000 - 32000) / 1000))) ** K

                        # layer 5 (47348 - 51411)
                        if 47348 <= self.altitude < 51411:
                            # TODO: add density
                            scale = 0.0
                            T_final = T_4 + scale * ((self.altitude_gp-47000) / 1000)
                            p_final = p_4 * np.exp(-self.g_0*self.M* (self.altitude_gp-47000)/(self.R_star*T_4))

                        elif 51411 <= self.altitude:
                            # TODO: add density
                            scale = 0.0
                            T_5 = T_4 + scale * ((51000 - 47000) / 1000)
                            p_5 = p_4 * np.exp(-self.g_0*self.M* (51000-47000)/(self.R_star*T_4))

                            # layer 6 (51411 - 71799)
                            if 51411 <= self.altitude < 71799:
                                # TODO: add density
                                scale = -2.8
                                T_final = T_5 + scale * ((self.altitude_gp-51000) / 1000)
                                K = (self.g_0 * self.M * 1000 / (self.R_star * scale))
                                p_final = p_5 * (T_5 / (T_5 + scale * ((self.altitude_gp - 51000) / 1000))) ** K

                            elif 71799 <= self.altitude:
                                # TODO: add density
                                scale = -2.8
                                T_6 = T_5 + scale * ((71000 - 51000) / 1000)
                                K = (self.g_0 * self.M * 1000 / (self.R_star * scale))
                                p_6 = p_5 * (T_5 / (T_5 + scale * ((71000 - 51000) / 1000))) ** K

                                # layer 7 (71799 - 86000)
                                if 71799 <= self.altitude < 86000:
                                    # TODO: add density
                                    scale = -2.0
                                    T_final = T_6 + scale * ((self.altitude_gp-71000) / 1000)
                                    K = (self.g_0 * self.M * 1000 / (self.R_star * scale))
                                    p_final = p_6 * (T_6 / (T_6 + scale * ((self.altitude_gp - 71000) / 1000))) ** K
                                elif 86000 <= self.altitude:
                                    # TODO: add density
                                    scale = -2.0
                                    T_7 = T_6 + scale * ((86000 - 71000) / 1000)
                                    K = (self.g_0 * self.M * 1000 / (self.R_star * scale))
                                    p_7 = p_6 * (T_6 / (T_6 + scale * ((86000 - 71000) / 1000))) ** K

                                    # layer 8 (86000 - 91000 89716)
                                    if 86000 <= self.altitude < 91000:
                                        # TODO: add pressure
                                        # TODO: add density
                                        scale = 0.0
                                        T_final = T_7 + scale * ((self.altitude-86000) / 1000)
                                        p_final = p_7 * np.exp(-self.g_0*self.M* ((self.altitude_gp-86000) / 1000)/(self.R_star*T_7))

                                    elif 91000 <= self.altitude:
                                        # TODO: add pressure
                                        # TODO: add density
                                        scale = 0.0
                                        T_8 = T_7 + scale * ((91000 - 86000) / 1000)
                                        p_8 = p_7 * np.exp(-self.g_0*self.M* ((91000 - 86000) / 1000)/(self.R_star*T_7))

                                        # layer 9 (91000 - 110000)
                                        if 91000 <= self.altitude < 110000:
                                            # TODO: add pressure
                                            # TODO: add density
                                            T_c = 263.1905
                                            A = -76.3232
                                            a = -19.9429
                                            T_final = T_c + A * (1-((((self.altitude-91000)/1000)/a)**2)) ** 0.5
                                        elif 110000 <= self.altitude:
                                            # TODO: add pressure
                                            # TODO: add density
                                            T_c = 263.1905
                                            A = -76.3232
                                            a = -19.9429
                                            T_9 = T_c + A * (1 - ((((110000 - 91000)/1000) / a) ** 2)) ** 0.5

                                            # layer 10 (110000 - 120000)
                                            if 110000 <= self.altitude < 120000:
                                                # TODO: add pressure
                                                # TODO: add density
                                                scale = 12.0
                                                T_final = T_9 + scale * ((self.altitude-110000) / 1000)

                                            elif 120000 <= self.altitude:
                                                # TODO: add pressure
                                                # TODO: add density
                                                scale = 12.0
                                                T_10 = T_9 + scale * ((120000 - 110000) / 1000)

                                                # layer 11 (120000 - 1000000)
                                                if 120000 <= self.altitude < 1000000:
                                                    # TODO: add pressure
                                                    # TODO: add density
                                                    lbd = scale / (self.T_inf - T_10)
                                                    xi = ((self.altitude - 120000) * (self.R_0 + 120000) / (self.R_0 + self.altitude))/1000
                                                    T_final = self.T_inf - (self.T_inf - T_10) * np.exp(- lbd * xi)
        print(T_final, p_final)
        return T_final, p_final

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


# TEST
# temps = []
# hs = []
# i = 0
# while i <= 1000:
#     atmosphere = StandardAtmosphere(i)
#     temp = atmosphere.T_std
#     hs.append(i)
#     temps.append(temp)
#     i += 0.01
#
# df = pd.DataFrame({'temp': temps, 'h': hs})
#
# plot = df.hvplot(x = 'temp', y = 'h').opts(width = 400, height = 800)
# hvplot.show(plot)



