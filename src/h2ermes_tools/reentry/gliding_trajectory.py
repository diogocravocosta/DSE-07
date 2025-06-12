# external
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from scipy.integrate import simpson

# internal
import data.constants as ct
from h2ermes_tools.reentry import atmosphere as atm

class GlidingEntry:
    def __init__(self,
                 planet,
                 altitude,
                 entry_speed,
                 vehicle_mass,
                 vehicle_surface,
                 nose_radius,
                 C_D,
                 C_L,
                 latitude_0,
                 longitude_0,
                 lift_parameter = None,
                 lift_drag_ratio = None,
                 atmosphere = "exponential",
                 boundary_layer = "laminar"
                 ):
        # ASSUMPTIONS: spherical non-rotating earth, atmosphere at rest, trajectory in flat plane, constant g
        # CONSTANTS
        # planetary
        self.planet = planet
        self.g_0, self.rho_0, self.R_e, self.scale_height, self.m = self.get_constants()


        self.beta = 1 / self.scale_height
        self.c_star = 1.1097e8

        # flight dynamics TODO
        self.altitude = altitude * 1000  # m
        self.alt_array = np.arange(0, self.altitude+1)
        self.longitude = longitude_0
        self.latitude = latitude_0
        self.V_c_0 = np.sqrt(self.g_0 * (self.R_e + self.altitude))
        if entry_speed > self.V_c_0:
            self.entry_speed = self.V_c_0  # m/s
        else:
            self.entry_speed = entry_speed
#         self.entry_flight_path_angle = flight_path_angle  # deg
#         self.entry_flight_path_angle_rad = np.radians(self.entry_flight_path_angle)  # rad

        # geometry
        self.mass = vehicle_mass
        self.S = vehicle_surface
        self.nose_radius = nose_radius

        # aerodynamics
        self.C_D = C_D
        self.C_L = C_L

        if boundary_layer == "laminar":
            self.n = 0.5
        elif boundary_layer == "turbulent":
            self.n = 0.2
        else:
            raise ValueError("boundary_layer must be either 'laminar' or 'turbulent'")
#
        if lift_drag_ratio is None:
            self.lift_drag_ratio = self.C_L / self.C_D
        else:
            self.lift_drag_ratio = lift_drag_ratio
#
        # CALCULATIONS

        self.V_c = self.circular_velocity()

        # atmosphere TODO
        self.atmosphere_model = atmosphere
        self.atmosphere = self.get_atmosphere()


        # lift parameter
        if lift_parameter is None:
            self.lift_parameter = (self.mass * self.g_0) / (self.S * self.C_D)
        else:
            self.lift_parameter = lift_parameter
#
        # velocity ratio
        self.v_ratio, self.v, self.v_min = self.velocity()
#
        # flight path angle
        self.flight_path_angle_equilibrium_rad = self.equilibrium_glide_angle()
        self.flight_path_angle_equilibrium_deg = self.flight_path_angle_equilibrium_rad * 180 / np.pi

        # flight range ratio
        self.flight_range_ratio = self.get_flight_range_ratio()

        # flight time
        self.t_h, self.t_min, self.t_sec = self.get_flight_time()                         # s

        # DECELERATION
        # maximum deceleration TODO
        self.a, self.a_max, self.a_ratio = self.get_gliding_deceleration()                                # g's
        self.v_ratio_a_max, self.V_a_max = self.get_max_deceleration_v_ratio()
        self.rho_a_max, self.h_a_max = self.get_max_deceleration_altitude()

        # self.deceleration_dataframe = pd.DataFrame({"loads": self.a,
        #                                             ""})

        # heat flux TODO
        self.qc, self.qc_max, self.q_normalized = self.get_stagnation_heatflux()

        # dataframe TODO
        self.dataframe = pd.DataFrame({
                                       "normalized_velocity_ratio": self.v_ratio,
                                        "equilibrium_glide_angle": self.flight_path_angle_equilibrium_deg,
                                        "loads": self.a,
                                        "deceleration": self.a*9.81,
                                        "deceleration_ratio": self.a / self.a_max,
                                        "heatflux_ratio": self.q_normalized,
                                        "heatflux": self.qc,
                                         "v": self.v
                                         }, index = [self.alt_array/1000])
    def get_constants(self):
        # TESTED
        if self.planet == "Earth" or self.planet == "earth":
            g = ct.g_0
            rho_0 = 1.225
            R = ct.earth_radius
            scale_height = 7200
            m = 3
        elif self.planet == "Mars" or self.planet == "mars":
            g = 3.711
            rho_0 = 0.699
            R = 3389500
            scale_height = 11000
            m = 0
        elif self.planet == "Jupiter" or self.planet == "jupiter":
            g = 24.79
            rho_0 = 1.267
            R = 69911000
            scale_height = 27000
            m = 0
        elif self.planet == "Saturn" or self.planet == "saturn":
            g = 10.44
            rho_0 = 0.622
            R = 58232000
            scale_height = 59500
            m = 0
        else:
            raise ValueError("planet must be either 'earth', 'mars', 'jupiter', 'saturn'")
        return g, rho_0, R, scale_height, m

    def equations_of_motion(self, h, V, gamma):
        # calculations
        rho = atm.ExponentialAtmosphere(h).rho_exp
        circular_velocity = np.sqrt(self.g_0 * (self.R_e + h * 1000))

        # equations
        dV_dt = (- 0.5 * self.C_D * rho * V**2 * self.S - self.mass * self.g_0 * np.sin(gamma)) / self.mass
        dgamma_dt = (0.5 * self.C_L * rho * V**2 * self.S - self.mass * self.g_0 * np.sin(gamma) * (1 - (V/circular_velocity)**2))/(self.mass * V)
        dh_dt = V * np.sin(gamma)

    def get_atmosphere(self):
        # TESTED
        dfs = []
        for i in self.alt_array:
            height = i / 1000
            if self.atmosphere_model == "exponential":
                atmos = atm.ExponentialAtmosphere(height)
                density = atmos.rho_exp
            elif self.atmosphere_model == "standard":
                atmos = atm.StandardAtmosphere(height)
                density = atmos.rho_std
            elif self.atmosphere_model == "nrlmsise-00":
                atmos = atm.NrlmsiseAtmosphere(height,
                                               time = "2025-07-01 22:18:45",
                                               latitude = self.latitude,
                                               longitude = self.longitude)
                density = atmos.rho_nrl
            else:
                raise ValueError("atmosphere_model must be either 'exponential' or 'standard' or 'nrlmsise-00'")
            df = pd.DataFrame({'density': density}, index = [height])
            dfs.append(df)

        density_dataframe = pd.concat(dfs)

        return density_dataframe


    def circular_velocity(self):
        # TESTED
        circular_velocity = self.V_c_0 * np.sqrt(1 / (1 + self.altitude /  self.R_e))

        return circular_velocity

    def velocity(self):
        # TESTED
        # normalized velocity ratio : with respect to circular velocity
        normalised_velocity_ratio = np.sqrt(self.lift_parameter / (0.5 * self.rho_0 * np.exp(-self.beta * self.alt_array) * self.V_c_0 ** 2 + self.lift_parameter))
        v = normalised_velocity_ratio * self.V_c_0
        v_min = (1 / self.V_c_0) * np.sqrt(2 * self.lift_parameter / self.rho_0) * self.V_c_0

        return normalised_velocity_ratio, v, v_min

    def equilibrium_glide_angle(self):
        # TESTED
        flight_path_angle = -(1 / (self.beta * self.R_e)) * (2 / self.lift_drag_ratio) * self.v_ratio ** (-2)

        return flight_path_angle

    def get_flight_range_ratio(self):
        # TESTED
        entry_circular_ratio = self.entry_speed / self.V_c_0
        Rf_Re = -0.5 * self.lift_drag_ratio * np.log(1 - (entry_circular_ratio ** 2))

        return Rf_Re

    def get_flight_time(self):
        # TESTED

        entry_circular_ratio = self.entry_speed / self.V_c_0
        t_flight_sec = ((0.5 * self.V_c_0 * self.lift_drag_ratio) / self.g_0) * np.log(
            (1 + entry_circular_ratio) / (1 - entry_circular_ratio))
        t_flight_min = t_flight_sec / 60
        t_flight_h = t_flight_min / 60

        return t_flight_h, t_flight_min, t_flight_sec

    def get_gliding_deceleration(self):
        # CORRECT
        a = (1/self.lift_drag_ratio) * (1 - (self.v_ratio ** 2))
                                        #- (2 / (self.beta * self.R_e * self.v_ratio ** 2)))

        v_ratio = (2 / (self.beta * self.R_e))**(1/4)
        # a_max_1 = self.lift_drag_ratio ** (-1) * (1 - 2 * np.sqrt(2 / (self.beta * self.R_e)))
        a_max = (1/self.lift_drag_ratio) * (1 - (v_ratio ** 2) )
                                            #- (2 / (self.beta * self.R_e * v_ratio ** 2)))

        a_ratio = a/a_max



        return a, a_max, a_ratio

    def get_max_deceleration_v_ratio(self):
        v_ratio = (2 / (self.beta * self.R_e))**(1/4)
        v_max_dec = v_ratio * self.V_c_0

        return v_ratio, v_max_dec

    def get_max_deceleration_altitude(self):
        rho = 2 * self.lift_parameter*((1/self.V_a_max**2)- (1/self.V_c_0**2))
        max_deceleration_altitude = (1 / self.beta) * np.log(self.rho_0/rho)

        return rho, max_deceleration_altitude
#
    def get_stagnation_heatflux(self):
        c1 = self.c_star * (1 / np.sqrt(self.rho_0)) * (1 / self.V_c_0 ** 3)
        c2 = c1 * (1 / self.nose_radius ** self.n) * (2 * self.lift_parameter / self.V_c_0 ** 2) ** (
                    1 - self.n) * self.V_c_0 ** self.m

        rho_max = 4 * (self.lift_parameter / (self.V_c_0 ** 2)) * ((1 - self.n) / (self.m + 2 * (self.n - 1)))
        V_max = np.sqrt((self.m - 2 * (1 - self.n)) / self.m) * self.V_c_0
        qc_max = c1 * (1 / self.nose_radius ** self.n) * rho_max ** (1 - self.n) * V_max ** self.m
        qc = c2 * (((1 / self.v_ratio ** 2) - 1) ** (
                    1 - self.n)) * self.v_ratio ** self.m
        normalised_heat_flux = qc/ qc_max


        return qc*1e-3, qc_max*1e-3, normalised_heat_flux

glide = GlidingEntry("earth",
                 100,
                 7700,
                 27000,
                 80,
                 5,
                 1.6,
                 0.205,
                 0,
                 0,
                 lift_parameter = None,
                 lift_drag_ratio = 0.13,
                 atmosphere = "exponential",
                 boundary_layer = "laminar")

print(glide.a_max, glide.qc_max, glide.V_a_max, glide.flight_range_ratio*glide.R_e)
