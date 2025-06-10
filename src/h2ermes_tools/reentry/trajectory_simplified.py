# external
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from scipy.integrate import simpson
import hvplot.pandas

# internal
import src.data.constants as ct
from h2ermes_tools import atmosphere as atm

class GlidingEntry:
    def __init__(self,
                 planet,
                 altitude,
                 vehicle_mass,
                 vehicle_surface,
                 C_D,
                 C_L,
#                 entry_speed,
#                 flight_path_angle,
#                 c_l,
#                 nose_radius=7,
#                 lift_drag_ratio = None,
#                 lift_parameter = None,
#                 boundary_layer="turbulent",
#                 range_to_cover = None
                 ):
        # ASSUMPTIONS: spherical non-rotating earth, atmosphere at rest, trajectory in flat plane, constant g
        # CONSTANTS
        # planetary
        self.planet = planet
        self.g, self.rho_0, self.R, self.scale_height, self.m = self.get_constants()

        self.beta = 1 / self.scale_height
        self.c_star = 1.1097e8

        # flight dynamics TODO
#         self.entry_speed = entry_speed  # m/s
#         self.entry_flight_path_angle = flight_path_angle  # deg
#         self.entry_flight_path_angle_rad = np.radians(self.entry_flight_path_angle)  # rad
        self.altitude = altitude * 1000  # m
#
        # geometry
        self.mass = vehicle_mass
        self.S = vehicle_surface
#         self.nose_radius = nose_radius

        # aerodynamics
        self.C_D = C_D
        self.C_L = C_L

#         if boundary_layer == "laminar":
#             self.n = 0.5
#         elif boundary_layer == "turbulent":
#             self.n = 0.2
#         else:
#             raise ValueError("boundary_layer must be either 'laminar' or 'turbulent'")
#
#         if lift_drag_ratio is None:
#             self.lift_drag_ratio = self.c_l / self.c_d
#         else:
#             self.lift_drag_ratio = lift_drag_ratio
#
        # CALCULATIONS
        self.circular_velocity_0 = np.sqrt(self.g * (self.R + self.altitude))
#         self.altitude = np.linspace(self.h_0 * 1000, 0, 10000)       # m

#
#         if range_to_cover is None:
#             self.lift_drag_ratio = self.lift_drag_ratio
#         else:
#             self.lift_drag_ratio = -2 * ((range_to_cover/self.Re)/(np.log(1-(self.entry_speed**2/self.v_c**2))))
#
        # atmosphere TODO
#         self.rho = self.rho_0 * np.exp(-self.beta * self.altitude)
#
        # lift parameter
#         if lift_parameter is None:
#             self.lift_parameter = (self.mass * self.g) / (self.S * self.c_l)
#         else:
#             self.lift_parameter = lift_parameter
#
        # normalized velocity ratio TODO
#         self.normalized_v_ratio, self.v, self.v_min = self.get_velocity()
#
        # flight path angle TODO
#         self.flight_path_angle_equilibrium = self.get_flight_path_angle()
#
        # flight range ratio TODO
#         self.flight_range_ratio = self.get_flight_range()
#
        # flight time TODO
#         self.t = self.get_flight_time()                         # s
#
        # maximum deceleration TODO
#         self.a, self.a_max = self.get_max_gliding_deceleration()
#
        # heat flux TODO
#         self.qc, self.qc_max, self.q_normalized = self.get_stagnation_heatflux()
#
        # dataframe TODO
#         self.dataframe = pd.DataFrame({"altitude": self.altitude / 1000,
#                                        "normalized_velocity_ratio": self.normalized_v_ratio * self.v_c,
#                                        "loads": self.a,
#                                        "deceleration": self.a*9.81,
#                                        "heatflux ratio": self.q_normalized,
#                                        "heatflux": self.qc,
#                                        "v": self.v
#                                        })
    def get_constants(self):
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
        rho = atm.Atmosphere(h).rho_exp
        circular_velocity = np.sqrt(self.g * (self.R + h * 1000))

        # equations
        dV_dt = (- 0.5 * self.C_D * rho * V**2 * self.S - self.mass * self.g * np.sin(gamma)) / self.mass
        dgamma_dt = (0.5 * self.C_L * rho * V**2 * self.S - self.mass * self.g * np.sin(gamma) * (1 - (V/circular_velocity)**2))/(self.mass * V)
        dh_dt = V * np.sin(gamma)


#     def get_velocity(self):
#         # CORRECT
#         v_c = np.sqrt(self.g * (self.Re + self.altitude))
#         normalised_velocity_ratio = np.sqrt(self.lift_parameter / (
#                     0.5 * self.rho_0 * np.exp(-self.beta * self.altitude) * v_c ** 2 + self.lift_parameter))
#         v = normalised_velocity_ratio * v_c
#         v_min = (1 / self.v_c) * np.sqrt(2 * self.lift_parameter / self.rho_0) * self.v_c
#
#         return normalised_velocity_ratio, v, v_min
#
#     def get_flight_path_angle(self):
#         # CORRECT
#
#         flight_path_angle = -(1 / (self.beta * self.Re)) * (2 / self.lift_drag_ratio) * self.normalized_v_ratio ** (-2)
#
#         return flight_path_angle
#
#     def get_flight_range(self):
#         # CORRECT
#
#         entry_circular_ratio = self.entry_speed / self.v_c
#         Rf_Re = -0.5 * self.lift_drag_ratio * np.log(1 - (entry_circular_ratio ** 2))
#
#         return Rf_Re
#
#     def get_flight_time(self):
#         # CORRECT
#
#         entry_circular_ratio = self.entry_speed / self.v_c
#         t_flight = ((0.5 * self.v_c * self.lift_drag_ratio) / (2 * self.g)) * np.log(
#             (1 + entry_circular_ratio) / (1 - entry_circular_ratio))
#
#         return t_flight/60
#
#     def get_max_gliding_deceleration(self):
#         # CORRECT
#         #a = 1/self.lift_drag_ratio * (1 - (self.normalized_v_ratio ** 2) - (2 / (self.beta * self.Re * self.normalized_v_ratio ** 2)))
#         a_max = self.lift_drag_ratio ** (-1) * (1 - 2 * np.sqrt(2 / (self.beta * self.Re)))
#
#         a = self.lift_drag_ratio ** (-1) * (1-(self.normalized_v_ratio ** 2))
#
#         return a, a_max
#
#     def get_max_deceleration_altitude(self):
#         v_ratio = (2 / (self.beta * self.Re))**(1/4)
#         v_max_dec = v_ratio * self.v_c
#         rho = 2 * self.lift_parameter*((1/v_max_dec**2)- 1/self.v_c**2)
#         max_deceleration_altitude = (1 / self.beta) * np.log(-self.rho_0/rho)
#
#         return max_deceleration_altitude
#
#     def get_stagnation_heatflux(self):
#         c1 = self.c_star * (1 / np.sqrt(self.rho_0)) * (1 / self.v_c ** 3)
#         c2 = c1 * (1 / self.nose_radius ** self.n) * (2 * self.lift_parameter / self.v_c ** 2) ** (
#                     1 - self.n) * self.v_c ** self.constant.m
#
#         rho_max = 4 * (self.lift_parameter / (self.v_c ** 2)) * ((1 - self.n) / (self.constant.m + 2 * (self.n - 1)))
#         V_max = np.sqrt((self.constant.m - 2 * (1 - self.n)) / self.constant.m) * self.v_c
#         qc_max = c1 * (1 / self.nose_radius ** self.n) * (rho_max) ** (1 - self.n) * (V_max) ** self.constant.m
#         qc = c2 * (((1 / self.normalized_v_ratio ** 2) - 1) ** (
#                     1 - self.n)) * self.normalized_v_ratio ** self.constant.m
#         normalised_heat_flux = qc/ qc_max
#
#
#         return qc*1e-3, qc_max*1e-3, normalised_heat_flux
#
#
# # surface = np.pi * 4.92**2
# #
# # # glide_1 = GlidingEntry(planet="Earth", entry_speed=7200,flight_path_angle=-1,h_0=100,m=30000,S=surface,c_d= 1.4, c_l =  range_to_cover=3124643)
# #
# # glide_2 = GlidingEntry(planet="Earth", entry_speed=7200,flight_path_angle=-1,h_0=100,m=20000,S=surface,c_d=1.3,c_l=0.3)
# #
# # print(glide_2.flight_range_ratio*glide_2.Re, glide_2.lift_drag_ratio, glide_2.lift_parameter, glide_2.qc_max, glide_2.a_max)
# #
# #
# # # glide_high = GlidingEntry(planet="Earth", entry_speed=7200,flight_path_angle=-1,h_0=120,m=30000,S=9,c_d=1.2,c_l=1.2, lift_drag_ratio=1)
# # #
# # plot_1 = glide_2.dataframe.hvplot(x="v", y="altitude")
# # # plot_2 = glide_2.dataframe.hvplot(x="v", y="altitude")
# # # plot = plot_1 * plot_2
# # hvplot.show(plot_1)
#
#
