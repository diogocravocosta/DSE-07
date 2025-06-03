# external
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from scipy.integrate import simpson
import hvplot.pandas

# internal
import src.data.constants as ct

class Constants:
    def __init__(self):
        self.g = ct.g_0       # m/s^2
        self.rho_0 = 1.225    # kg/m^3
        self.R = ct.radius    # meters
        self.Hs = 7200        # m

        self.beta = 1 / self.Hs
        self.c_star = 1.1097e8
        self.m = 3

class GlidingEntry:
    def __init__(self,
                 planet,
                 entry_speed,
                 flight_path_angle,
                 h_0,
                 m,
                 S,
                 c_d,
                 c_l,
                 nose_radius=7,
                 lift_drag_ratio = 2,
                 lift_parameter = None,
                 boundary_layer="laminar",
                 range_to_cover = None
                 ):
        # CONSTANT
        # planetary
        self.planet = planet
        self.constant = Constants()
        self.g = self.constant.g
        self.rho_0 = self.constant.rho_0
        self.Re = self.constant.R

        self.h_s = self.constant.Hs
        self.beta = self.constant.beta
        self.c_star = self.constant.c_star

        # flight dynamics
        self.entry_speed = entry_speed  # m/s
        self.entry_flight_path_angle = flight_path_angle  # deg
        self.entry_flight_path_angle_rad = np.radians(self.entry_flight_path_angle)  # rad
        self.h_0 = h_0  # km

        # geometry
        self.mass = m
        self.S = S
        self.nose_radius = nose_radius

        # aerodynamics
        self.c_d = c_d
        self.c_l = c_l


        if boundary_layer == "laminar":
            self.n = 0.5
        elif boundary_layer == "turbulent":
            self.n = 0.2
        else:
            raise ValueError("boundary_layer must be either 'laminar' or 'turbulent'")

        # CALCULATIONS
        self.altitude = np.linspace(self.h_0 * 1000, 0, 10000)       # m

        # circular velocity
        self.v_c = np.sqrt(self.g * self.Re)

        if range_to_cover is None:
            self.lift_drag_ratio = lift_drag_ratio
        else:
            self.lift_drag_ratio = -2 * ((range_to_cover/self.Re)/(np.log(1-(self.entry_speed**2/self.v_c**2))))

        # atmosphere
        self.rho = self.rho_0 * np.exp(-self.beta * self.altitude)

        # lift parameter
        if lift_parameter is None:
            self.lift_parameter = (self.mass * self.g) / (self.S * self.c_l)
        else:
            self.lift_parameter = lift_parameter

        # normalized velocity ratio
        self.normalized_v_ratio, self.v, self.v_min = self.get_velocity()

        # flight path angle
        self.flight_path_angle_equilibrium = self.get_flight_path_angle()

        # flight range ratio
        self.flight_range_ratio = self.get_flight_range()

        # flight time
        self.t = self.get_flight_time()                         # s

        # maximum deceleration
        self.a, self.a_max = self.get_max_gliding_deceleration()

        # heat flux
        self.qc, self.qc_max, self.q_normalized = self.get_stagnation_heatflux()

        # dataframe
        self.dataframe = pd.DataFrame({"altitude": self.altitude / 1000,
                                       "normalized_velocity_ratio": self.normalized_v_ratio,
                                       "loads": self.a,
                                       "deceleration": self.a*9.81,
                                       "heatflux ratio": self.q_normalized,
                                       "heatflux": self.qc,
                                       })

    def get_velocity(self):
        # CORRECT

        normalised_velocity_ratio = np.sqrt(self.lift_parameter / (
                    0.5 * self.rho_0 * np.exp(-self.beta * self.altitude) * self.v_c ** 2 + self.lift_parameter))
        v = normalised_velocity_ratio * self.v_c
        v_min = (1 / self.v_c) * np.sqrt(2 * self.lift_parameter / self.rho_0) * self.v_c

        return normalised_velocity_ratio, v, v_min

    def get_flight_path_angle(self):
        # CORRECT

        flight_path_angle = -(1 / (self.beta * self.Re)) * (2 / self.lift_drag_ratio) * self.normalized_v_ratio ** (-2)

        return flight_path_angle

    def get_flight_range(self):
        # CORRECT

        entry_circular_ratio = self.entry_speed / self.v_c
        Rf_Re = -0.5 * self.lift_drag_ratio * np.log(1 - (entry_circular_ratio ** 2))

        return Rf_Re

    def get_flight_time(self):
        # CORRECT

        entry_circular_ratio = self.entry_speed / self.v_c
        t_flight = ((0.5 * self.v_c * self.lift_drag_ratio) / (2 * self.g)) * np.log(
            (1 + entry_circular_ratio) / (1 - entry_circular_ratio))

        return t_flight/60

    def get_max_gliding_deceleration(self):
        # CORRECT
        #a = 1/self.lift_drag_ratio * (1 - (self.normalized_v_ratio ** 2) - (2 / (self.beta * self.Re * self.normalized_v_ratio ** 2)))
        a_max = self.lift_drag_ratio ** (-1) * (1 - 2 * np.sqrt(2 / (self.beta * self.Re)))

        a = self.lift_drag_ratio ** (-1) * (1-(self.normalized_v_ratio ** 2))

        return a, a_max

    def get_max_deceleration_altitude(self):
        v_ratio = (2 / (self.beta * self.Re))**(1/4)
        v_max_dec = v_ratio * self.v_c
        rho = 2 * self.lift_parameter*((1/v_max_dec**2)- 1/self.v_c**2)
        max_deceleration_altitude = (1 / self.beta) * np.log(-self.rho_0/rho)

        return max_deceleration_altitude

    def get_stagnation_heatflux(self):
        c1 = self.c_star * (1 / np.sqrt(self.rho_0)) * (1 / self.v_c ** 3)
        c2 = c1 * (1 / self.nose_radius ** self.n) * (2 * self.lift_parameter / self.v_c ** 2) ** (
                    1 - self.n) * self.v_c ** self.constant.m

        rho_max = 4 * (self.lift_parameter / (self.v_c ** 2)) * ((1 - self.n) / (self.constant.m + 2 * (self.n - 1)))
        V_max = np.sqrt((self.constant.m - 2 * (1 - self.n)) / self.constant.m) * self.v_c
        qc_max = c1 * (1 / self.nose_radius ** self.n) * (rho_max) ** (1 - self.n) * (V_max) ** self.constant.m
        qc = c2 * (((1 / self.normalized_v_ratio ** 2) - 1) ** (
                    1 - self.n)) * self.normalized_v_ratio ** self.constant.m
        normalised_heat_flux = qc/ qc_max


        return qc*1e-3, qc_max*1e-3, normalised_heat_flux


glide = GlidingEntry(planet="Earth", entry_speed=7200,flight_path_angle=-10,h_0=120,m=30000,S=9,c_d=1.5,c_l=0.5,lift_drag_ratio=0.3)


print(glide.dataframe, glide.qc_max)

plot = glide.dataframe.hvplot.line(x="normalized_velocity_ratio",y="heatflux ratio",title="Heat Flux",ylabel="kW/m^2",xlabel="V ratio",width=600,height=400)
hvplot.show(plot)