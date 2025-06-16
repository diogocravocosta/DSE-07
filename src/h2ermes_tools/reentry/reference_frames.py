import numpy as np
from astropy.units import deg_C


# # internal
# import data.constants as ct
# from h2ermes_tools.reentry import atmosphere as atm

class ReferenceFrame:
    def __init__(self, frames, frame_1_dict, frame_2_dict):

        frame_1 = frames[0]
        frame_2 = frames[1]

        # inertial planetocentric (I)

        if frame_1 == "inertial" or frame_2 == "inertial":
            x_I, y_I, z_I = 0, 0, 0
            x_dot_I, y_dot_I, z_dot_I = 0, 0, 0

        # orbital (O) & local orbital (L)
        # ORBITAL
        if frame_1 == "orbital" or frame_2 == "orbital":
            if frame_1 == "orbital":
                orbital_dictionary = frame_1_dict
            else:
                orbital_dictionary = frame_2_dict
            self.e = orbital_dictionary["eccentricity"]  # eccentricity, (0 ≤ e < 1)
            self.a = orbital_dictionary["semi_major_axis"]  # semi_major_axis (a > Re)
            self.i = orbital_dictionary["inclination"]  # inclination (0deg ≤ i ≤ 180deg)
            self.omega_pericenter = orbital_dictionary["pericenter_argument"]  # pericenter argument (0deg ≤ omega_pericenter ≤ 360deg)
            self.omega_asc_node = orbital_dictionary["ascending_node_longitude"]  # ascending node longitude (0deg ≤ omega_asc_node ≤ 360deg)
            self.M = orbital_dictionary["mean_anomaly"]  # mean anomaly (0deg ≤ M ≤ 360deg)
            self.E = orbital_dictionary["eccentric_anomaly"]  # eccentric anomaly
            self.tau = orbital_dictionary["pericenter_passage_time"]  # time of pericenter passage
            self.t_0 = orbital_dictionary["epoch_time"]  # epoch
            self.t = 0
            self.theta = orbital_dictionary["true_anomaly"] # true anomaly

            self.mu = orbital_dictionary["gravitation_parameter"]  # gravitation parameter

            self.n = np.sqrt(self.mu ** 3 / self.a)

            self.M_0 = self.n * (self.t_0 - self.tau)

            self.M = self.M_0 + self.n * (self.t - self.t_0)

            self.M = self.E - self.e * np.sin(self.E)

        # rotating planetocentric (R)

        # body (B)

        # vertical (V)

        # trajectory - groundspeed (TG)

        # trajectory - airspeed (TA)

        # aerodynamic - groundspeed (AG)

        # aerodynamics - airspeed (AA)

        # wind (W)

        # propulsion (P)

    def position_and_velocity(self):

        # rotational
        x_R, y_R, z_R = 0, 0, 0
        u_R, v_R, w_R = 0, 0, 0

        # spherical
        R = 0 # distance
        tau = 0 # longitude (0deg ≤ tau < 360deg)
        delta = 0 # latitude (-90deg ≤ delta < 90deg)
        V_g = 0 # groundspeed
        gamma_g = 0 # flight path angle
        chi_g = 0 # heading

    def attitude(self):
        # CLASSICAL
        phi = 0 # roll angle
        theta = 0 # pitch angle
        psi = 0 # yaw angle

        # AERODYNAMIC
        alpha = 0 # angle of attack (–180deg ≤ alpha < 180deg)
        beta = 0 # angle of sideslip (–90deg ≤ beta < 90deg)
        sigma = 0 # bank angle (–180deg ≤ sigma < 180deg)

        # ANGULAR RATES
        p = 0 # roll rate
        q = 0 # pitch rate
        r = 0 # yaw rate

    def C_x(self, alpha):
        C_x = np.array([[1, 0, 0],
                        [0, np.cos(alpha), np.sin(alpha)],
                        [0, -np.sin(alpha), np.cos(alpha)]])

        return C_x

    def C_y(self, alpha):
        C_y = np.array([[np.cos(alpha), 0, -np.sin(alpha)],
                        [0, 1, 0],
                        [np.sin(alpha), 0, np.cos(alpha)]])

        return C_y

    def C_z(self, alpha):
        C_z = np.array([[np.cos(alpha), np.sin(alpha), 0],
                        [-np.sin(alpha), np.cos(alpha), 0],
                        [0, 0, 1]])

        return C_z



    def rotate_unit_axis(self, alpha_1, alpha_2, alpha_3):
        # transformation matrix
        C_BA = self.C_y(alpha_2) * self.C_z(-alpha_3) * self.C_x(-alpha_1)

        # inverse
        C_AB = self.C_x(alpha_1) * self.C_z(alpha_3) * self.C_y(-alpha_2)

    def inertial_to_orbital(self):
        # I to L
        C_LI = self.C_z((-self.omega_pericenter+self.theta)) * self.C_x(self.i) * self.C_z(self.omega_asc_node)

        # I to O
        a = np.array([np.cos(self.omega_asc_node), np.sin(self.omega_asc_node),0])

        C_OI = np.array([[np.cos(self.i)+((np.cos(self.omega_asc_node))**2)*(1-np.cos(self.i)),
                          np.cos(self.omega_asc_node)*np.sin(self.omega_asc_node)*(1-np.cos(self.i)),
                          ]])







orbital_dict = {"eccentricity": 0,
                "semi_major_axis":0,
                "inclination": 0,
                "pericenter_argument": 0,
                "ascending_node_longitude": 0,
                "mean_anomaly": 0,
                "eccentric_anomaly": 0,
                "pericenter_passage_time": 0,
                "epoch_time": 0,
                "true_anomaly": 0,
                "gravitation_parameter": 0,
                }

# inertial_dict = {}
#
# ref = ReferenceFrame(["orbital", "inertial"],
#                      orbital_dict,
#                      inertial_dict)
angle_deg = int(90)
angle = np.radians(angle_deg)

cosine = np.cos(angle)
print(cosine)

calc = np.cos(angle)*np.sin(angle)*(1-np.cos(angle))
print(calc)







