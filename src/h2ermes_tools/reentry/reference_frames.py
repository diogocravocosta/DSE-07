import numpy as np
from astropy.units import deg_C


# # internal
# import data.constants as ct
# from h2ermes_tools.reentry import atmosphere as atm

class ReferenceFrame:
    def __init__(self, frames, frame_dicts):

        self.frame_1 = frames[0]
        self.frame_2 = frames[1]

        self.frame_1_dict = frame_dicts[0]
        self.frame_2_dict = frame_dicts[1]

        # aerodynamics - airspeed (AA)
        if self.frame_1 == "aero_a" or self.frame_2 == "aero_a":
            if self.frame_1 == "aero_a":
                aero_a_dictionary = self.frame_1_dict
            else:
                aero_a_dictionary = self.frame_2_dict
            self.sigma_a = aero_a_dictionary["bank_angle"]
            self.alpha_a = aero_a_dictionary["angle_of_attack"]
            self.beta_a = aero_a_dictionary["angle_of_sideslip"]

        # aerodynamic - groundspeed (AG)
        if self.frame_1 == "aero_g" or self.frame_2 == "aero_g":
            if self.frame_1 == "aero_g":
                aero_g_dictionary = self.frame_1_dict
            else:
                aero_g_dictionary = self.frame_2_dict
            self.alpha_g = aero_g_dictionary["angle_of_attack"]
            self.beta_g = aero_g_dictionary["angle_of_sideslip"]
            self.sigma_g = aero_g_dictionary["bank_angle"]

        # body (B)

        # inertial planetocentric (I)

        if self.frame_1 == "inertial" or self.frame_2 == "inertial":
            if self.frame_1 == "inertial":
                inertial_dictionary = self.frame_1_dict
            else:
                inertial_dictionary = self.frame_2_dict

            self.x_I, self.y_I, self.z_I = inertial_dictionary["x"], inertial_dictionary["y"], inertial_dictionary["z"]
            self.x_dot_I, self.y_dot_I, self.z_dot_I = inertial_dictionary["x_dot"], inertial_dictionary["y_dot"], inertial_dictionary["z_dot"]

        # orbital (O) & local orbital (L)
        # ORBITAL
        if self.frame_1 == "orbital" or self.frame_2 == "orbital":
            if self.frame_1 == "orbital":
                orbital_dictionary = self.frame_1_dict
            else:
                orbital_dictionary = self.frame_2_dict
            self.e = orbital_dictionary["eccentricity"]  # eccentricity, (0 ≤ e < 1)
            self.a = orbital_dictionary["semi_major_axis"]  # semi_major_axis (a > Re)
            self.i = orbital_dictionary["inclination"]  # inclination (0deg ≤ i ≤ 180deg)
            self.omega_pericenter = orbital_dictionary["pericenter_argument"]  # pericenter argument (0deg ≤ omega_pericenter ≤ 360deg)
            self.omega_asc_node = orbital_dictionary["ascending_node_longitude"]  # ascending node longitude (0deg ≤ omega_asc_node ≤ 360deg)
            self.M = orbital_dictionary["mean_anomaly"]  # mean anomaly (0deg ≤ M ≤ 360deg)
            self.E = orbital_dictionary["eccentric_anomaly"]  # eccentric anomaly
            self.tau_o = orbital_dictionary["pericenter_passage_time"]  # time of pericenter passage
            self.t_0 = orbital_dictionary["epoch_time"]  # epoch
            self.t_o = 0
            self.theta = orbital_dictionary["true_anomaly"] # true anomaly

            self.mu = orbital_dictionary["gravitation_parameter"]  # gravitation parameter

            self.n = np.sqrt(self.mu ** 3 / self.a)

            self.M_0 = self.n * (self.t_0 - self.tau_o)

            self.M = self.M_0 + self.n * (self.t_o - self.t_0)

            self.M = self.E - self.e * np.sin(self.E)

        # rotating planetocentric (R)
        if self.frame_1 == "rotating" or self.frame_2 == "rotating":
            if self.frame_1 == "rotating":
                rotating_dictionary = self.frame_1_dict
            else:
                rotating_dictionary = self.frame_2_dict
            self.omega_cb = rotating_dictionary["rotational_rate_central_body"]
            self.t_r = rotating_dictionary["epoch_time"]




        # vertical (V)
        if self.frame_1 == "vertical" or self.frame_2 == "vertical":
            if self.frame_1 == "vertical":
                vertical_dictionary = self.frame_1_dict
            else:
                vertical_dictionary = self.frame_2_dict

            self.tau_v = vertical_dictionary["longitude"]
            self.delta_v = vertical_dictionary["latitude"]

        # trajectory - groundspeed (TG)
        if self.frame_1 == "traj_g" or self.frame_2 == "traj_g":
            if self.frame_1 == "traj_g":
                traj_g_dictionary = self.frame_1_dict
            else:
                traj_g_dictionary = self.frame_2_dict

            self.gamma_g = traj_g_dictionary["flight_path_angle"]
            self.chi_g = traj_g_dictionary["heading"]

        # trajectory - airspeed (TA)
        if self.frame_1 == "traj_a" or self.frame_2 == "traj_a":
            if self.frame_1 == "traj_a":
                traj_a_dictionary = self.frame_1_dict
            else:
                traj_a_dictionary = self.frame_2_dict

            self.gamma_a = traj_a_dictionary["flight_path_angle"]
            self.chi_a = traj_a_dictionary["heading"]




        # wind (W)
        if self.frame_1 == "wind" or self.frame_2 == "wind":
            if self.frame_1 == "wind":
                wind_dictionary = self.frame_1_dict
            else:
                wind_dictionary = self.frame_2_dict

            self.gamma_w = wind_dictionary["flight_path_angle"]
            self.chi_w = wind_dictionary["heading"]

        # propulsion (P)
        if self.frame_1 == "propulsion" or self.frame_2 == "propulsion":
            if self.frame_1 == "propulsion":
                propulsion_dictionary = self.frame_1_dict
            else:
                propulsion_dictionary = self.frame_2_dict
            self.epsilon_p = propulsion_dictionary["elevation_angle_thrust"]
            self.psi_p = propulsion_dictionary["azimuth_angle_thrust"]

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

        # ANGULAR RATES
        p = 0 # roll rate
        q = 0 # pitch rate
        r = 0 # yaw rate


    # MATRIX CALCULATIONS
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

    # FRAME TRANSFORMATIONS
    def aerodynamic_to_trajectory_a(self):
        # AA to TA
        C_TAAA = self.C_x(self.sigma_a)

        return C_TAAA

    def aerodynamic_to_vertical(self):
        # AA to V
        C_VAA = self.C_z(-self.chi_a) * self.C_y(-self.gamma_a) * self.C_x(self.sigma_a)

        return C_VAA

    def aerodynamic_to_trajectory_g(self):
        # AA to TG
        C_TGAA = self.C_y(self.gamma_g)*self.C_z(self.chi_g)*self.C_z(-self.chi_a)* self.C_y(-self.gamma_a)*self.C_x(self.sigma_a)

        return C_TGAA

    def inertial_to_orbital(self):
        # I to L
        C_LI = self.C_z((-self.omega_pericenter+self.theta)) * self.C_x(self.i) * self.C_z(self.omega_asc_node)

        # I to O
        a = np.array([np.cos(self.omega_asc_node), np.sin(self.omega_asc_node),0])

        C_OI = np.array([[np.cos(self.i)+((np.cos(self.omega_asc_node))**2)*(1-np.cos(self.i)),
                            np.cos(self.omega_asc_node)*np.sin(self.omega_asc_node)*(1-np.cos(self.i)),
                            -np.sin(self.omega_asc_node)*np.sin(self.i)],
                         [np.cos(self.omega_asc_node)*np.sin(self.omega_asc_node)*(1-np.cos(self.i)),
                            np.cos(self.i)+((np.sin(self.omega_asc_node))**2)*(1-np.cos(self.i)),
                            np.cos(self.omega_asc_node)*np.sin(self.i)],
                         [np.sin(self.omega_asc_node)*np.sin(self.i),
                            -np.cos(self.omega_asc_node)*np.sin(self.i),
                            np.cos(self.i)]])
        lbd = self.omega_asc_node + (self.omega_pericenter + self.theta)

        # L to O
        C_LO = np.array([[np.cos(lbd), np.sin(lbd), 0],
                         [-np.sin(lbd), np.cos(lbd), 0],
                         [0, 0, 1]])

        return C_LI, C_OI, C_LO

    def rotating_to_inertial(self):
        # R to I
        C_IR = self.C_z(self.omega_cb * self.t_r)

        return C_IR

    def vertical_to_rotating(self):
        # V to R
        C_RV = self.C_z(-self.tau_v) * self.C_y((np.pi/2) + self.delta_v)

        return C_RV

    def wind_to_vertical(self):
        # W to V
        C_VW = self.C_z(-self.chi_w) * self.C_y(-self.gamma_w)

        return C_VW

    def trajectory_g_to_vertical(self):
        # TG to V
        C_VTG = self.C_z(-self.chi_g) * self.C_y(-self.gamma_g)

        return C_VTG

    def trajectory_a_to_vertical(self):
        # TA to V
        C_VTA = self.C_z(-self.chi_a) * self.C_y(-self.gamma_a)

        return C_VTA



    def body_to_aerodynamic_air(self):
        # B to AA
        C_AAB = self.C_z(self.beta_a) * self.C_y(-self.alpha_a)

        return C_AAB

    def body_to_aerodynamic_ground(self):
        # B to AG
        C_AGB = self.C_z(self.beta_g) * self.C_y(-self.alpha_g)

        return C_AGB

    def body_to_propulsion(self):
        # B to P
        C_PB = self.C_y(self.epsilon_p) * self.C_z(self.psi_p)

        return C_PB

    def body_to_trajectory(self):
        if self.frame_1 == "traj_a" or self.frame_2 == "traj_a":
            alpha = self.alpha_a
            beta = self.beta_a
            sigma = self.sigma_a

        else:
            alpha = self.alpha_g
            beta = self.beta_g
            sigma = self.sigma_g

        C_TB = self.C_x(sigma) * self.C_z(beta) * self.C_y(-alpha)

        return C_TB

    def vertical_to_inertial(self):
        # V to I
        C_IV = self.C_z(self.omega_cb * self.t_r) * self.C_y((np.pi / 2) + self.delta_v)

        return C_IV



        














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







