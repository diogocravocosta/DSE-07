class H2ermes:
    """
    Class to initiate the state of the H2ermes vehicle at the beginning of re-entry

    """
    def __init__(self,
                 mass: float,
                 surface_area: float,
                 length: float,
                 drag_coefficient: float,
                 lift_coefficient: float,
                 initial_altitude: float,
                 initial_velocity: float,
                 final_velocity: float,
                 flight_path_angle: float,
                 angle_of_attack:float,
                 heading_angle: float = 0,
                 bank_angle: float = 110,
                 side_slip_angle: float = 0,
                 angular_rate: float = 0,
                 initial_longitude: float = -90,
                 initial_latitude: float = 0,
                 target_longitude: float = 0,
                 target_latitude: float = 0,



                 ):
        # INITIAL STATES
        self.t_0 = 0 # initial time
        self.h_0 = initial_altitude # initial altitude
        self.V_0 = initial_velocity # initial velocity

        self.long_0 = initial_longitude # initial longitude
        self.lat_0 = initial_latitude # initial latitude
        self.gamma_0 = flight_path_angle # initial flight path angle
        self.xi_0 = heading_angle # initial heading angle
        self.sigma_0 = bank_angle # initial bank angle
        self.alpha_0 = angle_of_attack # initial angle of attack
        self.beta_0 = side_slip_angle # initial side-slip angle
        self.angular_rate_0 = angular_rate # initial angular rate angle

        # VEHICLE DETAILS
        self.mass = mass # vehicle mass
        self.J = 0 # moment of inertia tensor #TODO: write code to calculate
        self.r_cg = 0 # COM vector #TODO: write code to calculate
        self.S_ref = surface_area # reference area
        self.L_ref = length # reference length
        self.C_D = drag_coefficient  # drag coefficient
        self.C_L = lift_coefficient  # lift coefficient

        # FINAL STATES
        self.long_t = target_longitude # target longitude
        self.lat_t = target_latitude # target latitude
        self.V_f = final_velocity # final velocity
