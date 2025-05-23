import pandas as pd
import hvplot.pandas

from src.trajectories import ballistic, gliding, skipping

hvplot.extension('bokeh')

class Rocket:
    def __init__(self,
                 entry_speed,
                 surface,
                 h_0,
                 mass,
                 flight_path_angle,
                 lift_drag_ratio,
                 scale_height=7200,
                 nose_radius=7,
                 c_d=0.5,
                 c_l=1.0,
                 boundary_layer="laminar",
                 skip_number=5,
                 lift_parameter = None,
                 ballistic_coefficient=None):
        self.planet = "Earth"
        self.entry_speed = entry_speed
        self.surface = surface
        self.h_0 = h_0
        self.mass = mass
        self.flight_path_angle = flight_path_angle
        self.lift_drag_ratio = lift_drag_ratio
        self.scale_height = scale_height
        self.nose_radius = nose_radius
        self.c_d = c_d
        self.c_l = c_l
        self.boundary_layer = boundary_layer
        self.skip_number = skip_number
        self.lift_parameter = lift_parameter
        self.ballistic_coefficient = ballistic_coefficient

        self.skip = skipping.SkippingEntry(self.planet, self.entry_speed,self.flight_path_angle,self.h_0,
                 self.mass,
                 self.surface,
                 self.c_d,
                 self.c_l,
                 self.nose_radius,
                 self.lift_drag_ratio,
                 self.lift_parameter,
                 self.boundary_layer,
                 self.skip_number,
                 )
        self.glide = gliding.GlidingEntry(self.planet, self.entry_speed,self.flight_path_angle,self.h_0,
                 self.mass,
                 self.surface,
                 self.c_d,
                 self.c_l,
                 self.nose_radius,
                 self.lift_drag_ratio,
                 self.lift_parameter,
                 self.boundary_layer
                 )
        self.ballistic = ballistic.Ballistic(self.planet, self.entry_speed,self.flight_path_angle,self.h_0,
                 self.mass,
                 self.surface,
                 self.c_d,
                 self.nose_radius, self.ballistic_coefficient,self.boundary_layer
                 )
        # qcs
        self.skip_qc_max, self.glide_qc_max, self.ballistic_qc_max = self.get_qc_max()

    def get_qc_max(self):
        skip_qc_max = self.skip.qc_max
        glide_qc_max = self.glide.qc_max
        ballistic_qc_max = self.ballistic.q_c_max

        return skip_qc_max, glide_qc_max, ballistic_qc_max

rocket = Rocket(entry_speed=7700,surface=36,h_0=120,mass=10000,flight_path_angle=-10,lift_drag_ratio=2,scale_height=7200,nose_radius=7,c_d=0.5,c_l=1.0,boundary_layer="laminar",skip_number=1, lift_parameter=None)

print(rocket.skip_qc_max, rocket.glide_qc_max, rocket.ballistic_qc_max)

