import numpy as np
import matplotlib.pyplot as plt
from data import constants as cs
from data import material as mat
import rainflow
from scipy.stats import linregress
import fatigue_calculation_tool as fct

phi = 6
                           
tank_radius = 5
thickness_tank = 0.004
material = mat.Material(
        density=7850,
        youngs_modulus=200e9,
        thermal_expansion_coefficient=1.5e-6,  # Changed from thermal_expansion_coefficient
        fracture_toughness=200,
        yield_strength=500e6)

fuel_reentry_LH2 = 3000
dry_mass = 20000
launch_mass = 220000
payload_mass = 15500
g_reentry_force_ratio = 6
g_launch_force_ratio = 6
max_thrust2weight=4.3


thickness_minimum_tank = fct.thickness_optimization_fatigue(phi,
                                    tank_radius,
                                    thickness_tank,
                                    material,
                                    fuel_reentry_LH2,
                                    dry_mass,
                                    launch_mass,
                                    payload_mass,
                                    g_reentry_force_ratio,
                                    g_launch_force_ratio,
                                    max_thrust2weight,
                                    const_miner = 188400,
                                    exp_coeff_miner = -0.155,
                                    number_of_launches=40,
                                    min_launches=25,
                                    safety_factor=2,
                                    time_mission=[0, 0.1, 0.3, 18, 21, 23, 24],
                                    a_crack_depth=0.001,
                                    pressure_launch_tank=3e5,
                                    pressure_dock_vent=5e5,
                                    pressure_vent=1e6,
                                    pressure_after_dock=4.5e5,
                                    T_lh2=20,
                                    T_space=4,
                                    T_ambient_earth=300,
                                    T_gh2=150,
                                    T_gh2_ext=200)

print("Min thickness of tank:", thickness_minimum_tank,"m")
