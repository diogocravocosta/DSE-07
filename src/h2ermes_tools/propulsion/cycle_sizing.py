import numpy as np
from pressure_simulation import pressure_profile_propulsion


required_pump_pressure_rise = 




#constant data
allowable_pressure_rise_per_stage = 16 # in MPa, maximum pressure rise per stage of the pump


def sizing_pump(required_pump_pressure_rise, propellant_vapor_pressure, pump_inlet_pressure):
    volume_flow_rate = mass_flow / propellant_density  # in m^3/s
    net_positive_suction_head = (pump_inlet_pressure - propellant_vapor_pressure) / (g_0 * propellant_density)
    pump_head_pressure_rise = required_pump_pressure_rise / (g_0 * propellant_density) # in Pa

    number_of_stages = required_pump_pressure_rise / allowable_pressure_rise_per_stage + 1
    pump_rotational_speed = (u_ss * (net_positive_suction_head ** 0.75)) / (np.sqrt(volume_flow_rate))
    stage_specific_speed = pump_rotational_speed * np.sqrt(volume_flow_rate) / ((pump_head_pressure_rise/number_of_stages)** 0.75)

