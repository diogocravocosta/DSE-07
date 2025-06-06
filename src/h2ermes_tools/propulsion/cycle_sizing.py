import numpy as np
from h2ermes_tools.propulsion.pressure_simulation import pressure_profile_propulsion

velocity_propellant = 10 #in m/s
chamber_pressure = 6100000 #in Pa
tank_pressure = 25000 #in Pa, lowest pressure at tank
propellant_density = 423
required_pump_pressure_rise = pressure_profile_propulsion(propellant_density, velocity_propellant, chamber_pressure, tank_pressure)
allowable_pressure_rise_per_stage = 16000000 # in Pa, maximum pressure rise per stage of the pump
mass_flow = 387.8498585536046


def sizing_pump(required_pump_pressure_rise, propellant_vapor_pressure, pump_inlet_pressure, propellant_density):

    g_0 = 9.81  # gravitational acceleration in m/s^2
    volume_flow_rate = mass_flow / propellant_density  # in m^3/s
    net_positive_suction_head = (pump_inlet_pressure - propellant_vapor_pressure) / (g_0 * propellant_density)
    pump_head_pressure_rise = required_pump_pressure_rise / (g_0 * propellant_density) # in Pa

    number_of_stages = required_pump_pressure_rise / allowable_pressure_rise_per_stage + 1
    stage_specific_speed = 2.0 #value for liquid hydrogen
    pump_efficiency = 0.79
    pump_rotational_speed = (stage_specific_speed * (pump_head_pressure_rise/number_of_stages)**0.75)/np.sqrt(volume_flow_rate)
    power_required_pump = (g_0*mass_flow*pump_head_pressure_rise)/pump_efficiency
    return number_of_stages, pump_rotational_speed, stage_specific_speed, power_required_pump


if __name__ == "__main__":
    number_of_stages, pump_rotational_speed, stage_specific_speed, power_required_pump = sizing_pump(required_pump_pressure_rise, propellant_vapor_pressure=0, pump_inlet_pressure=tank_pressure, propellant_density=propellant_density)
    print(f"Number of stages: {number_of_stages}")
    print(f"Pump rotational speed: {pump_rotational_speed} rad/s")
    print(f"Stage specific speed: {stage_specific_speed}")
    print(f"Power required for pump: {power_required_pump} W")  # Output in MW



def sizing_turbine():
    gamma = 
    turbine_pressure_ratio =  1.75 # pressure ratio across the turbine
    isentropic_spouting_velocity = np.sqrt(2*c_p*turbine_inlet_temperature*(1-(1/turbine_pressure_ratio)**(gamma/(gamma-1))))

    A = 1.5 #empirical coefficient
    B = 0.6 #empirical exponent
    turbopump_mass = A * (pump_shaft_torque ** B)
    