import numpy as np
from h2ermes_tools.propulsion.pressure_simulation import pressure_profile_propulsion
from pyfluids import Fluid, FluidsList, Input

velocity_propellant = 10 #in m/s
chamber_pressure = 6100000 #in Pa
tank_pressure = 300000 #in Pa, lowest pressure at tank
hydrogen = Fluid(FluidsList.Hydrogen).with_state(
    Input.temperature(77), Input.pressure(300000)
)
propellant_density = 70.85
x, required_pump_pressure_rise, pump_inlet_pressure = pressure_profile_propulsion(propellant_density, velocity_propellant, chamber_pressure, tank_pressure)
allowable_pressure_rise_per_stage = 16000000 # in Pa, maximum pressure rise per stage of the pump
mass_flow = 514.49 
propellant_vapor_pressure = 133322


def sizing_pump(required_pump_pressure_rise, propellant_vapor_pressure, pump_inlet_pressure, propellant_density, mass_flow):
    print('propellant density is', propellant_density)
    g_0 = 9.81  # gravitational acceleration in m/s^2
    volume_flow_rate = mass_flow / propellant_density  # in m^3/s
    print ('pump inlet pressure is ', pump_inlet_pressure, 'Pa')
    print('propellant vapor pressure is', propellant_vapor_pressure, 'Pa')
    net_positive_suction_head = (pump_inlet_pressure - propellant_vapor_pressure) / (g_0 * propellant_density)
    pump_head_pressure_rise = required_pump_pressure_rise / (g_0 * propellant_density) # in m
    print(pump_head_pressure_rise, ' pump head pressure rise in m')
    print(required_pump_pressure_rise, ' required pump pressure rise in Pa')
    number_of_stages = (required_pump_pressure_rise / allowable_pressure_rise_per_stage) 
    number_of_stages = int(np.ceil(number_of_stages)) + 1 # round up to the nearest whole number
    print(number_of_stages, 'number of stages')
    stage_specific_speed = 2.0 #value for liquid hydrogen
    pump_efficiency = 0.79
    #pump_rotational_speed = (stage_specific_speed * (pump_head_pressure_rise/number_of_stages)**0.75)/np.sqrt(volume_flow_rate)
    pump_rotational_speed = (130 * net_positive_suction_head**0.75) / np.sqrt(volume_flow_rate)
    power_required_pump = (g_0*mass_flow*pump_head_pressure_rise)/pump_efficiency
    pump_shaft_torque = power_required_pump / pump_rotational_speed  
    print(pump_shaft_torque, ' pump shaft torque in Nm')
    A = 1.5 #empirical coefficient
    B = 0.6 #empirical exponent
    turbopump_mass = A * (pump_shaft_torque ** B)
    
    turbine_pressure_ratio =  1.75 # pressure ratio across the turbine
    #isentropic_spouting_velocity = np.sqrt(2*c_p*turbine_inlet_temperature*(1-(1/turbine_pressure_ratio)**(gamma/(gamma-1))))
    #turbine_mean_pitch_diameter = (2 * u_m) / pump_rotational_speed
    return number_of_stages, pump_rotational_speed, stage_specific_speed, power_required_pump, net_positive_suction_head, turbopump_mass


if __name__ == "__main__":
    number_of_stages, pump_rotational_speed, stage_specific_speed, power_required_pump, net_positive_suction_head, turbopump_mass = sizing_pump(required_pump_pressure_rise, propellant_vapor_pressure, pump_inlet_pressure, propellant_density, mass_flow)
    print(f"Pump rotational speed: {pump_rotational_speed} rad/s")
    print(f"Stage specific speed: {stage_specific_speed}")
    print(f"Power required for pump: {power_required_pump} W")  # Output in MW
    print(f"Net positive suction head: {net_positive_suction_head} m")
    print(f"Turbopump mass: {turbopump_mass} kg")  # Output in kg

    