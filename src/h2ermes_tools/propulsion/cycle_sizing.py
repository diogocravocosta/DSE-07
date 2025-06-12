import numpy as np
from h2ermes_tools.propulsion.pressure_simulation import pressure_profile_propulsion
from h2ermes_tools.propulsion.general_characteristics_calculation import calculate_mass_flow


def sizing_pump(required_pump_pressure_rise,
                propellant_vapor_pressure,
                pump_inlet_pressure,
                propellant_density,
                mass_flow,
                allowable_pressure_rise_per_stage):
    # print('propellant density is', propellant_density)
    g_0 = 9.81  # gravitational acceleration in m/s^2
    volume_flow_rate = mass_flow / propellant_density  # in m^3/s
    # print ('pump inlet pressure is ', pump_inlet_pressure, 'Pa')
    # print('propellant vapor pressure is', propellant_vapor_pressure, 'Pa')
    net_positive_suction_head = (pump_inlet_pressure - propellant_vapor_pressure) / (g_0 * propellant_density)
    pump_head_pressure_rise = required_pump_pressure_rise / (g_0 * propellant_density) # in m
    # print(pump_head_pressure_rise, ' pump head pressure rise in m')
    # print(required_pump_pressure_rise, ' required pump pressure rise in Pa')
    number_of_stages = (required_pump_pressure_rise / allowable_pressure_rise_per_stage) 
    number_of_stages = int(np.ceil(number_of_stages)) + 1 # round up to the nearest whole number
    # print(number_of_stages, 'number of stages')
    stage_specific_speed = 2.0 #value for liquid hydrogen
    pump_efficiency = 0.79
    #pump_rotational_speed = (stage_specific_speed * (pump_head_pressure_rise/number_of_stages)**0.75)/np.sqrt(volume_flow_rate)
    pump_rotational_speed = (130 * net_positive_suction_head**0.75) / np.sqrt(volume_flow_rate)
    power_required_pump = (g_0*mass_flow*pump_head_pressure_rise)/pump_efficiency
    pump_shaft_torque = power_required_pump / pump_rotational_speed  
    # print(pump_shaft_torque, ' pump shaft torque in Nm')
    A = 1.5 #empirical coefficient
    B = 0.6 #empirical exponent
    turbopump_mass = A * (pump_shaft_torque ** B)
    
    turbine_pressure_ratio =  1.75 # pressure ratio across the turbine
    #isentropic_spouting_velocity = np.sqrt(2*c_p*turbine_inlet_temperature*(1-(1/turbine_pressure_ratio)**(gamma/(gamma-1))))
    #turbine_mean_pitch_diameter = (2 * u_m) / pump_rotational_speed
    return number_of_stages, pump_rotational_speed, stage_specific_speed, power_required_pump, net_positive_suction_head, turbopump_mass

def size_turbopump(
        tank_pressure: float,
        thrust: float,
        propellant_density: float = 70.85,
        velocity_propellant: float = 10,
        chamber_pressure: float = 300000,
        allowable_pressure_rise_per_stage: float = 16000000,
        propellant_vapor_pressure: float = 133322
        ):
    """
    Wrapper for sizing the turbopumps
    """
    (x,
     required_pump_pressure_rise,
     pump_inlet_pressure
     ) = pressure_profile_propulsion(
        propellant_density,
        velocity_propellant,
        chamber_pressure,
        tank_pressure
    )
    mass_flow = calculate_mass_flow(thrust)

    (
        number_of_stages,
        pump_rotational_speed,
        stage_specific_speed,
        power_required_pump,
        net_positive_suction_head,
        turbopump_mass
    ) = sizing_pump(
        required_pump_pressure_rise=required_pump_pressure_rise,
        propellant_vapor_pressure=propellant_vapor_pressure,
        pump_inlet_pressure=pump_inlet_pressure,
        propellant_density=propellant_density,
        mass_flow=mass_flow,
        allowable_pressure_rise_per_stage=allowable_pressure_rise_per_stage
    )

    # print(f'Rotational Speed: {pump_rotational_speed} rad/s\n'
    #       f'Stage Specific Speed: {stage_specific_speed}\n'
    #       f'Power Required: {power_required_pump} W\n'
    #       f'Net Positive Suction Head: {net_positive_suction_head} m\n'
    #       f'Turbopump Mass: {turbopump_mass} kg')

    return turbopump_mass


if __name__ == "__main__":
    size_turbopump(3e5, 2.17e6)

    