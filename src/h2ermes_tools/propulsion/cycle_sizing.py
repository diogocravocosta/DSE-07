import numpy as np
from h2ermes_tools.propulsion.pressure_simulation import pressure_profile_propulsion
from h2ermes_tools.propulsion.general_characteristics_calculation import mass_flow_rate
from pyfluids import Fluid, FluidsList, Input


chamber_pressure = 6100000 #in Pa
#tank_pressure = 270000 #in Pa, lowest pressure at tank
prop_tank_pressure = 270000
oxidizer_tank_pressure = 200000


propellant_density = 70.85
oxidizer_density = 1141

mass_flow = 514.49



def sizing_pump(prop_tank_pressure, oxidizer_tank_pressure, maximum_thrust, specific_impulse, chamber_pressure=61e5):

    mass_flow = mass_flow_rate(maximum_thrust, specific_impulse)
    #constants
    g_0 = 9.81  # gravitational acceleration in m/s^2
    velocity_propellant = 10 #in m/s
    velocity_oxidizer = 10 #m/s
    propellant_density = 70.85
    oxidizer_density = 1141

    x, required_pump_pressure_rise_prop, pump_inlet_pressure_prop = pressure_profile_propulsion(propellant_density, velocity_propellant, chamber_pressure, prop_tank_pressure)
    y, required_pump_pressure_rise_oxidizer, pump_inlet_pressure_oxidizer = pressure_profile_propulsion(oxidizer_density, velocity_oxidizer, chamber_pressure, oxidizer_tank_pressure)
    print(pump_inlet_pressure_oxidizer)
    #SIZING FOR PROPELLANT TURBOPUMP
    propellant_vapor_pressure = 84000 #found for LH2 at 20K
    allowable_pressure_rise_per_stage_prop = 16000000 #for LH2 from SPAD
    mass_flow_propellant = mass_flow/7
    volume_flow_rate_propellant = mass_flow_propellant / propellant_density  # in m^3/s
    net_positive_suction_head_prop = (pump_inlet_pressure_prop - propellant_vapor_pressure) / (g_0 * propellant_density)
    pump_head_pressure_rise_prop = required_pump_pressure_rise_prop / (g_0 * propellant_density) # in m
    number_of_stages_prop = (required_pump_pressure_rise_prop / allowable_pressure_rise_per_stage_prop) 
    number_of_stages_prop = int(np.ceil(number_of_stages_prop)) + 1 # round up to the nearest whole number
    stage_specific_speed_prop = 2.0 #value for liquid hydrogen
    pump_efficiency_prop = 0.79
    pump_rotational_speed_prop = (130 * net_positive_suction_head_prop**0.75) / np.sqrt(volume_flow_rate_propellant)
    power_required_pump_prop = (g_0*mass_flow_propellant*pump_head_pressure_rise_prop)/pump_efficiency_prop
    pump_shaft_torque_prop = power_required_pump_prop / pump_rotational_speed_prop  
    A = 1.5 #empirical coefficient
    B = 0.6 #empirical exponent
    turbopump_mass_prop = A * (pump_shaft_torque_prop ** B)


    #SIZING FOR OXIDIZER TURBOPUMP
    oxidizer_density = 1141
    oxidizer_vapor_pressure = 19707
    allowable_pressure_rise_per_stage_oxidizer = 47000000 #for all others from SPAD
    mass_flow_oxidizer = (mass_flow*6)/7
    volume_flow_rate_oxidizer= mass_flow_oxidizer / oxidizer_density  # in m^3/s
    net_positive_suction_head_oxidizer = (pump_inlet_pressure_oxidizer - oxidizer_vapor_pressure) / (g_0 * oxidizer_density)
    pump_head_pressure_rise_oxidizer = required_pump_pressure_rise_oxidizer / (g_0 * oxidizer_density) # in m
    number_of_stages_oxidizer = (required_pump_pressure_rise_oxidizer / allowable_pressure_rise_per_stage_oxidizer) 
    number_of_stages_oxidizer = int(np.ceil(number_of_stages_oxidizer)) + 1 # round up to the nearest whole number
    stage_specific_speed_oxidizer = 2.0 #value for liquid hydrogen
    pump_efficiency_oxidizer = 0.79
    pump_rotational_speed_oxidizer = (130 * net_positive_suction_head_oxidizer**0.75) / np.sqrt(volume_flow_rate_oxidizer)
    power_required_pump_oxidizer = (g_0*mass_flow_oxidizer*pump_head_pressure_rise_oxidizer)/pump_efficiency_oxidizer
    pump_shaft_torque_oxidizer = power_required_pump_oxidizer / pump_rotational_speed_oxidizer  
    A = 1.5 #empirical coefficient
    B = 0.6 #empirical exponent
    turbopump_mass_oxidizer = A * (pump_shaft_torque_oxidizer ** B)
    print(turbopump_mass_prop, turbopump_mass_oxidizer)

    

    mass_turbopumps = turbopump_mass_prop + turbopump_mass_oxidizer

    #turbine_mean_pitch_diameter = (2 * u_m) / pump_rotational_speed
    return mass_turbopumps

if __name__ == "__main__":
    mass_turbopumps = sizing_pump(mass_flow, chamber_pressure, prop_tank_pressure, oxidizer_tank_pressure)
    print(f"Mass of turbopumps: {mass_turbopumps:.2f} kg")
    # Output the mass of the turbopumps in kg

