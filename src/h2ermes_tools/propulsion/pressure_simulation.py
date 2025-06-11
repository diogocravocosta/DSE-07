import matplotlib.pyplot as plt


#data
velocity_propellant = 10 #in m/s
chamber_pressure = 6100000 #in Pa
prop_tank_pressure = 300000 #in Pa, lowest pressure at tank
# Propellant properties
# Propellant density in kg/m^3
propellant_density = 70.85



def pressure_profile_propulsion(propellant_density, velocity_propellant, chamber_pressure, prop_tank_pressure):
    # From tank to pump inlet
    dynamic_pressure_increase = 0.5 * propellant_density * velocity_propellant**2 
    line_losses = 50000  # conservative estimate, in Pa
    pump_inlet_pressure = prop_tank_pressure - dynamic_pressure_increase - line_losses  

    if pump_inlet_pressure < 0:
        pump_inlet_pressure = 0

    # From chamber pressure to pump discharge pressure 
    injector_losses = 0.3 * chamber_pressure  # throttled injector assumption
    pump_discharge_pressure = chamber_pressure + injector_losses + line_losses  

    required_pump_pressure_rise = (pump_discharge_pressure - pump_inlet_pressure) * 1.1  

    return pump_discharge_pressure, required_pump_pressure_rise, pump_inlet_pressure


if __name__ == "__main__":
    required_pump_pressure_rise = pressure_profile_propulsion(propellant_density, velocity_propellant, chamber_pressure, prop_tank_pressure)
    print(f"Required pump pressure rise: {required_pump_pressure_rise / 1e6} MPa")
    # Output the required pump pressure rise in MPa
    print(f"Pump inlet pressure: {required_pump_pressure_rise[2] / 1e6} MPa")

