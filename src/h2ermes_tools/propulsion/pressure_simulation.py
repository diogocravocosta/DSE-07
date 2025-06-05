import matplotlib.pyplot as plt


#data
velocity_propellant = 10 #in m/s
chamber_pressure = 6100000 #in Pa
tank_pressure = 25000 #in Pa, lowest pressure at tank
# Propellant properties
# Propellant density in kg/m^3
propellant_density = 423



def pressure_profile_propulsion(propellant_density, velocity_propellant, chamber_pressure, tank_pressure):
    # From tank to pump inlet
    dynamic_pressure_increase = 0.5 * propellant_density * velocity_propellant**2  # in Pa
    line_losses = 50000  # conservative estimate, in Pa
    pump_inlet_pressure = tank_pressure - dynamic_pressure_increase - line_losses  # in Pa

    if pump_inlet_pressure < 0:
        pump_inlet_pressure = 0

    # From chamber pressure to pump discharge pressure 
    injector_losses = 0.3 * chamber_pressure  # throttled injector assumption
    pump_discharge_pressure = chamber_pressure + injector_losses + line_losses  # in Pa

    required_pump_pressure_rise = (pump_discharge_pressure - pump_inlet_pressure) * 1.1  # in Pa

    return required_pump_pressure_rise  # in Pa


if __name__ == "__main__":
    required_pump_pressure_rise = pressure_profile_propulsion(propellant_density, velocity_propellant, chamber_pressure, tank_pressure)
    print(f"Required pump pressure rise: {required_pump_pressure_rise / 1e6} MPa")
    # Output the required pump pressure rise in MPa

