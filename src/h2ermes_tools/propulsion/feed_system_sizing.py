import numpy as np
from pyfluids import Fluid, FluidsList, Input
from h2ermes_tools.propulsion.pressure_simulation import pressure_profile_propulsion


head_loss_coefficient= 1.2 #takes non-reversible pressure losses into account, 1.2 for radiused edges of injector inlet, 1.7 for sharp edges

def obtain_cross_sectional_areas(total_mass_flow, prop_tank_pressure, oxidizer_tank_pressure, chamber_pressure=6100000): 
    #properties
    propellant_density = 65 #kg/m^3 #taken from BoilFast
    LOX = Fluid(FluidsList.Oxygen).with_state(
    Input.temperature(77), Input.pressure(300000)
)
    #Velocities
    velocity_propellant_tank_exit = 10 #m/s 
    velocity_oxidizer_tank_exit = 25 #m/s 
    velocity_propellant_turbopump_exit = 10 #m/s
    velocity_oxidizer_turbopump_exit = 25 #m/s

    #mass flow calculations based on O/F ratio of 6
    total_mass_flow_propellant = total_mass_flow * 1/7 
    mass_flow_propellant_thruster = total_mass_flow_propellant/24 #in kg/s
    total_mass_flow_oxidizer = total_mass_flow * 6/7
    mass_flow_oxidizer_thruster = total_mass_flow_propellant/24 #in kg/s


    #propellant and oxidizer channel areas, from fuel tank to turbopump
    main_propellant_channel_area = total_mass_flow_propellant/(propellant_density * velocity_propellant_tank_exit) #in m2
    main_propellant_channel_diameter = ((main_propellant_channel_area/np.pi) ** 0.5) * 2 #in m
    main_oxidizer_channel_area = total_mass_flow_oxidizer/(LOX.density * velocity_oxidizer_tank_exit) #in m2
    main_oxidizer_channel_diameter = ((main_oxidizer_channel_area/np.pi) ** 0.5) * 2 #in m

    #propellant and oxidizer channel area from turpopump to injector
    propellant_channel_area_to_thrusters  = mass_flow_propellant_thruster/(propellant_density * velocity_propellant_turbopump_exit) #in m2
    small_propellant_diameter = ((propellant_channel_area_to_thrusters/np.pi) ** 0.5) * 2 #in m
    oxidizer_channel_area_to_thrusters = mass_flow_oxidizer_thruster/(LOX.density * velocity_oxidizer_turbopump_exit) #in m2
    small_oxidizer_diameter = ((oxidizer_channel_area_to_thrusters/np.pi) ** 0.5) * 2 #in m
    
    # Iterative calculation of thickness using the Barlow equation
    def calculate_thickness(inner_diameter, pressure, stress): 
        tolerance = 1e-6  # Convergence tolerance
        thickness = 0.001  # Initial guess for thickness in meters
        while True:
            outer_diameter = inner_diameter + thickness
            new_thickness = (pressure * outer_diameter) / (2 * stress)
            if abs(new_thickness - thickness) < tolerance:
                break
            thickness = new_thickness
        return thickness

    small_pressure_prop, o = pressure_profile_propulsion(propellant_density, velocity_propellant_turbopump_exit, chamber_pressure, prop_tank_pressure)
    internal_pressure_small_oxidizer, u = pressure_profile_propulsion(LOX.density, velocity_oxidizer_turbopump_exit, chamber_pressure, oxidizer_tank_pressure)
    chamber_pressure = 6100000  # in Pa, chamber pressure
    thickness_big_propellant = calculate_thickness(main_propellant_channel_diameter, prop_tank_pressure, 210 * 10 **6)
    thickness_big_oxidizer = calculate_thickness(main_oxidizer_channel_diameter, oxidizer_tank_pressure, 210 * 10 **6)
    thickness_small_propellant = calculate_thickness(small_propellant_diameter, small_pressure_prop, 210 * 10 **6)
    thickness_small_oxidizer = calculate_thickness(small_oxidizer_diameter, internal_pressure_small_oxidizer, 210 * 10 **6)

    return (main_propellant_channel_area, main_propellant_channel_diameter, main_oxidizer_channel_area,
            main_oxidizer_channel_diameter, propellant_channel_area_to_thrusters, small_propellant_diameter,
            oxidizer_channel_area_to_thrusters, small_oxidizer_diameter, thickness_big_propellant,
            thickness_big_oxidizer, thickness_small_propellant, thickness_small_oxidizer)
    


total_mass_flow = 514.49
prop_tank_pressure = 300000
oxidizer_tank_pressure = 270000

main_propellant_channel_area, main_propellant_channel_diameter, main_oxidizer_channel_area, \
main_oxidizer_channel_diameter, propellant_channel_area_to_thrusters, small_propellant_diameter, \
oxidizer_channel_area_to_thrusters, small_oxidizer_diameter, thickness_big_propellant, \
thickness_big_oxidizer, thickness_small_propellant, thickness_small_oxidizer = obtain_cross_sectional_areas(
    total_mass_flow, prop_tank_pressure, oxidizer_tank_pressure, chamber_pressure=6100000)


if __name__ == "__main__":
    # Example usage
    # Print the results
    print(f"Channel area for propellant from tanks to turboprop: {main_propellant_channel_area:.6f} m^2")
    print(f"Channel diameter for propellant from tanks to turboprop: {main_propellant_channel_diameter:.6f} m")
    print(f"Channel area for oxidizer from tanks to turboprop: {main_oxidizer_channel_area:.6f} m^2")
    print(f"Channel diameter for oxidizer from tanks to turboprop: {main_oxidizer_channel_diameter:.6f} m")

    print(f"Channel area for propellant from turboprop to thrusters: {propellant_channel_area_to_thrusters:.6f} m^2")
    print(f"Channel diameter for propellant from turboprop to thrusters: {small_propellant_diameter:.6f} m")
    print(f"Channel area for oxidizer from turboprop to thrusters: {oxidizer_channel_area_to_thrusters:.6f} m^2")
    print(f"Channel diameter for oxidizer from turboprop to thrusters: {small_oxidizer_diameter:.6f} m")

    print(f"Thickness of propellant channel from tanks to turboprop: {thickness_big_propellant:.6f} m")
    print(f"Thickness of oxidizer channel from tanks to turboprop: {thickness_big_oxidizer:.6f} m")
    print(f"Thickness of propellant channel from turboprop to thrusters: {thickness_small_propellant:.6f} m")
    print(f"Thickness of oxidizer channel from turboprop to thrusters: {thickness_small_oxidizer:.6f} m")




