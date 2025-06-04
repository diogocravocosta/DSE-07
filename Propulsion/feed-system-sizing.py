import numpy as np

mass_flow = 15.2 #in kg/s
velocity_propellant = 10 #m/s
propellant_density = 70.85 #kg/m^3
head_loss_coefficient= 1.2 #takes non-reversible pressure losses into account, 1.2 for radiused edges of injector inlet, 1.7 for sharp edges

def obtain_size_propellant_channel(mass_flow, propellant_density, velocity_propellant, head_loss_coefficient): 
    propellant_channel_area = mass_flow/(propellant_density * velocity_propellant) #in m2
    propellant_channel_diameter = ((propellant_channel_area/np.pi) ** 0.5) * 2 #in m
    #pressure_drop_injector = head_loss_coefficient * 0.5 * propellant_density * velocity_propellant**2
    return propellant_channel_area, propellant_channel_diameter


propellant_channel_area, propellant_channel_diameter = obtain_size_propellant_channel(mass_flow, propellant_density, velocity_propellant, head_loss_coefficient)

# Print the results
print(f"Propellant Channel Area: {propellant_channel_area:.6f} m^2")
print(f"Propellant Channel Diameter: {propellant_channel_diameter:.6f} m")