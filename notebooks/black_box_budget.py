import numpy as np

mass_LH2 = 10000 # mass of LH2 in kg
transfer_time = 12 # time to transfer in hours
transfer_time = transfer_time * 3600 # convert to seconds
mass_flow = mass_LH2 / transfer_time # mass flow in kg/s
pipe_material_density = 8000 # density of pipe material in kg/m^3ç, stainless steel was used
flow_velocity = 10 # flow velocity in m/s
density_LH2 = 70.9 # density of LH2 in kg/m^3, source: https://demaco-cryogenics.com/blog/energy-density-of-hydrogen/#:~:text=In%20liquid%20form%20and%20at,density%20of%2070.9%20kg%2Fm%C2%B3.

#add margins
mass_flow_collision = mass_flow * 1.2 # 20% margin for collision


area = mass_flow_collision / (density_LH2*flow_velocity) # area in m^2
mass = area * pipe_material_density # mass of the pipe in kg
diameter = (4*area/np.pi)**0.5 # diameter in m
final_diameter = diameter * 2 # due to second wall and vacuum insulation
final_area = (4*final_diameter/np.pi)**0.5 # final area in m^2

print('Area of the pipe:', final_area, 'm^2')
print('Diameter of the pipe:', final_diameter, 'm') 