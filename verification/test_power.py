from src.h2ermes_tools.power_mass import power_plant_mass, lb_to_kg

shuttle_power = 14000 # kW

# Nominal Shuttle consumption is 8 tanks of oxygen and hydrogen over 8 days
print(power_plant_mass(average_power=shuttle_power, mission_length_days=8, safety_factor=1))
# Mass of 2 power cells and contents of 3 oxygen and hydrogen tanks
print(2 * lb_to_kg(255) + 3 * lb_to_kg(781 + 92))