import numpy as np
import matplotlib.pyplot as plt

# High power battery
# source: Burke
battery_P = {
    # Specific power density
    'rho_W' : 3000,
    # Specific energy density
    'rho_Wh' : 100,
}

# High energy battery
# source: Burke
battery_E = {
    # Specific power density
    'rho_W' : 800,
    # Specific energy density
    'rho_Wh' : 170,
}

# Fuel cells
# source: P&P reader, Shuttle cells
fuel_cell = {
    # Specific power density of cell
    'rho_W' : 132, # 132/103.7
    # Specific energy density of reactant
    'rho_Wh' : 3000,
}

# Solar panel array
# At beginning of life, 28% efficiency
# source: New SMAD
# Deployment penalty is an estimate
solar_panel = {
    # Mass per panel area
    'rho_A' : 2.8, # kg/m^2
    # Power per panel area
    'rho_W' : 245,  # W/m^2
    # Mass of panel deployment/retraction system (necessary for reentry)
    # Defined by multiplying panel mass with this factor
    'deployment_system_penalty' : 2,
}

def size_fuel_cell(max_power, avg_power, mission_length, fuel_cell):
    # Assumptions:
    # - Power output linearly scales fuel cell weight
    # - Energy required linearly scales reactant mass
    m_cell = max_power / fuel_cell['rho_W']
    m_reactant = avg_power * mission_length / fuel_cell['rho_Wh']
    m_cell = np.full_like(m_reactant, m_cell) # Same shape as m_reactant
    return m_cell, m_reactant

def size_battery(max_power, avg_power, mission_length, battery):
    # Battery mass driven by power
    m_battery_W = max_power / battery['rho_W']
    # Battery mass driven by energy
    m_battery_E = avg_power * mission_length / battery['rho_Wh']
    m_battery_W = np.full_like(m_battery_E, m_battery_W)
    # Return the larger one of the two masses
    return np.maximum(m_battery_E, m_battery_W)

def size_solar_panels(max_power, avg_power, solar_panel, battery, h):
    # Assumptions:
    # - no solar panel power degradation (only on longer timescales)
    T_orbit = orbit_period(h) / 3600 # in hours
    T_eclipse = time_in_eclipse(h) / 3600 # in hours
    # Energy the battery needs to provide during eclipse
    energy_in_eclipse = T_eclipse * avg_power
    # Time spent in sun
    T_sun = T_orbit - T_eclipse
    m_battery_orbit = size_battery(max_power, avg_power, T_eclipse, battery) # Eclipse time
    m_battery_launch = size_battery(max_power, avg_power, 2, battery) # 45 min orbital insertion, 1 hour on ground
    m_battery = max(m_battery_orbit, m_battery_launch) # Take the larger of the two

    # Power needed to charge the battery during sun
    battery_charging_power = energy_in_eclipse / T_sun
    total_power = avg_power + battery_charging_power
    sp_area = total_power / solar_panel['rho_W']
    print(f"Solar panel area: {sp_area:.2f} m^2")
    sp_mass = sp_area * solar_panel['rho_A']

    return sp_mass * solar_panel['deployment_system_penalty'], m_battery

def orbit_period(h):
    # Orbital period for circular orbit around Earth at altitude h
    r_E = 6378 # km
    mu = 3.986e5 # km3/s2
    a = r_E + h # Semi-major axis
    T = 2 * np.pi * np.sqrt(a ** 3 / mu)
    return T

def time_in_eclipse(h):
    # Time spent in eclipse for circular orbit around Earth at altitude h, assuming perpendicular Sun rays
    r_E = 6378
    a = r_E + h
    # Angle corresponding to the part of orbit in eclipse
    l = 2 * np.arcsin(r_E / a)
    # Ratio of eclipse to full orbit
    orbit_eclipse_percent = l / (2 * np.pi)
    return orbit_eclipse_percent * orbit_period(h)

# Tests
#print(orbit_period(430)/60) # ISS
#print(time_in_eclipse(430)/60)

# Define mission parameters
avg_power = 14e3 # Shuttle parameters, no payload
max_power = 21e3
# New SMAD says 2-3x average power
mission_length = np.geomspace(1/4, 14 * 24, 5000) # in hours, 15 min to 5 days
h = 600 # altitude of depot

specific_mission_length = 5 * 24 # 3 days nominal, 5 days contingency

# Size fuel cells
w_fuel_cell, w_fuel_cell_reactant = size_fuel_cell(max_power, avg_power, mission_length, fuel_cell)

# Size batteries
w_battery_P = size_battery(max_power, avg_power, mission_length, battery_P)
w_battery_E = size_battery(max_power, avg_power, mission_length, battery_E)

# Size solar panels
w_solar_panel_P, w_sp_battery_P = size_solar_panels(max_power, avg_power, solar_panel, battery_P, h)
w_solar_panel_E, w_sp_battery_E = size_solar_panels(max_power, avg_power, solar_panel, battery_E, h)

# Reshape to match mission length
w_solar_panel_P = np.full_like(mission_length, w_solar_panel_P)
w_solar_panel_E = np.full_like(mission_length, w_solar_panel_E)
w_sp_battery_P = np.full_like(mission_length, w_sp_battery_P)
w_sp_battery_E = np.full_like(mission_length, w_sp_battery_E)

# Calculate total masses
w_fuel_cell_total = w_fuel_cell + w_fuel_cell_reactant
w_solar_panel_P_total = w_solar_panel_P + w_sp_battery_P
w_solar_panel_E_total = w_solar_panel_E + w_sp_battery_E

# Find the minimum mass system for all mission lengths
w_min = np.minimum.reduce([w_fuel_cell_total,
                           w_battery_E,
                           w_battery_P,
                           w_solar_panel_E_total,
                           w_solar_panel_P_total])

# If a system is minimum mass for a certain mission length, extract those masses from the arrays
w_min_fuel_cell = np.where(w_min == w_fuel_cell_total, w_fuel_cell_total, np.nan)
w_min_battery_E = np.where(w_min == w_battery_E, w_battery_E, np.nan)
w_min_battery_P = np.where(w_min == w_battery_P, w_battery_P, np.nan)
w_min_solar_panel_E = np.where(w_min == w_solar_panel_E_total, w_solar_panel_E_total, np.nan)
w_min_solar_panel_P = np.where(w_min == w_solar_panel_P_total, w_solar_panel_P_total, np.nan)


y1 = w_min_battery_E[~np.isnan(w_min_battery_E)][-1]
y2 = w_min_fuel_cell[~np.isnan(w_min_fuel_cell)][-1]

x1 = mission_length[w_min_battery_E == y1][0]
x2 = mission_length[w_min_fuel_cell == y2][0]

w0 = size_battery(max_power, avg_power, specific_mission_length, battery_E)
print(f"Battery mass at 5 days: {w0:.1f}")

w1, w2 = size_fuel_cell(max_power, avg_power, specific_mission_length, fuel_cell)
print(f"Fuel cell mass at 5 days: {w1:.1f}, reactant mass {w2:.1f}, total {w1 + w2:.1f}")

w3, w4 = size_solar_panels(max_power, avg_power, solar_panel, battery_E, h)
print(f"Solar panel mass at 5 days: {w3:.1f}, battery mass: {w4:.1f}, total {w3 + w4:.1f}")

# Plot, with mission length on log x-axis
plt.figure()
plt.plot(mission_length, w_fuel_cell, label='Fuel cell', color='cadetblue')
plt.plot(mission_length, w_fuel_cell_reactant, label='Fuel cell reactant', color='royalblue', alpha=0.7)
plt.plot(mission_length, w_fuel_cell_total, label='Fuel cell total', color='deepskyblue', alpha=0.7)
plt.plot(mission_length, w_min_fuel_cell, color='deepskyblue', linewidth=3, linestyle='dashed')

plt.plot(mission_length, w_battery_P, label='High power Li-ion battery', color='limegreen')
plt.plot(mission_length, w_min_battery_P, color='blue', linewidth=3)
plt.plot(mission_length, w_battery_E, label='High energy Li-ion battery', color='green')
plt.plot(mission_length, w_min_battery_E, color='green', linewidth=3, linestyle='dashed')

#plt.plot(mission_length, w_solar_panel_P, label='Solar panel using high power Li-ion battery')
#plt.plot(mission_length, w_solar_panel_E, label='Solar panel using high energy Li-ion battery')
#plt.plot(mission_length, w_sp_battery_P, label='High power Li-ion battery of solar panel')
#plt.plot(mission_length, w_sp_battery_E, label='High energy Li-ion battery of solar panel')
#plt.plot(mission_length, w_solar_panel_P + w_sp_battery_P, label='Solar panel (high power Li-ion)')
plt.plot(mission_length, w_solar_panel_E_total, label='Solar panel (with high energy Li-ion)', color='orange')
plt.plot(mission_length, w_min_solar_panel_E, color='orange', linewidth=3, linestyle='dashed')

# Plot crossover points
plt.scatter(x1, y1, color='red', linewidths=3, zorder=10)
plt.scatter(x2, y2, color='red', linewidths=3, zorder=10)
# Annotate crossover points, put annotation below the point
plt.annotate(f'{x1:.1f} hours', xy=(x1, y1), xytext=(2.5, -15), textcoords='offset points', color='red')
plt.annotate(f'{x2/24:.1f} days', xy=(x2, y2), xytext=(2.5, -15), textcoords='offset points', color='red')

# Plot a vertical line at the mission length of 5 days
plt.axvline(x=specific_mission_length, color='gray', linestyle='dotted')

#plt.plot(mission_length, w_min, label='Minimum mass system', color='black', linewidth=2)
plt.xscale('log')
#plt.yscale('log')
# Set max y-axis to 1000 kg
plt.ylim(0, 1000)
plt.xlabel('Mission length (h)')
plt.ylabel('Mass (kg)')
plt.legend()
plt.grid()
plt.show()
