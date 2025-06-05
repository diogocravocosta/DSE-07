"""This file contains all constants used throughout our codebase"""

# General physics
gravitational_constant = 6.67430e-11  # N*m2/kg2
speed_of_light = 299792458  # m/s, speed of light in vacuum
magnetic_permeability_constant = 7.8 * 10**15 # Teslas * m**3 # magnetic constant (vacuum permeability)

# Earth parameters
g_0 = 9.80665  # m/s2, standard gravity
earth_mass = 5.9722e24  # kg
gravitational_parameter = gravitational_constant * earth_mass  # m3/s2
earth_radius = 6378137  # m
sidereal_rotation_period = 23 * 3600 + 56 * 60 + 4.1  # s, period of one Earths rotation around its axis with respect to distant stars
density_sea_level = 1.225  # kg/m3
orbit_altitude_density = 1.03 * 10**-14  # kg/m3, density at 600 km altitude
solar_constant = 1361  # W/m2, solar constant at an altitude of 600 km

# Atmospheric constants
R_star = 8314.32  # J/(kmol K), universal gas constant
molecular_mass = 28.9644  # kg/kmol, molecular mass of Earth's atmosphere
gamma = 1.4  # specific heat ratio
