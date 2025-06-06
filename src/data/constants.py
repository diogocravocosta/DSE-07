"""This file contains all constants used throughout our codebase"""

# General physics
gravitational_constant = 6.67430e-11  # N*m2/kg2

# Earth parameters
g_0 = 9.80665  # m/s2, standard gravity
earth_mass = 5.9722e24  # kg
gravitational_parameter = gravitational_constant * earth_mass  # m3/s2
earth_radius = 6378137  # m
sidereal_rotation_period = 23 * 3600 + 56 * 60 + 4.1  # s, period of one Earths rotation around its axis with respect to distant stars
density_sea_level = 1.225  # kg/m3

# Atmospheric constants
R_star = 8314.32  # J/(kmol K), universal gas constant
molecular_mass = 28.9644  # kg/kmol, molecular mass of Earth's atmosphere
gamma = 1.4  # specific heat ratio
