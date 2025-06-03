"""This file contains all constants used throughout our codebase"""

# General physics
gravitational_constant = 6.67430e-11 #N*m2/kg2
boltzmann = 5.67e-8  # stefan boltzmann constant



# Earth parameters
g_0 = 9.81 # m/s2, standard gravity
earth_mass = 5.9722e24 # kg
gravitational_parameter = gravitational_constant * earth_mass # m3/s2
earth_radius = 6378137 # m
sidereal_rotation_period = 23 * 3600 + 56 * 60 + 4.1 # s, period of one Earths rotation around its axis with respect to distant stars
