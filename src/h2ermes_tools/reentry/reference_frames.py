from scipy.optimize import dual_annealing
import data.constants as ct
import numpy as np
from pyswarm import pso
from scipy.optimize import fsolve

# Define the symbols
scale_height = 7200
g = ct.g_0
earth_radius = ct.earth_radius
lift_drag_ratio = 0.13
v_entry = 7974.762214
beta = 1/scale_height
entry_angle = -0.023604017
entry_height = 100000
drag_coefficient = 1.6
surface = np.pi * 4.92**2
mass = 33000
density_0 = 1.225
# Define the constants and parameters
constant = g * earth_radius
v_entry_ratio = v_entry/np.sqrt(constant)

c1 = 0.5*lift_drag_ratio
c2 = beta*earth_radius*(np.cos(entry_angle)-1)
print(v_entry_ratio)

# Define the objective function
# Define the functions for numerical solving
def equations(vars):
    v_final_ratio = vars
    eq1 = -entry_angle - (c1 * np.log((v_entry_ratio**2) / (v_final_ratio**2))) / (((1/v_final_ratio**2 - 1) / (c2) - 1))
    # eq2 = -np.exp(-beta * height_final) + np.exp(-beta * entry_height) - (1 + (1 / (beta * earth_radius)) * ((1/v_final_ratio**2) - 1) - np.cos(entry_angle)) / (0.5 * lift_drag_ratio * (drag_coefficient * surface / mass * beta) * density_0)
    return eq1

# Initial guesses for v_final and height_final
initial_guess = [v_entry_ratio]

# Solve the equations numerically
solutions = fsolve(equations, initial_guess)

print(solutions)
