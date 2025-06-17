# external
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from scipy.integrate import simpson
from scipy.optimize import fsolve
import sympy as sp

# internal
import data.constants as ct

# Define the symbols
scale_height = 7200
g = ct.g_0
earth_radius = ct.earth_radius
lift_drag_ratio = 0.13
# v_entry = 7974.762214
beta = 1/scale_height
entry_angle = 0.023604017
entry_height = 100000
drag_coefficient = 1.6
surface = np.pi * 4.92**2
mass = 39000
density_0 = 1.225
# Define the constants and parameters
v_circular = 7844.195565

v_entry_ratio = 1.01657
print("V_entry / V_circular = ",v_entry_ratio)

v_entry_ratio_sq = v_entry_ratio**2
print("V_entry_ratio_squared = ",v_entry_ratio_sq)

c1 = 0.5 * lift_drag_ratio
print("0.5 * lift_drag_ratio = ",c1)

cosine = np.cos(entry_angle)
print("cosine entry_angle = ",cosine)

cosine_minus_1 = cosine - 1
print("cosine - 1 = ",cosine_minus_1)

c2 = beta * earth_radius * cosine_minus_1
print("beta * earth_radius * cosine_minus_1 = ",c2)

beta_height_entry = beta * entry_height
print("beta * entry_height = ",beta_height_entry)

c3 = 1 /(beta * earth_radius)
print("1 / beta * earth_radius = ",c3)

c4 = drag_coefficient * surface / (mass * beta)
print("drag_coefficient * S / m * beta = ",c4)

c5 = c4 * density_0 * c1
print("denominator =  ",c5)



# v_final_ratio,height_final = sp.symbols(
#     'v_final height_final'
# )

def equation_1(v_final_ratio):
    eq_1 = entry_angle + (c1 * (np.log(v_entry_ratio**2 / v_final_ratio**2) /
                                   ((((1 / v_final_ratio**2) - 1) / c2) - 1)))

    return eq_1

# Define the equation
# equation_1 = sp.Eq(-entry_angle, c1 * (sp.log(v_entry_ratio**2 / v_final_ratio**2) /
#                                   ((((1 / v_final_ratio**2) - 1) / c2) - 1)))


# # Solve for v_final
# solutions = sp.solve(equation_1, v_final_ratio, warn = True)

# Display the solution
# print(solutions)

sol_1 = fsolve(equation_1, v_entry_ratio)

v_final_ratio = sol_1[0]
#
# def equation_2(beta_height_final):
#     eq_2 = np.exp(-beta_height_entry) - np.exp(-beta_height_final) + (1 + (c3 * ((1/v_final_ratio**2)-1)) - cosine) / c5
#     return eq_2
#
# sol_2 = fsolve(equation_2, beta_height_entry)

beta_height_final = sp.symbols(
    'beta_height_final')

equation_2 = sp.Eq(
    sp.exp(-beta_height_final) - sp.exp(-beta_height_entry),(1 + (c3 * ((1/v_final_ratio**2)-1)) - cosine) / c5)

sol_2 = sp.solve(equation_2, beta_height_final)

beta_height_final = np.float64(sol_2[0])
height_final = beta_height_final /  beta

N = c1 * c4 - c3 * (1/density_0 * np.exp(-beta_height_final))*((1/v_final_ratio**2)-1)





print(v_final_ratio, height_final, N)
