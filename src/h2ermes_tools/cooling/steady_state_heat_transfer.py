"""
Steady-state 1-D heat conduction in a stainless-steel wall with:
• Constant incident heat flux at the hot face
• T⁴ radiation from the hot face (linearized for steady-state)
• Constant cooling flux at the cold face
Material: AISI 304L (properties assumed T-independent)
Scheme  : Finite difference, uniform grid
"""

import numpy as np
import matplotlib.pyplot as plt
from h2ermes_tools.cooling.material import SS310
import time

# ---------------- USER INPUT ----------------
wall_thickness = 4e-3  # [m]
num_nodes = 11
node_spacing = wall_thickness / (num_nodes - 1)

density = SS310.density  # [kg/m³]
specific_heat = SS310.specific_heat  # [J/(kg·K)]
thermal_conductivity = 25  # [W/(m·K)]

incident_heat_flux = 100_000.0  # [W/m²]
cooling_heat_flux = 50_000.0    # [W/m²]
emissivity = 0.9
stefan_boltzmann = 5.670374419e-8  # [W/(m²·K⁴)]

# For steady-state, we linearize the radiation term around an estimated surface temperature
# (or use a fixed-point iteration for more accuracy)
T_guess = 1000.0  # [K] initial guess for surface temperature

# Set up the coefficient matrix and right-hand side
A = np.zeros((num_nodes, num_nodes))
b = np.zeros(num_nodes)

# Internal nodes (finite difference: d²T/dx² = 0)
for i in range(1, num_nodes - 1):
    A[i, i - 1] = 1
    A[i, i] = -2
    A[i, i + 1] = 1
    b[i] = 0

# Boundary node 0 (hot face):
# -k (T0 - T1)/dx = q_incident - q_radiation
# Linearize q_radiation ≈ eps*sigma*4*T_guess^3*(T0 - T_guess) + eps*sigma*T_guess^4
rad_coeff = emissivity * stefan_boltzmann * 4 * T_guess**3
rad_const = emissivity * stefan_boltzmann * T_guess**4 - rad_coeff * T_guess
A[0, 0] = -thermal_conductivity / node_spacing - rad_coeff
A[0, 1] = thermal_conductivity / node_spacing
b[0] = -(incident_heat_flux - rad_const)

# Boundary node N-1 (cold face):
# -k (TN-1 - TN-2)/dx = -q_cooling
A[-1, -2] = -thermal_conductivity / node_spacing
A[-1, -1] = thermal_conductivity / node_spacing
b[-1] = -cooling_heat_flux

start_time = time.time()

# Solve the linear system
T = np.linalg.solve(A, b)

elapsed = time.time() - start_time

# Plot the result
# x = np.linspace(0, wall_thickness, num_nodes)
plt.figure(figsize=(9, 6))
# plt.plot(x, T, marker='o', label='Steady-state')
# plt.xlabel("x [m]")
# plt.ylabel("Temperature [K]")
# plt.title(f"1-D {num_nodes}-node wall - steady-state solution")
# plt.grid(True)
# plt.legend()
# plt.tight_layout()
# plt.show()

print("Steady-state solution finished")
# Estimate physical time to steady state using the slowest (lowest) eigenmode
# For a slab, tau ≈ L^2 / (π^2 * alpha)
thermal_diffusivity = thermal_conductivity / (density * specific_heat)
slab_tau = wall_thickness**2 / (np.pi**2 * thermal_diffusivity)
print(f"Estimated physical time to steady state (1/e): {slab_tau:.2f} s")
print(f"Δx = {node_spacing * 1e3:5.3f} mm")
print(f"Surface temperature (node 0): {T[0]:.1f} K")
print(f"Centre temperature (node (N-1)/2):  {T[int((num_nodes - 1) / 2)]:.1f} K")
print(f"Cooled face temperature (node N-1): {T[num_nodes - 1]:.1f} K")
