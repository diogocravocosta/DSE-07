"""
1-D transient heat conduction in a stainless-steel wall
• Constant incident heat flux at the hot face
• T⁴ radiation from the hot face
• Constant cooling flux at the cold face
Material: AISI 304L (properties assumed T-independent)
Scheme  : Explicit FTCS, uniform grid
"""

import numpy as np
import matplotlib.pyplot as plt

# ---------------- USER INPUT ----------------
# Geometry & discretization
wall_thickness = 0.010  # total wall thickness [m]
num_nodes = 3  # number of nodes
node_spacing = wall_thickness / (num_nodes - 1)  # node spacing [m]

# Material properties (304L, room-temperature values)
density = 8030.0  # [kg/m³]
specific_heat = 500.0  # [J/(kg·K)]
thermal_conductivity = 14.6  # [W/(m·K)]
thermal_diffusivity = thermal_conductivity / (density * specific_heat)

# Boundary heat fluxes ( + = heat INTO the wall, – = heat OUT )
incident_heat_flux = 100_000.0  # constant aerodynamic heating on node 0 [W/m²]
cooling_heat_flux = 30_000.0  # constant cooling on node N-1 [W/m²]

# Radiation
emissivity = 0.78  # typical oxidised 304L
stefan_boltzmann = 5.670374419e-8  # [W/(m²·K⁴)]

# Initial & run parameters
initial_temperature = 20.0  # uniform start temperature [K]
total_time = 900.0  # total simulated time [s]

# Stability-limited time-step (Fo = α dt / dx² ≤ 0.5)
max_dt = 0.5 * node_spacing ** 2 / thermal_diffusivity

time_step = 0.9 * max_dt  # 90% of the limit for safety
num_time_steps = int(np.ceil(total_time / time_step))
time_axis = np.linspace(0, num_time_steps * time_step, num_time_steps + 1)
# --------------------------------------------

# ---------------- ARRAYS -------------------
temperature = np.full(num_nodes, initial_temperature, dtype=float)  # current temperatures
temperature_next = temperature.copy()  # next-step temperatures
fourier_number = thermal_diffusivity * time_step / node_spacing ** 2  # scalar, <= 0.5

# For plotting
plot_every = max(1, num_time_steps // 10)  # plot 10 curves
x = np.linspace(0, wall_thickness, num_nodes)

plt.figure(figsize=(9, 6))
plt.xlabel("x [m]")
plt.ylabel("Temperature [K]")
plt.title(f"1-D {num_nodes}-node wall - explicit FTCS")
plt.grid(True)

# ---------------- TIME LOOP ----------------
for step in range(num_time_steps):
    # Internal nodes (1 … N-2)
    for idx in range(1, num_nodes - 1):
        temperature_next[idx] = (
            temperature[idx]
            + fourier_number * (temperature[idx + 1] - 2 * temperature[idx] + temperature[idx - 1])
        )

    # Node 0: incident flux – radiation – conduction to node 1
    radiative_loss = emissivity * stefan_boltzmann * temperature[0] ** 4
    net_heat_flux_0 = (
        incident_heat_flux
        - radiative_loss
        - thermal_conductivity * (temperature[0] - temperature[1]) / node_spacing
    )
    control_volume_0 = node_spacing / 2.0
    temperature_next[0] = temperature[0] + time_step * net_heat_flux_0 / (
        density * specific_heat * control_volume_0
    )

    # Node N-1: conduction from node N-2 – cooling flux
    net_heat_flux_last = (
        thermal_conductivity * (temperature[num_nodes - 2] - temperature[num_nodes - 1]) / node_spacing
        - cooling_heat_flux
    )
    control_volume_last = node_spacing / 2.0
    temperature_next[num_nodes - 1] = temperature[num_nodes - 1] + time_step * net_heat_flux_last / (
        density * specific_heat * control_volume_last
    )

    # Advance in time
    temperature[:] = temperature_next[:]

    # Plot a handful of profiles
    if step % plot_every == 0 or step == num_time_steps - 1:
        plt.plot(x, temperature, label=f"t = {step * time_step:>5.0f} s")

# ---------------- RESULTS ------------------
plt.legend()
plt.tight_layout()
plt.show()

print("Simulation finished")
print(f"Δx = {node_spacing * 1e3:5.3f} mm,   dt = {time_step:7.4f} s  (Fo = {fourier_number:5.3f})")
print(f"Final surface temperature (node 0): {temperature[0]:.1f} K")
print(f"Final centre temperature (node (N-1)/2):  {temperature[int((num_nodes - 1) / 2)]:.1f} K")
print(f"Final cooled face temperature (node N-1): {temperature[num_nodes - 1]:.1f} K")
