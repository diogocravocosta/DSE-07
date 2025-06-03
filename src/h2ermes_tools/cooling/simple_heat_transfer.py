"""
1-D 10-node stainless-steel wall with
 • constant incident flux at the hot face
 • T⁴ radiation from the hot face
 • constant cooling flux at the cold face
Material: AISI 304L (properties assumed T-independent)
Scheme  : Explicit FTCS, uniform grid
"""

import numpy as np
import matplotlib.pyplot as plt

# ----------  USER INPUT  ----------------------------------------------------
# Geometry & discretisation
L = 0.010  # [m] total wall thickness
N = 31  # number of nodes
dx = L / (N - 1)  # [m] node spacing

# Material (304L, room-temperature values)
rho = 8_030.0  # [kg m-3]
cp = 500.0  # [J kg-1 K-1]
k = 14.6  # [W m-1 K-1]
alpha = k / (rho * cp)  # [m² s-1] thermal diffusivity

# Boundary heat fluxes ( + = heat INTO the wall, – = heat OUT )
q_inc = 100e3  # [W m-2] constant aerodynamic heating on node 0
q_c = 30e3  # [W m-2] constant cooling on node 9  (drawn out)

# Radiation
eps = 0.78  # emissivity, typical oxidised 304L
sigma = 5.670374419e-8  # [W m-2 K-4] Stefan–Boltzmann constant

# Initial & run parameters
T_init = 20.0  # [K] uniform start temperature
t_end = 900.0  # [s] total simulated time

# Stability-limited time-step (Fo = α dt / dx² ≤ 0.5)
dt_max = 0.5 * dx**2 / alpha
dt = 0.9 * dt_max  # 90 % of the limit for safety
Nt = int(np.ceil(t_end / dt))  # number of steps (integer)
time = np.linspace(0, Nt * dt, Nt + 1)  # time axis for plotting
# ---------------------------------------------------------------------------


# ----------  ARRAYS ---------------------------------------------------------
T = np.full(N, T_init, dtype=float)  # current temperatures
Tn1 = T.copy()  # next-step temperatures
Fo = alpha * dt / dx**2  # Fourier number (scalar, <= 0.5)

# For plotting
plot_every = max(1, Nt // 10)  # plot 10 curves
x = np.linspace(0, L, N)

plt.figure(figsize=(9, 6))
plt.xlabel("x [m]")
plt.ylabel("Temperature [K]")
plt.title(f"1-D {N}-node wall - explicit FTCS")
plt.grid(True)

# ----------  TIME LOOP ------------------------------------------------------
for n in range(Nt):
    # --- internal nodes (1 … 8) ---
    for i in range(1, N - 1):
        Tn1[i] = T[i] + Fo * (T[i + 1] - 2 * T[i] + T[i - 1])

    # --- node 0 : incident flux – radiation – conduction to node 1 -----------
    q_r = eps * sigma * T[0] ** 4  # Stefan-Boltzmann
    net_q_0 = q_inc - q_r - k * (T[0] - T[1]) / dx  # W m-2, into CV
    vol_0 = dx / 2.0  # [m] half-CV thickness
    Tn1[0] = T[0] + dt * net_q_0 / (rho * cp * vol_0)

    # --- node N-1 : conduction from node N-2 – cooling flux ----------------------
    net_q_last = k * (T[N - 2] - T[N - 1]) / dx - q_c  # +ve into CV
    vol_last = dx / 2.0
    Tn1[N - 1] = T[N - 1] + dt * net_q_last / (rho * cp * vol_last)

    # advance in time
    T[:] = Tn1[:]

    # plot a handful of profiles
    if n % plot_every == 0 or n == Nt - 1:
        plt.plot(x, T, label=f"t = {n * dt:>5.0f} s")

# ----------  RESULTS --------------------------------------------------------
plt.legend()
plt.tight_layout()
plt.show()

print("Simulation finished")
print(f"Δx = {dx * 1e3:5.3f} mm,   dt = {dt:7.4f} s  (Fo = {Fo:5.3f})")
print(f"Final surface temperature (node 0): {T[0]:.1f} K")
print(f"Final centre temperature (node (N-1)/2):  {T[int((N - 1) / 2)]:.1f} K")
print(f"Final cooled face temperature (node N-1): {T[N - 1]:.1f} K")
