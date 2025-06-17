import numpy as np
import matplotlib.pyplot as plt
from pyfluids import Fluid, FluidsList, Input


def simulate_channel(mdot, p_in, T_in, L, D, q_flux, n_segments):
    fluid = Fluid(FluidsList.Hydrogen)
    dz = L / n_segments
    A_c = np.pi * (D / 2) ** 2
    A_h = np.pi * D * dz
    G = mdot / A_c

    z = np.linspace(0, L, n_segments + 1)
    p = np.zeros_like(z)
    T = np.zeros_like(z)
    rho = np.zeros_like(z)

    p[0] = p_in
    T[0] = T_in

    for i in range(n_segments):
        try:
            state_i = fluid.with_state(Input.pressure(p[i]), Input.temperature(T[i]))
        except Exception as e:
            print(
                f"Error at segment {i}, P={p[i] / 1e6:.2f} MPa, T={T[i]:.2f} K. Error: {e}"
            )
            break

        T[i + 1] = T[i] + (q_flux * A_h) / (mdot * state_i.specific_heat)

        state_i_plus_1 = fluid.with_state(
            Input.pressure(p[i]), Input.temperature(T[i + 1])
        )

        T_avg = (T[i] + T[i + 1]) / 2
        state_avg = fluid.with_state(Input.pressure(p[i]), Input.temperature(T_avg))

        Re = (G * D) / state_avg.dynamic_viscosity

        f_D = (0.79 * np.log(Re) - 1.64) ** -2 if Re >= 2300 else 64 / Re

        dp_friction = (f_D * G**2) / (2 * state_avg.density * D) * dz
        dp_momentum = G**2 * (1 / state_i_plus_1.density - 1 / state_i.density)

        p[i + 1] = p[i] - (dp_friction + dp_momentum)
        rho[i] = state_i.density

    rho[-1] = fluid.with_state(Input.pressure(p[-1]), Input.temperature(T[-1])).density

    total_dp = p[0] - p[-1]
    print(f"Inlet:  P = {p[0] / 1e6:.2f} MPa, T = {T[0]:.2f} K")
    print(f"Outlet: P = {p[-1] / 1e6:.2f} MPa, T = {T[-1]:.2f} K")
    print(f"Total Pressure Drop: {total_dp / 1e6:.4f} MPa")

    return {"z": z, "p": p, "T": T, "rho": rho}


def plot_results(results, inlet_pressure, channel=None, mass_flow=None):
    z, p, T, rho = results["z"], results["p"], results["T"], results["rho"]

    fig, axs = plt.subplots(2, 2, figsize=(12, 9), sharex=False)
    fig.suptitle(r"1D Supercritical Hydrogen Channel Simulation")
    axs = axs.flatten()

    axs[0].plot(z, p / 1e5, label=r"$P$")
    axs[0].set_ylabel(r"Pressure $P$ (bar)")
    axs[0].set_title(r"Pressure vs. Channel Length $z$")
    axs[0].legend()
    axs[0].grid(True)
    # No x-label for top row

    axs[1].plot(z, T, label=r"$T$")
    try:
        T_pc = FluidsList.ParaHydrogen.with_state(
            Input.pressure(inlet_pressure), Input.quality(0)
        ).temperature
        axs[1].axhline(
            y=T_pc,
            linestyle="--",
            label=fr"Pseudo-critical $T$ at {inlet_pressure / 1e6} MPa",
        )
        axs[1].legend()
    except Exception:
        pass  # Fails if pressure is above critical
    axs[1].set_ylabel(r"Temperature $T$ (K)")
    axs[1].set_title(r"Temperature vs. Channel Length $z$")
    axs[1].grid(True)
    # No x-label for top row

    axs[2].plot(z, rho, label=r"$\rho$")
    axs[2].set_ylabel(r"Density $\rho$ (kg/m$^3$)")
    axs[2].set_xlabel(r"Channel Length $z$ (m)")
    axs[2].set_title(r"Density vs. Channel Length $z$")
    axs[2].legend()
    axs[2].grid(True)

    # Nusselt number vs. channel length (Taylor relationship)
    if channel is not None and mass_flow is not None:
        nusselts = []
        for Ti in T:
            fluid = Fluid(FluidsList.ParaHydrogen).with_state(
                Input.temperature(Ti), Input.pressure(inlet_pressure)
            )
            rho_i = fluid.density
            mu = fluid.dynamic_viscosity
            D = channel.hydraulic_diameter
            A = channel.cross_sectional_area
            v = mass_flow / (rho_i * A)
            Re = rho_i * v * D / mu
            Pr = fluid.prandtl
            Nu = 0.023 * (Re ** 0.8) * (Pr ** 0.4) * (0.55 ** (-0.57 - (1.59 / 1.0)))
            nusselts.append(Nu)
        axs[3].plot(z, nusselts, label=r"Nusselt Number $Nu$")
        axs[3].set_ylabel(r"Nusselt Number $Nu$")
        axs[3].set_xlabel(r"Channel Length $z$ (m)")
        axs[3].set_title(r"Nusselt Number vs. Channel Length $z$")
        axs[3].legend()
        axs[3].grid(True)
    else:
        axs[3].set_visible(False)

    plt.tight_layout(rect=[0, 0, 1, 0.96])
    plt.show()


if __name__ == "__main__":
    p_inlet = 4e6  # Pascals
    mass_flow = 0.002  # kg/s
    channel_diameter = 0.01  # m
    channel_length = 20.0  # m
    channel_roughness = 1e-5  # m (example value)

    from h2ermes_tools.cooling.channel import CircularChannel

    channel = CircularChannel(channel_diameter, channel_length, channel_roughness)

    sim_results = simulate_channel(
        mdot=mass_flow,  # kg/s
        p_in=p_inlet,  # Pa
        T_in=20,  # K
        L=channel_length,  # m
        D=channel_diameter,  # m
        q_flux=50000,  # W/m^2
        n_segments=200,
    )

    if sim_results:
        plot_results(sim_results, p_inlet, channel=channel, mass_flow=mass_flow)
