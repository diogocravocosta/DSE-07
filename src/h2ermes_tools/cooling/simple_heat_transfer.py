"""
1-D transient heat conduction in a stainless-steel wall
• Constant incident heat flux at the hot face
• T⁴ radiation from the hot face
• Constant cooling flux at the cold face
Material: AISI 310 (properties assumed T-independent)
Scheme  : Explicit FTCS, uniform grid
"""

import numpy as np
import matplotlib.pyplot as plt
from matplotlib import cm
from matplotlib.colors import Normalize
import matplotlib
from matplotlib.animation import FuncAnimation

from h2ermes_tools.cooling.material import SS310


def run_simulation(
    wall_thickness=0.010,
    num_nodes=3,
    incident_heat_flux=100_000.0,
    cooling_heat_flux=30_000.0,
    emissivity=0.9,
    initial_temperature=20.0,
    total_time=900.0,
    thermal_conductivity=25,
    density=None,
    specific_heat=None,
    plot_every=10,
):
    # Material properties
    if density is None:
        density = SS310.density
    if specific_heat is None:
        specific_heat = SS310.specific_heat
    thermal_diffusivity = thermal_conductivity / (density * specific_heat)
    node_spacing = wall_thickness / (num_nodes - 1)
    max_dt = 0.5 * node_spacing**2 / thermal_diffusivity
    time_step = 0.9 * max_dt
    num_time_steps = int(np.ceil(total_time / time_step))
    fourier_number = thermal_diffusivity * time_step / node_spacing**2
    stefan_boltzmann = 5.670374419e-8

    temperature = np.full(num_nodes, initial_temperature, dtype=float)
    temperature_next = temperature.copy()
    x = np.linspace(0, wall_thickness, num_nodes)

    # For animation: store all profiles
    all_temperatures = []
    all_times = []

    for step in range(num_time_steps):
        # Internal nodes
        for idx in range(1, num_nodes - 1):
            temperature_next[idx] = temperature[idx] + fourier_number * (
                temperature[idx + 1] - 2 * temperature[idx] + temperature[idx - 1]
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
            thermal_conductivity
            * (temperature[num_nodes - 2] - temperature[num_nodes - 1])
            / node_spacing
            - cooling_heat_flux
        )
        control_volume_last = node_spacing / 2.0
        temperature_next[num_nodes - 1] = temperature[
            num_nodes - 1
        ] + time_step * net_heat_flux_last / (
            density * specific_heat * control_volume_last
        )
        temperature[:] = temperature_next[:]
        # Store for animation
        all_temperatures.append(temperature.copy())
        all_times.append(step * time_step)
    return (
        np.array(all_temperatures),
        np.array(all_times),
        x,
        fourier_number,
        time_step,
        node_spacing,
    )


def plot_profiles(
    all_temperatures,
    all_times,
    x,
    fourier_number,
    time_step,
    node_spacing,
    num_profiles=10,
):
    min_temp = 13.8
    max_temp = getattr(SS310, "maximum_temperature", 1500)
    norm = Normalize(vmin=min_temp, vmax=max_temp)
    cmap = matplotlib.colormaps["plasma"]
    sm = cm.ScalarMappable(norm=norm, cmap=cmap)
    plt.figure(figsize=(9, 6))
    plt.xlabel("x [m]")
    plt.ylabel("Temperature [K]")
    plt.title(f"1-D {len(x)}-node wall - explicit FTCS")
    plt.grid(True)
    plot_every = max(1, len(all_times) // num_profiles)
    for i in range(0, len(all_times), plot_every):
        color = cmap(norm(np.max(all_temperatures[i])))
        plt.plot(
            x, all_temperatures[i], label=f"t = {all_times[i]:>5.0f} s", color=color
        )
    # Always plot the last profile
    color = cmap(norm(np.max(all_temperatures[-1])))
    plt.plot(x, all_temperatures[-1], label=f"t = {all_times[-1]:>5.0f} s", color=color)
    # Add a red dashed line for the material maximum temperature
    plt.axhline(
        max_temp,
        color="red",
        linestyle="--",
        linewidth=2,
        label=f"{SS310.name} Maximum Peak Temperature: ({max_temp} K)",
    )
    _cbar = plt.colorbar(sm, label="Temperature [K]", ax=plt.gca())
    plt.legend()
    plt.tight_layout()
    plt.show()


def animate_profiles(
    all_temperatures,
    all_times,
    x,
    fourier_number,
    time_step,
    node_spacing,
    interval=50,
    save_path=None,
    save_format=None,
):
    min_temp = 13.8
    max_temp = getattr(SS310, "maximum_temperature", 1500)
    norm = Normalize(vmin=min_temp, vmax=max_temp)
    cmap = matplotlib.colormaps["plasma"]
    sm = cm.ScalarMappable(norm=norm, cmap=cmap)
    fig, ax = plt.subplots(figsize=(9, 6))
    ax.set_xlabel("x [m]")
    ax.set_ylabel("Temperature [K]")
    ax.set_title(f"1-D {len(x)}-node wall - explicit FTCS (Animated)")
    ax.grid(True)
    # Set axis limits for visibility
    ax.set_xlim(np.min(x), np.max(x))
    ax.set_ylim(min_temp - 10, max_temp + 110)
    # Initialize line with x and initial temperature
    (line,) = ax.plot(x, all_temperatures[0], lw=2)
    # Add a red dashed line for the material maximum temperature
    ax.axhline(
        max_temp,
        color="red",
        linestyle="--",
        linewidth=2,
        label=f"{SS310.name} Maximum Peak Temperature: ({max_temp} K)",
    )
    _cbar = plt.colorbar(sm, label="Temperature [K]", ax=ax)
    time_text = ax.text(0.02, 0.95, "", transform=ax.transAxes)

    def init():
        line.set_data(x, all_temperatures[0])
        time_text.set_text(f"t = {all_times[0]:.1f} s")
        return line, time_text

    def update(frame):
        y = all_temperatures[frame]
        color = cmap(norm(np.max(y)))
        line.set_data(x, y)
        line.set_color(color)
        time_text.set_text(f"t = {all_times[frame]:.1f} s")
        return line, time_text

    _ani = FuncAnimation(
        fig, update, frames=len(all_times), init_func=init, blit=True, interval=interval
    )
    ax.legend()
    plt.tight_layout()
    if save_path:
        # Determine writer based on format
        if save_format is None:
            # Infer from file extension
            if save_path.lower().endswith(".mp4"):
                save_format = "mp4"
            elif save_path.lower().endswith(".gif"):
                save_format = "gif"
            else:
                save_format = "mp4"  # default
        if save_format == "mp4":
            _ani.save(save_path, writer="ffmpeg", dpi=150)
        elif save_format == "gif":
            _ani.save(save_path, writer="pillow", dpi=100)
        else:
            raise ValueError(f"Unsupported save_format: {save_format}")
    plt.show()


def print_summary(temperature, node_spacing, time_step, fourier_number):
    num_nodes = len(temperature)
    print("Simulation finished")
    print(
        f"Δx = {node_spacing * 1e3:5.3f} mm,   dt = {time_step:7.4f} s  (Fo = {fourier_number:5.3f})"
    )
    print(f"Final surface temperature (node 0): {temperature[0]:.1f} K")
    print(
        f"Final centre temperature (node (N-1)/2):  {temperature[int((num_nodes - 1) / 2)]:.1f} K"
    )
    print(
        f"Final cooled face temperature (node N-1): {temperature[num_nodes - 1]:.1f} K"
    )


if __name__ == "__main__":
    # ---------------- USER INPUT ----------------
    do_animation = False  # Set to True for animation, False for static plot
    wall_thickness = 0.010
    num_nodes = 3
    incident_heat_flux = 200_000.0
    cooling_heat_flux = 30_000.0
    emissivity = 0.90
    initial_temperature = 20.0
    total_time = 900.0
    # --------------------------------------------
    all_temperatures, all_times, x, fourier_number, time_step, node_spacing = (
        run_simulation(
            wall_thickness=wall_thickness,
            num_nodes=num_nodes,
            incident_heat_flux=incident_heat_flux,
            cooling_heat_flux=cooling_heat_flux,
            emissivity=emissivity,
            initial_temperature=initial_temperature,
            total_time=total_time,
        )
    )
    if do_animation:
        animate_profiles(
            all_temperatures,
            all_times,
            x,
            fourier_number,
            time_step,
            node_spacing,
            # save_path="./heat_transfer_animation.gif",
            # save_format="gif",
        )
    else:
        plot_profiles(
            all_temperatures, all_times, x, fourier_number, time_step, node_spacing
        )
    print_summary(all_temperatures[-1], node_spacing, time_step, fourier_number)
