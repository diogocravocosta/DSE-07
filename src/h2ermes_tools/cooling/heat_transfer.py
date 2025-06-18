"""
1-D transient heat conduction in a stainless-steel wall
- Variable incident heat flux at the hot face
- T⁴ radiation from the hot face
- Variable cooling flux at the cold face
Material: AISI 310 (properties assumed T-independent)
Scheme  : Explicit FTCS, uniform grid
"""

# plotting
import matplotlib.pyplot as plt
from matplotlib import cm
from matplotlib.colors import Normalize
import matplotlib
from matplotlib.animation import FuncAnimation

# calculations
import numpy as np
from data.constants import boltzmann
from scipy.interpolate import interp1d
from pyfluids import Fluid, FluidsList, Input

# h2ermes_tools imports
from h2ermes_tools.cooling.coolant import Coolant, hydrogen
from h2ermes_tools.cooling.material import Material, SS310
from h2ermes_tools.cooling.channel import RectangularChannel


class HeatShield:
    """
    Provides a simulation of the following scenario:

    1-D transient heat conduction in a stainless-steel wall
    - Variable incident heat flux at the hot face
    - T⁴ radiation from the hot face
    - Variable cooling flux at the cold face
    Material: AISI 310 (properties assumed T-independent)
    Scheme  : Explicit FTCS, uniform grid
    """

    def __init__(
        self,
        wall_thickness: float = 10e-3,
        num_nodes: int = 3,
        incident_heat_flux: np.ndarray = np.array([[0, 100_000]]),
        initial_temperature: float = 20.0,
        total_time: float = 900.0,
        coolant: Coolant = hydrogen,
        material: Material = SS310,
        heat_shield_diameter: float = 10.0,
        sphere_height: float = 2.0,
    ):
        self.wall_thickness = wall_thickness
        self.num_nodes = num_nodes
        self.incident_heat_flux = incident_heat_flux
        self.initial_temperature = initial_temperature
        self.total_time = total_time
        self.coolant = coolant
        self.material = material
        self.heat_shield_density = material.density
        self.specific_heat = material.specific_heat
        self.thermal_conductivity = material.thermal_conductivity
        self.stefan_boltzmann = boltzmann
        self.diameter = heat_shield_diameter
        self.height = sphere_height
        self.surface_area = np.pi * self.diameter * self.height

        # Validate parameters
        if self.num_nodes < 3:
            raise ValueError("num_nodes must be at least 3 for a valid simulation.")
        if self.wall_thickness <= 0:
            raise ValueError("wall_thickness must be a positive value.")
        if self.initial_temperature < 13.8:
            raise ValueError("initial_temperature must be above 13.8 K.")
        if self.total_time <= 0:
            raise ValueError("total_time must be a positive value.")
        if not isinstance(self.incident_heat_flux, np.ndarray):
            raise ValueError(
                "incident_heat_flux must be a numpy array with shape (N, 2) for time and value."
            )
        if self.incident_heat_flux.shape[1] != 2:
            raise ValueError(
                "incident_heat_flux must have two columns: time and heat flux value."
            )
        if not isinstance(self.coolant, Coolant):
            raise ValueError("coolant must be an instance of the Coolant class.")
        if not isinstance(self.material, Material):
            raise ValueError("material must be an instance of the Material class.")
        if self.diameter <= 0:
            raise ValueError("heat_shield_diameter must be a positive value.")
        if self.height <= 0:
            raise ValueError("sphere_height must be a positive value.")
        print("HeatShield class initialized. All parameters are valid.")

    def run_1d_simulation(self):
        """
        The following code will look similar to the run_simulation method,
        except it will calculate the area over which the incident heat flux is applied
        and the area over which the cooling flux is applied.

        TODO:
        - why I am copying all all the numbers into their variables?
        - make the initial temperature an array of temperature distribution
        """
        # 0. Initial Conditions, Constants & Setup
        T_i = self.initial_temperature
        temperature = np.full(self.num_nodes, T_i, dtype=float)
        temperature_next = temperature.copy()
        x = np.linspace(0, self.wall_thickness, self.num_nodes)

        node_spacing = self.wall_thickness / (self.num_nodes - 1)
        stefan_boltzmann = self.stefan_boltzmann

        # For animation: store all profiles
        all_temperatures = []
        all_times = []

        # 1. Heat shield material properties
        density = self.heat_shield_density
        specific_heat = self.specific_heat
        thermal_conductivity = self.thermal_conductivity

        # 2. Time step calculation (for stability)
        T_ref = 30  # Reference temperature for material properties
        # Below this reference temperature, the thermal diffusivity spikes
        # and a much lower time step is needed for stability
        cp_ref = specific_heat(T_ref)
        k_ref = thermal_conductivity(T_ref)
        thermal_diffusivity = k_ref / (density * cp_ref)
        max_dt = 0.5 * node_spacing**2 / thermal_diffusivity

        # If the coolant temperature drops below T_ref, manual adjustment to the time step
        # is needed to maintain numerical stability, 0.05 is generally sufficient
        # time_step = 0.05 * max_dt
        time_step = 0.9 * max_dt

        num_time_steps = int(np.ceil(self.total_time / time_step))
        fourier_number = thermal_diffusivity * time_step / node_spacing**2

        # 3. Incident heat flux input
        incident_heat_flux = self.incident_heat_flux
        times = incident_heat_flux[:, 0]  # expecting shape (N, 2): time, value
        values = incident_heat_flux[:, 1]
        incident_flux_func = interp1d(times, values)

        # 4. Start temporal simulation loop
        for step in range(num_time_steps):
            t_now = step * time_step
            # Evaluate temperature-dependent properties at current node temperatures
            cp_nodes = specific_heat(temperature)
            k_nodes = thermal_conductivity(temperature)

            # Internal nodes
            for idx in range(1, self.num_nodes - 1):
                # Use average k for conduction between nodes
                k_left = k_nodes[idx - 1]
                k_right = k_nodes[idx]
                k_avg = 0.5 * (k_left + k_right)

                # Update temperature using FTCS scheme
                temperature_next[idx] = temperature[idx] + (
                    time_step
                    * k_avg
                    / (density * cp_nodes[idx] * node_spacing**2)
                    * (
                        temperature[idx + 1]
                        - 2 * temperature[idx]
                        + temperature[idx - 1]
                    )
                )

            # Node 0: incident flux – radiation – conduction to node 1
            radiative_loss = (
                self.material.emissivity * stefan_boltzmann * temperature[0] ** 4
            )
            q_incident = incident_flux_func(
                t_now
            )  # TODO: this is what I want to adjust
            net_heat_flux_0 = (
                q_incident
                - radiative_loss
                - k_nodes[0] * (temperature[0] - temperature[1]) / node_spacing
            )
            control_volume_0 = node_spacing / 2.0
            temperature_next[0] = temperature[0] + time_step * net_heat_flux_0 / (
                density * cp_nodes[0] * control_volume_0
            )

            # Node N-1: conduction from node N-2 – cooling flux
            kN = k_nodes[-1]
            cpN = cp_nodes[-1]

            # Coolant cooling logic
            cold_wall_temp = temperature[self.num_nodes - 1]

            # Coolant bulk temperature
            T_bulk = self.coolant.fluid.temperature

            # Calculate heat transfer coefficient (update every time step!)
            h_cool = self.coolant.get_heat_transfer_coefficient()

            # Cooling heat flux (Fourier's law at interface)
            cooling_flux = h_cool * (cold_wall_temp - T_bulk)

            # Net heat flux at the last node
            net_heat_flux_last = (
                kN
                * (temperature[self.num_nodes - 2] - temperature[self.num_nodes - 1])
                / node_spacing
                - cooling_flux
            )
            control_volume_last = node_spacing / 2.0
            temperature_next[self.num_nodes - 1] = temperature[
                self.num_nodes - 1
            ] + time_step * net_heat_flux_last / (density * cpN * control_volume_last)
            temperature_next[self.num_nodes - 1] -= (
                cooling_flux * time_step / (density * cpN * control_volume_last)
            )

            # --- Coolant heating logic ---
            # segment_length = coolant.channel.length
            # area = coolant.channel.get_contact_area(segment_length)  # m^2
            # energy_to_coolant = cooling_flux * area * time_step  # Joules
            # coolant.add_energy(energy_to_coolant, time_step)

            # --- Pressure drop calculation: update coolant pressure ---
            # Use Darcy-Weisbach or other method in coolant class
            # This should update coolant.fluid.pressure and store history
            # coolant.update_pressure_drop(time_step)


            # Update temperature array
            temperature[:] = temperature_next[:]

            # Store for animation
            all_temperatures.append(temperature.copy())
            all_times.append(step * time_step)

            return (all_times, all_temperatures)

    def run_simulation(self):
        density = self.heat_shield_density
        specific_heat = self.specific_heat
        thermal_conductivity = self.thermal_conductivity

        # Use initial temperature for constants if property is callable
        T_ref = self.initial_temperature
        cp_ref = specific_heat(T_ref) if callable(specific_heat) else specific_heat
        k_ref = (
            thermal_conductivity(T_ref)
            if callable(thermal_conductivity)
            else thermal_conductivity
        )
        thermal_diffusivity = k_ref / (density * cp_ref)
        node_spacing = self.wall_thickness / (self.num_nodes - 1)
        max_dt = 0.5 * node_spacing**2 / thermal_diffusivity

        # manual adjustment for stability at lower temperatures
        # time_step needs to be smaller for higher thermal diffusivity
        # which occurs at lower temperatures
        time_step = 0.05 * max_dt

        num_time_steps = int(np.ceil(self.total_time / time_step))
        fourier_number = thermal_diffusivity * time_step / node_spacing**2
        stefan_boltzmann = (
            self.stefan_boltzmann
        )  # Stefan-Boltzmann constant [W/(m²·K⁴)]

        temperature = np.full(self.num_nodes, self.initial_temperature, dtype=float)
        temperature_next = temperature.copy()
        x = np.linspace(0, self.wall_thickness, self.num_nodes)

        # For animation: store all profiles
        all_temperatures = []
        all_times = []

        # Interpolate incident_heat_flux if it is an array (time, value)
        incident_heat_flux = self.incident_heat_flux

        # Expecting shape (N, 2): time, value
        times = incident_heat_flux[:, 0]
        values = incident_heat_flux[:, 1]

        # TODO: add try-except for interpolation errors etc. here
        incident_flux_func = interp1d(times, values)

        for step in range(num_time_steps):
            t_now = step * time_step
            # Evaluate temperature-dependent properties at current node temperatures
            cp_nodes = specific_heat(temperature)
            k_nodes = thermal_conductivity(temperature)

            # Internal nodes
            for idx in range(1, self.num_nodes - 1):
                # Use average k for conduction between nodes
                k_left = k_nodes[idx - 1]
                k_right = k_nodes[idx]
                k_avg = 0.5 * (k_left + k_right)

                # Update temperature using FTCS scheme
                temperature_next[idx] = temperature[idx] + (
                    time_step
                    * k_avg
                    / (density * cp_nodes[idx] * node_spacing**2)
                    * (
                        temperature[idx + 1]
                        - 2 * temperature[idx]
                        + temperature[idx - 1]
                    )
                )

            # Node 0: incident flux – radiation – conduction to node 1
            radiative_loss = (
                self.material.emissivity * stefan_boltzmann * temperature[0] ** 4
            )
            k0 = k_nodes[0]
            cp0 = cp_nodes[0]
            q_incident = incident_flux_func(t_now)
            net_heat_flux_0 = (
                q_incident
                - radiative_loss
                - k0 * (temperature[0] - temperature[1]) / node_spacing
            )
            control_volume_0 = node_spacing / 2.0
            temperature_next[0] = temperature[0] + time_step * net_heat_flux_0 / (
                density * cp0 * control_volume_0
            )

            # Node N-1: conduction from node N-2 – cooling flux
            kN = k_nodes[-1]
            cpN = cp_nodes[-1]

            # --- Coolant cooling logic ---
            # Determine coolant-side interface temperature (wall cold face)
            cold_wall_temp = temperature[self.num_nodes - 1]

            # Coolant bulk temperature
            T_bulk = self.coolant.fluid.temperature

            # Calculate heat transfer coefficient (update every time step!)
            h_cool = self.coolant.get_heat_transfer_coefficient()
            # Cooling heat flux (Fourier's law at interface)
            cooling_flux = h_cool * (cold_wall_temp - T_bulk)

            net_heat_flux_last = (
                kN
                * (temperature[self.num_nodes - 2] - temperature[self.num_nodes - 1])
                / node_spacing
                - cooling_flux
            )
            control_volume_last = node_spacing / 2.0
            temperature_next[self.num_nodes - 1] = temperature[
                self.num_nodes - 1
            ] + time_step * net_heat_flux_last / (density * cpN * control_volume_last)

            # --- Coolant heating logic ---
            # segment_length = coolant.channel.length
            # area = coolant.channel.get_contact_area(segment_length)  # m^2
            # energy_to_coolant = cooling_flux * area * time_step  # Joules
            # coolant.add_energy(energy_to_coolant, time_step)

            # --- Pressure drop calculation: update coolant pressure ---
            # Use Darcy-Weisbach or other method in coolant class
            # This should update coolant.fluid.pressure and store history
            # coolant.update_pressure_drop(time_step)

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

    def run_spacetime_simulation(self, channel_length=5.0, num_segments=4):
        """
        Marches in both time and space: for each time step, solves the 1-D wall heat transfer at each space segment,
        updating coolant properties and passing them to the next segment.
        """
        all_temperatures_segments = []
        all_times = None
        x = None
        fourier_number = None
        time_step = None
        node_spacing = None
        coolant_states = []
        # Initial coolant state for the first segment
        coolant = self.coolant
        for seg in range(num_segments):
            # Create a new HeatShield for this segment, with the current coolant state
            segment_sim = HeatShield(
                wall_thickness=self.wall_thickness,
                num_nodes=self.num_nodes,
                incident_heat_flux=self.incident_heat_flux,
                initial_temperature=self.initial_temperature
                if seg == 0
                else all_temperatures_segments[-1][0, -1],
                total_time=self.total_time,
                coolant=coolant,
                material=self.material,
            )
            # Run the 1-D simulation for this segment
            all_temperatures, all_times, x, fourier_number, time_step, node_spacing = (
                segment_sim.run_simulation()
            )
            all_temperatures_segments.append(all_temperatures)
            # After each time step, update coolant properties for the next segment
            # (Mockup: only update at the end of the segment for now)
            # Use the appropriate methods from coolant.py:
            #   - coolant.add_energy(...)
            #   - coolant.update_pressure_drop(...)
            # For a real implementation, this should be done at each time step and passed to the next segment
            # Here, we simply copy the coolant object (in reality, you would update its state)
            coolant_states.append(coolant)
            # For the next segment, pass the updated coolant (deepcopy if needed)
            # coolant = deepcopy(coolant)  # If needed, import deepcopy
            # (For now, just pass the same object for mockup)
        return (
            all_temperatures_segments,
            all_times,
            x,
            fourier_number,
            time_step,
            node_spacing,
            coolant_states,
        )

    @staticmethod
    def plot_profiles(
        all_temperatures,
        all_times,
        x,
        num_profiles=10,
        save=False,
    ):
        min_temp = 13.8
        max_temp = getattr(SS310, "maximum_temperature", 1500)
        norm = Normalize(vmin=min_temp, vmax=max_temp)
        cmap = matplotlib.colormaps["plasma"]
        sm = cm.ScalarMappable(norm=norm, cmap=cmap)
        plt.figure(figsize=(9, 6))
        plt.xlabel("x [m]")
        plt.ylabel("Temperature [K]")
        # plt.title(f"1-D {len(x)}-node wall - explicit FTCS")
        plt.grid(True)
        # Select evenly spaced indices to plot exactly num_profiles profiles (including first and last)
        if num_profiles > 1:
            indices = np.linspace(0, len(all_times) - 1, num_profiles, dtype=int)
        else:
            indices = [len(all_times) - 1]
        for i in indices:
            color = cmap(norm(np.max(all_temperatures[i])))
            plt.plot(
                x, all_temperatures[i], label=f"t = {all_times[i]:>5.0f} s", color=color
            )
        # Always plot the last profile, not used at the moment
        # color = cmap(norm(np.max(all_temperatures[-1])))
        # plt.plot(x, all_temperatures[-1], label=f"t = {all_times[-1]:>5.0f} s", color=color)

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
        plt.savefig("heat_shield_sim.pdf")
        plt.show()

    @staticmethod
    def animate_profiles(
        all_temperatures,
        all_times,
        x,
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
            fig,
            update,
            frames=len(all_times),
            init_func=init,
            blit=True,
            interval=interval,
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

    @staticmethod
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

    def estimate_heat_shield_spherical_mass(self) -> float:
        """
        Estimate the mass of a spherical cap heat shield as a function of diameter.
        Args:
            h (float): Height of the spherical cap [m] (vertical distance from base to top of cap)
            diameter (float, optional): Inner diameter of the sphere [m]. If None, uses self.diameter.
        Returns:
            float: Mass [kg]
        """
        mass = self.surface_area * self.wall_thickness * self.heat_shield_density
        return mass

    def estimate_heat_shield_mass(self) -> float:
        """
        Estimate the heat shield mass. This mass is composed of two spherical masses.
        The first mass is the solid, full heat shield wall and the second mass is the
        portion of the wall that has coolant channels inside it.

        The second part is essentially hollowed out by the cooling channels. The coolant
        channels are assumed to be layed through 50% of the surface.

        Note that the 1e-3 is added as an estimated wall thickness of the coolant channel
        on the cold side of the coolant.
        """
        coolant_cold_side_wall_thickness = 1e-3
        first_mass = self.estimate_heat_shield_spherical_mass()
        radius = (
            self.diameter / 2
            - self.coolant.channel.height
            - coolant_cold_side_wall_thickness
        )
        area = 2 * np.pi * radius * self.height * 0.5  # 50% of the surface
        second_mass = (
            area
            * (self.coolant.channel.height + coolant_cold_side_wall_thickness)
            * self.heat_shield_density
        )
        return first_mass + second_mass

    def estimate_coolant_loop_length(self, fraction) -> float:
        """
        This estimation assumes the coolant channels cover a certain portion of the
        surface area of the heat shield. Their length is estimated as if there was
        one long coolant channel whose total contact area is equal to the portion of the
        surface area that is covered by the coolant channels.

        Returns:
            float: Estimated length of the coolant channel [m]
        """
        contact_area = self.surface_area * fraction
        return contact_area / self.coolant.channel.width

    def calculate_coolant_loop_volume(self, fraction) -> float:
        """
        Calculate the volume of the coolant loop based on the coolant channel dimensions.

        Returns:
            float: Volume of the coolant loop [m^3]
        """
        # Assuming a rectangular channel for simplicity
        return (
            self.coolant.channel.cross_sectional_area
            * self.estimate_coolant_loop_length(fraction)
        )

    def estimate_coolant_mass_of_filled_loop(self, fraction) -> float:
        """
        Estimate the mass of the coolant in the loop based on the coolant channel volume
        and the fluid density.

        Returns:
            float: Mass of the coolant in the loop [kg]
        """
        coolant_volume = self.calculate_coolant_loop_volume(fraction)  # Full loop
        return coolant_volume * self.coolant.fluid.density


if __name__ == "__main__":
    # Example usage
    hydrogen = Fluid(FluidsList.Hydrogen).with_state(
        Input.pressure(50e5),  # Pa
        Input.temperature(13.8),  # K
    )
    channel = RectangularChannel(width=10e-3, height=2e-3, length=5.0, roughness=1e-5)
    coolant = Coolant(fluid=hydrogen, channel=channel, mass_flow=0.001)

    hs = HeatShield(
        wall_thickness=3e-3,
        num_nodes=3,
        incident_heat_flux=np.array(
            [[0, 100_000], [150, 666_000], [450, 666_000], [900, 0]]
        ),
        initial_temperature=30.0,
        total_time=900.0,
        coolant=coolant,
        material=SS310,
        heat_shield_diameter=10.0,
        sphere_height=2.5,
    )
    all_temperatures, all_times, x, fourier_number, time_step, node_spacing = hs.run_simulation()
    hs.plot_profiles(all_temperatures, all_times, x)
    # hs.animate_profiles(all_temperatures, all_times, x)
    # hs.print_summary(all_temperatures[-1], node_spacing, time_step, fourier_number)

    mass = hs.estimate_heat_shield_mass()
    print(f"Estimated Heat Shield Mass: {mass:.2f} kg")

    # Estimate coolant channel length for 50% surface area coverage
    coolant_channel_length = hs.estimate_coolant_loop_length(0.5)
    print(
        f"Estimated Coolant Channel Length (50% coverage): {coolant_channel_length:.2f} m"
    )

    # Estimate coolant mass in the loop
    coolant_mass = hs.estimate_coolant_mass_of_filled_loop(0.5)
    print(f"Estimated Coolant Mass in Loop (50% coverage): {coolant_mass:.2f} kg")

    # Run a 1-D simulation
    all_times, all_temperatures = hs.run_1d_simulation()
