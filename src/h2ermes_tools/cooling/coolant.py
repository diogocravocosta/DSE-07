import matplotlib.pyplot as plt
import numpy as np
from pyfluids import Fluid, FluidsList, Input
from h2ermes_tools.variables import coolant_inlet_pressure, coolant_inlet_temperature
from h2ermes_tools.cooling.channel import CircularChannel, RectangularChannel


class Coolant:
    """
    Coolant class for managing coolant properties and calculations.

    Attributes:
        fluid (Fluid): The Fluid object from pyfluids representing the coolant.
        channel (Channel): The channel through which the coolant flows.
        mass_flow (float): The mass flow rate of the coolant in kg/s.
    """

    def __init__(self, fluid, channel, mass_flow):
        self.fluid = fluid
        self.channel = channel
        self.mass_flow = mass_flow
        self.reynolds_number = self.get_reynolds_number()
        self.nusselt_number = self.get_nusselt_number_taylor()
        self.fluid_speed = self.get_fluid_speed()
        self.heat_transfer_coefficient = self.get_heat_transfer_coefficient()

    def plot_coolant_properties(
        self, temperatures=np.linspace(13.8, 300, 100), pressure=10e5
    ):
        enthalpies = [
            self.fluid.with_state(
                Input.temperature(T), Input.pressure(pressure)
            ).enthalpy
            for T in temperatures
        ]
        density = [
            self.fluid.with_state(
                Input.temperature(T), Input.pressure(pressure)
            ).density
            for T in temperatures
        ]

        specific_heats = [
            self.fluid.with_state(
                Input.temperature(T), Input.pressure(pressure)
            ).specific_heat
            for T in temperatures
        ]
        dynamic_viscosities = [
            self.fluid.with_state(
                Input.temperature(T), Input.pressure(pressure)
            ).dynamic_viscosity
            for T in temperatures
        ]

        fig, axs = plt.subplots(4, 1, figsize=(15, 7))
        axs[0].plot(temperatures, enthalpies)
        axs[0].set_ylabel(r"Enthalpy (J/kg)")
        axs[0].grid()
        # No legend

        axs[1].plot(temperatures, density)
        axs[1].set_ylabel(r"Density (kg/m$^3$)")
        axs[1].grid()
        # No legend

        axs[2].plot(temperatures, specific_heats)
        axs[2].set_ylabel(r"Specific Heat $c_p$ (J/kg$\cdot$K)")
        axs[2].grid()
        # No legend

        axs[3].plot(temperatures, dynamic_viscosities)
        axs[3].set_xlabel(r"Temperature $T$ (K)")
        axs[3].set_ylabel(r"Dynamic Viscosity $\mu$ (Pa$\cdot$s)")
        axs[3].grid()
        # No legend

        plt.tight_layout()
        plt.show()

    def get_nusselt_number_taylor(
        self, temperature_ratio=0.55, dimensionless_length=1.0
    ) -> float:
        """
        Calculate the Nusselt number based on the Taylor relationship.

        Parameters:
        - reynolds_number: Reynolds number (dimensionless)
        - temperature_ratio: Ratio of surface temperature to bulk temperature (dimensionless)
        - dimensionless_length: Dimensionless length (x/D, where x is the length and D is the diameter)

        Returns:
        - Nusselt number calculated using the Taylor empirical relationship
        """
        return (
            0.023
            * (self.reynolds_number**0.8)
            * (self.fluid.prandtl**0.4)
            * (temperature_ratio ** (-0.57 - (1.59 / dimensionless_length)))
        )

    def get_reynolds_number(self) -> float:
        """
        Calculate the Reynolds number.

        Returns:
        - Reynolds number (dimensionless)
        """
        fluid_speed = self.get_fluid_speed()
        return (
            self.fluid.density * fluid_speed * self.channel.hydraulic_diameter
        ) / self.fluid.dynamic_viscosity

    def get_fluid_speed(self):
        """
        Calculate the fluid speed based on mass flow rate.

        Returns:
        - Fluid speed (m/s)
        """
        return self.mass_flow / (self.fluid.density * self.channel.cross_sectional_area)

    def get_heat_transfer_coefficient(self):
        """
        Calculate the heat transfer coefficient.

        Returns:
        - Heat transfer coefficient (W/m²·K)
        """
        return (
            self.nusselt_number
            * self.fluid.conductivity
            / self.channel.hydraulic_diameter
        )

    def add_energy(self, energy):
        """
        Add energy to the coolant and update its temperature.
        Parameters:
        - energy: Energy to be added (J)
        - dt: Time step (s), optional (for history tracking)
        """
        # Calculate temperature rise: dT = Q / (m * cp)
        cp = self.fluid.specific_heat  # J/kg/K
        dT = energy / (self.mass_flow * cp)
        new_temp = self.fluid.temperature + dT
        # Update fluid state (keep pressure constant for now)
        self.fluid = self.fluid.with_state(
            Input.pressure(self.fluid.pressure),
            Input.temperature(new_temp),
        )

    def calculate_pressure_drop(self, segment_length=None):
        """
        Calculate and update the coolant pressure using Darcy-Weisbach equation
        on a section of the channel of a limited length.
        """
        # Darcy-Weisbach: dP = f * (L/D) * (rho*v^2/2)
        L = segment_length
        D = self.channel.hydraulic_diameter
        v = self.get_fluid_speed()
        rho = self.fluid.density
        # Estimate friction factor (f) using Blasius for turbulent, else laminar
        Re = self.get_reynolds_number()
        if Re < 2300:
            f = 64.0 / Re
        else:
            f = 0.3164 * Re**-0.25
        dP = f * (L / D) * (rho * v**2 / 2)

        return dP

    def reduce_pressure(self, pressure_drop):
        """
        Reduce the coolant pressure by a specified amount.
        Parameters:
        - pressure_drop: Pressure drop to be applied (Pa)
        """
        new_pressure = self.fluid.pressure - pressure_drop
        # if new_pressure < 0:
        #     raise ValueError("Pressure cannot be negative.")
        self.fluid = self.fluid.with_state(
            Input.pressure(new_pressure), Input.temperature(self.fluid.temperature)
        )


hydrogen = Fluid(FluidsList.Hydrogen).with_state(
    Input.pressure(coolant_inlet_pressure.value),  # Pa
    Input.temperature(coolant_inlet_temperature.value),  # K
)  # Initialize hydrogen coolant with inlet conditions, useful for testing

if __name__ == "__main__":
    channel = RectangularChannel(width=10e-3, height=5e-3, length=5.0, roughness=1e-5)
    channel = CircularChannel(diameter=10e-3, length=5.0, roughness=1e-5)
    coolant = Coolant(fluid=hydrogen, channel=channel, mass_flow=0.1)

    # --- Quick Tests ---
    coolant.plot_coolant_properties()
    reynolds_number = coolant.get_reynolds_number()
    nusselt_number = coolant.get_nusselt_number_taylor(reynolds_number)
    heat_transfer_coefficient = coolant.get_heat_transfer_coefficient()
    fluid_speed = coolant.get_fluid_speed()
    print(f"Reynolds Number: {reynolds_number}")
    print(f"Nusselt Number: {nusselt_number}")
    print(f"Heat Transfer Coefficient: {heat_transfer_coefficient} W/m²·K")
    print(f"Fluid Speed: {fluid_speed} m/s")

    # Add energy to the coolant
    energy_added = 670_000.0  # J
    print(f"Old Coolant Temperature: {coolant.fluid.temperature} K")
    coolant.add_energy(energy_added)
    print(f"New Coolant Temperature: {coolant.fluid.temperature} K")

    # Pressure drop calculation
    section_length = 1.0  # m
    pressure_drop = coolant.calculate_pressure_drop(section_length)
    print(f"Pressure Drop over {section_length} m: {pressure_drop} Pa")

    # Reduce coolant pressure
    print(f"Old Coolant Pressure: {coolant.fluid.pressure * 1e-5} bar")
    coolant.reduce_pressure(pressure_drop)
    print(f"New Coolant Pressure: {coolant.fluid.pressure * 1e-5} bar")
