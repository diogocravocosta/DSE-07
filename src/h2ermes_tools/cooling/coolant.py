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
        conductivities = [
            self.fluid.with_state(
                Input.temperature(T), Input.pressure(pressure)
            ).conductivity
            for T in temperatures
        ]

        fig, axs = plt.subplots(4, 1, figsize=(15, 7))
        axs[0].plot(temperatures, enthalpies, label="Enthalpy")
        # axs[0].set_title("Coolant Enthalpy vs Temperature")
        # axs[0].set_xlabel("Temperature (K)")
        axs[0].set_ylabel("Enthalpy (J/kg)")
        axs[0].grid()
        axs[0].legend()

        axs[1].plot(temperatures, density, label="Density")
        # axs[1].set_title("Coolant Density vs Temperature")
        # axs[1].set_xlabel("Temperature (K)")
        axs[1].set_ylabel("Density (kg/m³)")
        axs[1].grid()
        axs[1].legend()

        axs[2].plot(temperatures, specific_heats, label="Specific Heat")
        # axs[2].set_title("Coolant Specific Heat vs Temperature")
        # axs[2].set_xlabel("Temperature (K)")
        axs[2].set_ylabel("Specific Heat (J/kg·K)")
        axs[2].grid()
        axs[2].legend()

        axs[3].plot(temperatures, conductivities, label="Thermal Conductivity")
        axs[3].set_xlabel("Temperature (K)")
        axs[3].set_ylabel("Thermal Conductivity (W/m·K)")
        axs[3].grid()
        axs[3].legend()

        plt.suptitle(
            f"Coolant Properties for Hydrogen at P = {round(coolant_inlet_pressure.value * 1e-5, 1)} bar"
        )

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
        return self.nusselt_number * self.fluid.conductivity / self.channel.diameter


if __name__ == "__main__":
    hydrogen = Fluid(FluidsList.Hydrogen).with_state(
        Input.pressure(coolant_inlet_pressure.value),  # Pa
        Input.temperature(coolant_inlet_temperature.value),  # K
    )
    channel = RectangularChannel(width=10e-3, height=5e-3, length=5.0, roughness=1e-5)
    channel = CircularChannel(diameter=10e-3, length=5.0, roughness=1e-5)
    coolant = Coolant(fluid=hydrogen, channel=channel, mass_flow=0.1)

    # --- Quick Tests ---
    # coolant.plot_coolant_properties()
    reynolds_number = coolant.get_reynolds_number()
    nusselt_number = coolant.get_nusselt_number_taylor(reynolds_number)
    heat_transfer_coefficient = coolant.get_heat_transfer_coefficient()
    fluid_speed = coolant.get_fluid_speed()
    print(f"Reynolds Number: {reynolds_number}")
    print(f"Nusselt Number: {nusselt_number}")
    print(f"Heat Transfer Coefficient: {heat_transfer_coefficient} W/m²·K")
    print(f"Fluid Speed: {fluid_speed} m/s")
