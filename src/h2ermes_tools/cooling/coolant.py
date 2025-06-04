import matplotlib.pyplot as plt
import numpy as np
from pyfluids import Fluid, FluidsList, Input
from h2ermes_tools.variables import coolant_inlet_pressure, coolant_inlet_temperature
from h2ermes_tools.cooling.channel import CircularChannel


class Coolant:
    """
    Represents a coolant with properties derived from the pyfluids library.

    Attributes:
        fluid (Fluid): The fluid object from pyfluids representing the coolant.
        channel (Channel): The channel through which the coolant flows.
        mass_flow (float): The mass flow rate of the coolant in kg/s.
    """

    def __init__(self, fluid, channel, mass_flow):
        self.fluid = fluid
        self.channel = channel
        self.mass_flow = mass_flow

    def enthalpy(self, temperature, pressure):
        return self.fluid.with_state(
            Input.temperature(temperature), Input.pressure(pressure)
        ).enthalpy

    def density(self, temperature, pressure):
        return self.fluid.with_state(
            Input.temperature(temperature), Input.pressure(pressure)
        ).density

    def specific_heat(self, temperature, pressure):
        return self.fluid.with_state(
            Input.temperature(temperature), Input.pressure(pressure)
        ).specific_heat

    def plot_coolant_properties(
        self, temperatures=np.linspace(13.8, 300, 100), pressure=1e6
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

        fig, axs = plt.subplots(3, 1, figsize=(15, 7))
        axs[0].plot(temperatures, enthalpies, label="Enthalpy")
        # axs[0].set_title("Coolant Enthalpy vs Temperature")
        # axs[0].set_xlabel("Temperature (K)")
        axs[0].set_ylabel("Enthalpy (J/kg)")
        axs[0].grid()
        axs[0].legend()

        axs[1].plot(temperatures, density, label="Density", color="orange")
        # axs[1].set_title("Coolant Density vs Temperature")
        # axs[1].set_xlabel("Temperature (K)")
        axs[1].set_ylabel("Density (kg/m³)")
        axs[1].grid()
        axs[1].legend()

        axs[2].plot(temperatures, specific_heats, label="Specific Heat", color="green")
        # axs[2].set_title("Coolant Specific Heat vs Temperature")
        axs[2].set_xlabel("Temperature (K)")
        axs[2].set_ylabel("Specific Heat (J/kg·K)")
        axs[2].grid()
        axs[2].legend()

        plt.suptitle(
            f"Coolant Properties for Hydrogen at P = {round(coolant_inlet_pressure.value * 1e-5, 1)} bar"
        )

        plt.tight_layout()
        plt.show()

    def get_nusselt_number_taylor(
        self, reynolds_number, temperature_ratio=0.55, dimensionless_length=1.0
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
            * (reynolds_number**0.8)
            * (self.fluid.prandtl**0.4)
            * (temperature_ratio ** (-0.57 - (1.59 / dimensionless_length)))
        )

    def get_reynolds_number(self, fluid_speed) -> float:
        """
        Calculate the Reynolds number.

        Parameters:
        - fluid: Fluid object from pyfluids
        - fluid_speed: Speed of the fluid (m/s)

        Returns:
        - Reynolds number (dimensionless)
        """
        return (
            self.fluid.density * fluid_speed * self.channel.diameter
        ) / self.fluid.dynamic_viscosity


if __name__ == "__main__":
    coolant = Fluid(FluidsList.Hydrogen).with_state(
        Input.pressure(coolant_inlet_pressure.value),  # Pa
        Input.temperature(coolant_inlet_temperature.value),  # K
    )
    coolant = Coolant(fluid=coolant, channel=CircularChannel(1e-3, 5, 1e-5))

    # quick tests
    # coolant.plot_coolant_properties()
    reynolds_number = coolant.get_reynolds_number(10.0)  # Example fluid speed of 1 m/s
    nusselt_number = coolant.get_nusselt_number_taylor(reynolds_number)
    print(f"Reynolds Number: {reynolds_number}")
    print(f"Nusselt Number: {nusselt_number}")
