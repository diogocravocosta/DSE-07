import matplotlib.pyplot as plt
import numpy as np
from pyfluids import Fluid, FluidsList, Input
from h2ermes_tools.variables import coolant_inlet_pressure, coolant_inlet_temperature

coolant = Fluid(FluidsList.Hydrogen).with_state(
    Input.pressure(coolant_inlet_pressure.value),  # Pa
    Input.temperature(coolant_inlet_temperature.value),  # K
)


def plot_coolant_properties():
    temperatures = np.linspace(13.8, 300, 100)  # K
    pressures = np.linspace(1e5, 1e7, 100)  # Pa

    enthalpies = [
        coolant.with_state(
            Input.temperature(T), Input.pressure(coolant_inlet_pressure.value)
        ).enthalpy
        for T in temperatures
    ]
    density = [
        coolant.with_state(
            Input.temperature(T), Input.pressure(coolant_inlet_pressure.value)
        ).density
        for T in temperatures
    ]

    specific_heats = [
        coolant.with_state(
            Input.temperature(T), Input.pressure(coolant_inlet_pressure.value)
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
        f"Coolant Properties for Liquid Hydrogen at P= {round(coolant_inlet_pressure.value * 1e-5, 1)} bar"
    )

    plt.tight_layout()
    plt.show()


if __name__ == "__main__":
    plot_coolant_properties()
    print(
        f"Coolant Enthalpy at {coolant_inlet_temperature.value} K and {coolant_inlet_pressure.value} Pa: {coolant.enthalpy(coolant_inlet_temperature.value, coolant_inlet_pressure.value)} J/kg"
    )
    print(
        f"Coolant Density at {coolant_inlet_temperature.value} K and {coolant_inlet_pressure.value} Pa: {coolant.density(coolant_inlet_temperature.value, coolant_inlet_pressure.value)} kg/m³"
    )
