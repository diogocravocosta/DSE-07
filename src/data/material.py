from scipy import interpolate
import numpy as np
from matplotlib import pyplot as plt


class Material:
    """
    Represents a material with thermal and mechanical properties.

    Attributes:
        name (str): Name of the material.
        thermal_diffusivity (float): Thermal diffusivity [m^2/s].
        thermal_conductivity (float or callable): Thermal conductivity [W/m/K] or interpolated function.
        density (float): Density [kg/m^3].
        specific_heat (float): Specific heat [J/kg/K].
        poisson_ratio (float): Poisson's ratio.
        youngs_modulus (float): Young's modulus [Pa].
        thermal_expansion_coefficient (float): Thermal expansion coefficient [1/K].
        emissivity (float): Thermal emissivity.
        ultimate_strength (float): Ultimate strength [Pa].
        yield_strength (float): Yield strength [Pa].
        roughness_height (float): Effective roughness height [m].
        fracture_toughness (float): Fracture strength [MPa m^0.5].
        absorptivity (float): absorptivity [-]
    """

    def __init__(
        self,
        name=None,
        thermal_diffusivity=None,
        thermal_conductivity=None,
        density=None,
        specific_heat=None,
        poisson_ratio=None,
        youngs_modulus=None,
        thermal_expansion_coefficient=None,
        emissivity=None,
        ultimate_strength=None,
        yield_strength=None,
        roughness_height=None,
        fracture_toughness=None,
        absorptivity = None
    ):
        """
        Initialize a Material instance with its properties.
        """
        self.name = name
        self.alpha = thermal_diffusivity
        self.k = thermal_conductivity
        self.rho = density
        self.Cp = specific_heat
        self.nu = poisson_ratio
        self.E = youngs_modulus
        self.cte = thermal_expansion_coefficient
        self.eps = emissivity
        self.uts = ultimate_strength
        self.ys = yield_strength
        self.Ra = roughness_height
        self.Kic = fracture_toughness
        self.abs = absorptivity

    def set_thermal_conductivity(self, conductivity_values, temperature_values):
        """
        Set temperature-dependent thermal conductivity as an interpolated function.

        Args:
            conductivity_values (np.ndarray): Array of thermal conductivity values [W/m/K].
            temperature_values (np.ndarray): Array of corresponding temperatures [K].
        """
        self.k = interpolate.interp1d(
            temperature_values,
            conductivity_values,
            kind="linear",
            fill_value="extrapolate",
        )

    def plot_property(self, property_function, ylabel):
        """
        Plot a material property as a function of temperature.

        Args:
            property_function (callable): Function that takes temperature [K] and returns property value.
            ylabel (str): Label for the y-axis.
        """
        temperature_range = np.linspace(288, 1000, 100)
        plt.plot(temperature_range, property_function(temperature_range))
        plt.xlabel("Temperature [K]")
        plt.ylabel(ylabel)
        plt.title(self.name)
        plt.grid(True)
        plt.show()


# Ti-6Al-4V at room temperature
k = 6.7  # thermal conductivity [W/m/K]
Cp = 930  # specific heat [J/kg/K] at 870 C
rho = 4430  # density [kg/m^3]
v = 0.342  # Poisson ratio
E = 114e9  # Youngs modulus [Pa]
cte = 9.7e-6  # coefficinent of thermal expansion [1/K]
uts = 550e6  # ultimate strength [Pa]
eps = 0.8  # thermal emissivity
Ra = 2e-5  # effective roughness height [m] from Materialise
alpha = k / rho / Cp  # thermal diffusivity [m^2/s]
Ti6Al4V = Material("Ti-6Al-4V", alpha, k, rho, Cp, v, E, cte, eps, uts, Ra)

# set variable material properties from https://www.researchgate.net/publication/299647114_Developments_in_cutting_tool_technology_in_improving_machinability_of_Ti6Al4V_alloy_A_review
Ti6Al4V.set_thermal_conductivity(
    np.array([6.7, 9, 12, 15, 18]), np.array([283, 477, 700, 922, 1144])
)

# SS 1.4404 (SS 316L)
k = 23  # thermal conductivity [W/m/K] at 500 C
Cp = 590  # specific heat [J/kg/K] at 500 C
rho = 7750  # density [kg/m^3] at 500 C
v = 1  # poiosson ratio
E = 1  # youngs modulus [Pa]
cte = 1  # thermal expasion coefficinent [1/K]
uts = 1  # ultimate strength [Pa]
eps = 0.8  # thermal emissivity
Ra = 0.9e-5  # effective roughness height [m] from Materialize
alpha = k / rho / Cp  # thermal diffusivity [m^2/s]
SS14404 = Material("SS14404", alpha, k, rho, Cp, v, E, cte, eps, uts, Ra)

# variable material properties from 'Transient thermo-mechanical modeling of stress evolution and re-melt volume fraction in electron beam additive manufacturing process' by R.K. Adhitan
SS14404.set_thermal_conductivity(
    np.array([13.5, 16, 17.5, 19.5, 22, 23.5, 24.5, 25.5, 27, 28, 29]),
    np.array([273, 373, 473, 573, 673, 773, 873, 973, 1073, 1173, 1273]),
)

## SS 1.4845 (SS 310S)
k = 16.2  # thermal conductivity [W/m/K] at 500 C
Cp = 500  # specific heat [J/kg/K] at 500 C
rho = 7900  # density [kg/m^3] at 500 C
v = 0.3  # Poisson ratio
E = 1  # Youngs modulus [Pa]
cte = 1  # thermal expansion coefficient [1/K]
uts = 1  # ultimate strength [Pa]
eps = 0.8  # thermal emissivity
Ra = 0.9e-5  # effective roughness height [m] from Materialise
alpha = k / rho / Cp  # thermal diffusivity [m^2/s]
SS14845 = Material("SS14845", alpha, k, rho, Cp, v, E, cte, eps, uts, Ra)


if __name__ == "__main__":
    # Example usage
    Ti6Al4V.plot_property(Ti6Al4V.k, "Thermal Conductivity [W/m/K]")
    SS14404.plot_property(SS14404.k, "Thermal Conductivity [W/m/K]")
    print(Ti6Al4V.k(500))  # Get thermal conductivity at 500 K
    print(SS14404.k(500))  # Get thermal conductivity at 500 K
