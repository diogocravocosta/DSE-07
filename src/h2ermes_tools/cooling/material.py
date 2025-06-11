from scipy import interpolate
from matplotlib import pyplot as plt
import numpy as np
import tomllib

# Load material properties from TOML files
with open("src/data/materials/ss310.toml", "rb") as f:
    ss310_data = tomllib.load(f)

with open("src/data/materials/ss304.toml", "rb") as f:
    ss304_data = tomllib.load(f)


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
        roughness_height (float): Effective roughness height [m].
        maximum_temperature (float): Maximum operating temperature for the material [K].
    """

    def __init__(
        self,
        name,
        thermal_diffusivity,
        thermal_conductivity,
        density,
        specific_heat,
        poisson_ratio,
        youngs_modulus,
        thermal_expansion_coeffient,
        emissivity,
        ultimate_strength,
        roughness_height,
        maximum_temperature,
    ):
        """
        Initialize a Material instance with its properties.
        """
        self.name = name
        self.thermal_diffusivity = thermal_diffusivity
        self.thermal_conductivity = thermal_conductivity
        self.density = density
        self.specific_heat = specific_heat
        self.poisson_ratio = poisson_ratio
        self.youngs_modulus = youngs_modulus
        self.thermal_expansion_coeffient = thermal_expansion_coeffient
        self.emissivity = emissivity
        self.ultimate_strength = ultimate_strength
        self.roughness_height = roughness_height
        self.maximum_temperature = maximum_temperature

    def set_specific_heat(self, temperature_values, specific_heat_values):
        """
        Set temperature-dependent specific heat as an interpolated function.

        Args:
            temperature_values (np.ndarray): Array of temperatures [K].
            specific_heat_values (np.ndarray): Array of specific heat values [J/kg/K].
        """
        self.specific_heat = interpolate.interp1d(
            temperature_values,
            specific_heat_values,
            kind="linear",
            fill_value="extrapolate",
        )

    def set_thermal_conductivity(self, temperature_values, conductivity_values):
        """
        Set temperature-dependent thermal conductivity as an interpolated function.

        Args:
            conductivity_values (np.ndarray): Array of thermal conductivity values [W/m/K].
            temperature_values (np.ndarray): Array of corresponding temperatures [K].
        """
        self.thermal_conductivity = interpolate.interp1d(
            temperature_values,
            conductivity_values,
            kind="linear",
            fill_value="extrapolate",
        )

    def set_thermal_expansion_coefficient(self, temperature_values, cte_values):
        """
        Set temperature-dependent thermal expansion coefficient as an interpolated function.

        Args:
            temperature_values (np.ndarray): Array of temperatures [K].
            cte_values (np.ndarray): Array of thermal expansion coefficient values [1/K].
        """
        self.thermal_expansion_coeffient = interpolate.interp1d(
            temperature_values,
            cte_values,
            kind="linear",
            fill_value="extrapolate",
        )

    def set_yield_strength(self, temperature_values, yield_strength_values):
        """
        Set temperature-dependent yield strength as an interpolated function.

        Args:
            temperature_values (np.ndarray): Array of temperatures [K].
            yield_strength_values (np.ndarray): Array of yield strength values [Pa].
        """
        self.ultimate_strength = interpolate.interp1d(
            temperature_values,
            yield_strength_values,
            kind="linear",
            fill_value="extrapolate",
        )

    def set_youngs_modulus(self, temperature_values, youngs_modulus_values):
        """
        Set temperature-dependent Young's modulus as an interpolated function.

        Args:
            temperature_values (np.ndarray): Array of temperatures [K].
            youngs_modulus_values (np.ndarray): Array of Young's modulus values [Pa].
        """
        self.youngs_modulus = interpolate.interp1d(
            temperature_values,
            youngs_modulus_values,
            kind="linear",
            fill_value="extrapolate",
        )

    def set_thermal_diffusivity(self, temperature_values):
        """
        Set thermal diffusivity based on the material's thermal conductivity, density, and specific heat.

        Args:
            temperature_values (np.ndarray): Array of temperatures [K].
        """
        k_values = self.thermal_conductivity(temperature_values)
        Cp_values = self.specific_heat(temperature_values)

        self.thermal_diffusivity = k_values / (self.density * Cp_values)
        self.thermal_diffusivity = interpolate.interp1d(
            temperature_values,
            self.thermal_diffusivity,
            kind="linear",
            fill_value="extrapolate",
        )

    def plot_property(self, temperature_range, property_function, ylabel):
        """
        Plot a material property as a function of temperature.

        Args:
            temperature_range (np.ndarray): Array of temperatures [K].
            property_function (callable): Function that takes temperature [K] and returns property value.
            ylabel (str): Label for the y-axis.
        """
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
maxT = 400 + 273.15  # maximum operating temperature for the material [K]
Ti6Al4V = Material("Ti-6Al-4V", alpha, k, rho, Cp, v, E, cte, eps, uts, Ra, maxT)

# set variable material properties from https://www.researchgate.net/publication/299647114_Developments_in_cutting_tool_technology_in_improving_machinability_of_Ti6Al4V_alloy_A_review
Ti6Al4V.set_thermal_conductivity(
    np.array([6.7, 9, 12, 15, 18]), np.array([283, 477, 700, 922, 1144])
)

# SS 316L (SS 1.4404)
k = 23  # thermal conductivity [W/m/K] at 500 C
Cp = 590  # specific heat [J/kg/K] at 500 C
rho = 7750  # density [kg/m^3] at 500 C
v = 0.3  # Poisson ratio
E = 193e9  # youngs modulus [Pa]
cte = 15.9e-6  # thermal expansion coefficient [1/K]
uts = 515e6  # ultimate strength [Pa]
eps = 0.5  # thermal emissivity
Ra = 0.9e-5  # effective roughness height [m] from Materialize
alpha = k / rho / Cp  # thermal diffusivity [m^2/s]
maxT = 925 + 273.15  # maximum operating temperature for the material [K]
SS316L = Material("SS316L", alpha, k, rho, Cp, v, E, cte, eps, uts, Ra, maxT)

# variable material properties from 'Transient thermo-mechanical modeling of stress evolution and re-melt volume fraction in electron beam additive manufacturing process' by R.K. Adhitan
SS316L.set_thermal_conductivity(
    np.array([13.5, 16, 17.5, 19.5, 22, 23.5, 24.5, 25.5, 27, 28, 29]),
    np.array([273, 373, 473, 573, 673, 773, 873, 973, 1073, 1173, 1273]),
)

# SS 310 (SS 1.4845)
k = ss310_data["thermal_conductivity"]
Cp = ss310_data["specific_heat"]
rho = ss310_data["density"]
v = ss310_data["poisson_ratio"]  # Poisson ratio
E = ss310_data["youngs_modulus"]  # Young's modulus [Pa]
cte = ss310_data["thermal_expansion_coefficient"]  # thermal expansion coefficient [1/K]
uts = ss310_data["ultimate_strength"]  # ultimate strength [Pa]
eps = ss310_data["emissivity"]  # thermal emissivity
Ra = ss310_data["roughness_height"]  # effective roughness height [m]
maxT = ss310_data[
    "maximum_temperature"
]  # maximum operating temperature for the material [K]
alpha = k / rho / Cp  # thermal diffusivity [m^2/s]
SS310 = Material("SS310", alpha, k, rho, Cp, v, E, cte, eps, uts, Ra, maxT)

cp_data = np.genfromtxt(
    "src/data/materials/ss310-specific-heat.csv", delimiter=",", skip_header=1
)
k_data = np.genfromtxt(
    "src/data/materials/ss310-thermal-conductivity.csv", delimiter=",", skip_header=1
)
cte_data = np.genfromtxt(
    "src/data/materials/ss310-thermal-expansion-coefficient.csv",
    delimiter=",",
    skip_header=1,
)
yield_strength_data = np.genfromtxt(
    "src/data/materials/ss310-yield-strength.csv", delimiter=",", skip_header=1
)
youngs_modulus_data = np.genfromtxt(
    "src/data/materials/ss310-youngs-modulus.csv", delimiter=",", skip_header=1
)

SS310.set_specific_heat(cp_data[:, 0] + 273.15, cp_data[:, 1])
SS310.set_thermal_conductivity(k_data[:, 0] + 273.15, k_data[:, 1])
SS310.set_thermal_expansion_coefficient(cte_data[:, 0] + 273.15, cte_data[:, 1])
SS310.set_yield_strength(yield_strength_data[:, 0] + 273.15, yield_strength_data[:, 1])
SS310.set_youngs_modulus(youngs_modulus_data[:, 0] + 273.15, youngs_modulus_data[:, 1])
SS310.set_thermal_diffusivity(np.linspace(10, 1273, 100))

# SS304 (SS 1.4301)
k = ss304_data["thermal_conductivity"]
Cp = ss304_data["specific_heat"]
rho = ss304_data["density"]
v = ss304_data["poisson_ratio"]  # Poisson ratio
E = ss304_data["youngs_modulus"]  # Young's modulus [Pa]
cte = ss304_data["thermal_expansion_coefficient"]  # thermal expansion coefficient [1/K]
uts = ss304_data["ultimate_strength"]  # ultimate strength [Pa]
eps = ss304_data["emissivity"]  # thermal emissivity
Ra = ss304_data["roughness_height"]  # effective roughness height [m]
maxT = ss304_data[
    "maximum_temperature"
]  # maximum operating temperature for the material [K]
alpha = k / rho / Cp  # thermal diffusivity [m^2/s]
SS304 = Material("SS304", alpha, k, rho, Cp, v, E, cte, eps, uts, Ra, maxT)

cp_data = np.genfromtxt(
    "src/data/materials/ss304-specific-heat.csv", delimiter=",", skip_header=1
)
k_data = np.genfromtxt(
    "src/data/materials/ss304-thermal-conductivity.csv", delimiter=",", skip_header=1
)
cte_data = np.genfromtxt(
    "src/data/materials/ss304-thermal-expansion-coefficient.csv",
    delimiter=",",
    skip_header=1,
)
yield_strength_data = np.genfromtxt(
    "src/data/materials/ss304-yield-strength.csv", delimiter=",", skip_header=1
)
youngs_modulus_data = np.genfromtxt(
    "src/data/materials/ss304-youngs-modulus.csv", delimiter=",", skip_header=1
)
SS304.set_specific_heat(cp_data[:, 0] + 273.15, cp_data[:, 1])
SS304.set_thermal_conductivity(k_data[:, 0] + 273.15, k_data[:, 1])
SS304.set_thermal_expansion_coefficient(cte_data[:, 0] + 273.15, cte_data[:, 1])
SS304.set_yield_strength(yield_strength_data[:, 0] + 273.15, yield_strength_data[:, 1])
SS304.set_youngs_modulus(youngs_modulus_data[:, 0] + 273.15, youngs_modulus_data[:, 1])
SS304.set_thermal_diffusivity(np.linspace(10, 1273, 100))


if __name__ == "__main__":
    # Example usage
    # Ti6Al4V.plot_property(Ti6Al4V.thermal_conductivity, "Thermal Conductivity [W/m/K]")
    # SS316L.plot_property(SS316L.thermal_conductivity, "Thermal Conductivity [W/m/K]")
    # print(Ti6Al4V.thermal_conductivity(500))  # Get thermal conductivity at 500 K
    # print(SS316L.thermal_conductivity(500))  # Get thermal conductivity at 500 K

    temperature_range = np.linspace(
        10, 1273, 100
    )  # Temperature range from 300 K to 1273 K
    SS310.plot_property(
        temperature_range, SS310.specific_heat, "Specific Heat [J/kg/K]"
    )
    SS310.plot_property(
        temperature_range, SS310.thermal_conductivity, "Thermal Conductivity [W/m/K]"
    )
    SS310.plot_property(
        temperature_range,
        SS310.thermal_expansion_coeffient,
        "Thermal Expansion Coefficient [1/K]",
    )
    SS310.plot_property(
        temperature_range, SS310.ultimate_strength, "Ultimate Strength [Pa]"
    )
    SS310.plot_property(temperature_range, SS310.youngs_modulus, "Young's Modulus [Pa]")
    SS310.plot_property(
        temperature_range, SS310.thermal_diffusivity, "Thermal Diffusivity [m^2/s]"
    )

    SS304.plot_property(
        temperature_range, SS304.specific_heat, "Specific Heat [J/kg/K]"
    )
    SS304.plot_property(
        temperature_range, SS304.thermal_conductivity, "Thermal Conductivity [W/m/K]"
    )
    SS304.plot_property(
        temperature_range,
        SS304.thermal_expansion_coeffient,
        "Thermal Expansion Coefficient [1/K]",
    )
    SS304.plot_property(
        temperature_range, SS304.ultimate_strength, "Ultimate Strength [Pa]"
    )
    SS304.plot_property(temperature_range, SS304.youngs_modulus, "Young's Modulus [Pa]")
