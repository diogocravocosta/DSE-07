"""
THIS IS THE MAIN FILE FOR VARIABLES
This file contains the main variables used in the H2ERMES vehicle.
It includes the definition of the variables, their types, units and default values.
It also tracks their confidence, constraints, and margins.
"""

from dataclasses import dataclass
from typing import Any


@dataclass
class Variable:
    """
    A class to represent a variable in the H2ERMES vehicle.

    Attributes:
        name (str): The name of the variable.
        value (Any): The value of the variable.
        unit (str): The unit of the variable.
        latex-symbol (str): The LaTeX symbol for the variable (optional).
        confidence (str): Qualitative assessment of the accuracy of the variable value (optional).
        constraints (str): Constraints on the variable (optional).
        margin (tuple): The bottom and upper margin represented as multipliers, default is (1.0, 1.0).
    """

    name: str
    value: Any
    unit: str
    latex_symbol: str
    confidence: str
    constraints: str
    margin: tuple = (1.0, 1.0)

    def __post_init__(self):
        # Set default values for confidence, constraints, and margin if not provided
        if not hasattr(self, "latex_symbol"):
            self.latex_symbol = self.name.replace(" ", "_").lower()
        if not hasattr(self, "confidence"):
            self.confidence = "N/A"
        if not hasattr(self, "constraints"):
            self.constraints = "N/A"


# Enter your variables here
coolant_mass = Variable(
    name="Coolant Mass",
    value=0.0,
    unit="kg",
    latex_symbol="M_c",
    confidence="poor",
    constraints="must be positive",
    margin=(1.0, 1.25),
)

heat_shield_thickness = Variable(
    name="Heat Shield Thickness",
    value=4e-3,
    unit="m",
    latex_symbol="t_{s}",
    confidence="good",
    constraints="must be greater than 0",
    margin=(1.0, 1.1),
)

launch_vehicle_dimensions = Variable(
    name="Launch Vehicle Dimensions",
    value=[15,10,5], # Rectangular dimensions in meters (length, width, height); Needs to be updated with actual shape
    unit="m",
    latex_symbol="lv_{dim}",
    confidence="highly approximate", #Update with actual dimensions
    constraints="must be positive",
    margin=(1.0, 1.05),
)

mass_vehicle= Variable(
    name="Vehicle Mass",
    value=150000, # Approximate value for H2ermes vehicle
    unit="kg",
    latex_symbol="M_{vehicle}",
    confidence="highly approximate", #Needs to be updated with actual values
    constraints="must be positive",
    margin=(1.0, 1.1),
)

MMOI_vehicle = Variable(
    name="Mass Moment of Inertia",
    value=[625000, 1250000, 1625000], # Approximate values for chosen shape
    unit="kg*m^2",
    latex_symbol="I_{mass}",
    confidence="highly approximate", #Needs to be updated with actual values
    constraints="must be positive",
    margin=(1.0, 1.1),
)

vehicle_surface_area = Variable(
    name="Vehicle Surface Area",
    value=150,  # Placeholder value, needs to be calculated based on actual dimensions
    unit="m^2",
    latex_symbol="A_{vehicle}",
    confidence="approximate",
    constraints="must be positive",
    margin=(1.0, 1.1),
)

if __name__ == "__main__":
    # Example of how to create a variable
    hydrogen_tank_diameter = Variable(
        name="Hydrogen Tank Diameter",
        value=0.0,
        unit="m/s",
        latex_symbol="D_{t_{LH2}}",
        confidence="good",
        constraints="only positive values",
        margin=(1.0, 1.05),
    )

    # Example of how to access the variable's attributes
    print(hydrogen_tank_diameter.name)
    print(hydrogen_tank_diameter.value)
