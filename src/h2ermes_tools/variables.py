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
        confidence (str): The confidence level of the variable.
        constraints (str): Constraints on the variable.
        margin (str): Margin of error for the variable.
    """

    name: str
    value: Any
    unit: str
    confidence: str
    constraints: str
    margin: str

    def __post_init__(self):
        # Set default values for confidence, constraints, and margin if not provided
        if not hasattr(self, "confidence"):
            self.confidence = "N/A"
        if not hasattr(self, "constraints"):
            self.constraints = "N/A"
        if not hasattr(self, "margin"):
            self.margin = "N/A"


# Enter your variables here
coolant_mass = Variable(
    name="Coolant Mass",
    value=0.0,
    unit="kg",
    confidence="poor",
    constraints="must be positive",
    margin="high",
)

heat_shield_thickness = Variable(
    name="Heat Shield Thickness",
    value=4e-3,
    unit="m",
    confidence="good",
    constraints="must be greater than 0",
    margin="10%",
)

if __name__ == "__main__":
    # Example of how to create a variable
    hydrogen_tank_diameter = Variable(
        name="Hydrogen Tank Diameter",
        value=0.0,
        unit="m/s",
        confidence="good",
        constraints="only positive values",
        margin="5%",
    )

    # Example of how to access the variable's attributes
    print(hydrogen_tank_diameter.name)
    print(hydrogen_tank_diameter.value)
