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

orbit_insertion_delta_v = Variable(
    name="Orbit Insertion Delta V",
    value=5800.,
    unit="m/s",
    latex_symbol="\Delta V_{oi}",
    confidence="fair",
    constraints="must be positive",
    margin=(0.95, 1.05),
)

orbit_raising_delta_v = Variable(
    name="Orbit Raising Delta V",
    value=114.,
    unit="m/s",
    latex_symbol="\Delta V_{or}",
    confidence="good",
    constraints="must be positive",
    margin=(1.0, 1.05),
)

target_orbit_circularization_delta_v = Variable(
    name="Target Orbit Circularization Delta V",
    value=113.,
    unit="m/s",
    latex_symbol="\Delta V_{toc}",
    confidence="good",
    constraints="must be positive",
    margin=(1.0, 1.05),
)

deorbit_delta_v = Variable(
    name="Deorbit Delta V",
    value=160.,
    unit="m/s",
    latex_symbol="\Delta V_{do}",
    confidence="good",
    constraints="must be positive",
    margin=(1.0, 1.05),
)

landing_delta_v = Variable(
    name="Landing Delta V",
    value=500.,
    unit="m/s",
    latex_symbol="\Delta V_{l}",
    confidence="poor",
    constraints="must be positive",
    margin=(0.5, 1.05),
)

drag_compensation_delta_v = Variable(
    name="Drag Compensation Delta V",
    value=30.,
    unit="m/s",
    latex_symbol="\Delta V_{dc}",
    confidence="poor",
    constraints="must be positive",
    margin=(0.5, 1.5),
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
