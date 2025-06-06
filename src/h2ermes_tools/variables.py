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

engine_chamber_pressure = Variable(
        name="Engine Chamber Pressure",
        value=6100000.0,
        unit="Pa",
        latex_symbol="P_{c}",
        confidence="design parameter, should not change",
        constraints="only positive values",
        margin=(1.0, 1.0),
    )

diameter_combustion_chamber_sl = Variable(
        name="Diameter of the Combustion Chamber for Seal Level Optimized Combustion Chambers",
        value=0.17038,
        unit="m",
        latex_symbol="D_{cc_{sl}}",
        confidence="good, design can change",
        constraints="only positive values",
        margin=(1.0, 1.0),
    )

diameter_combustion_chamber_vac = Variable(
        name="Diameter of the Combustion Chamber for Vacuum Optimized Combustion Chambers",
        value=0.15966,
        unit="m",
        latex_symbol="D_{cc_{vac}}",
        confidence="good, design can change",
        constraints="only positive values",
        margin=(1.0, 1.0),
    )

UTS_Inconel_718 = Variable(
        name="Ultimate Tensile Strength of Inconel 718, the material used for the combustion chambers, at 900 k, the max temperature for the inner walls in the firebolt design report",
        value=1151424468,
        unit="Pa",
        latex_symbol="\sigma_{UTS_{Inconel718}}",
        confidence="good, depends on treatment",
        constraints="only positive values",
        margin=(0.95, 1.05),
    )

density_Inconel_718 = Variable(
        name="Desnity of Inconel 718, the material used for the combustion chambers, at 900 k, the max temperature for the inner walls in the firebolt design report",
        value=8220.9316989,
        unit="Kg/m^3",
        latex_symbol="\rho_{Inconel718}",
        confidence="good, depends on treatment",
        constraints="only positive values",
        margin=(1.0, 1.01),
    )

diameter_throat_sl = Variable(
        name="Diameter of the Throat for Seal Level Optimized Combustion Chambers",
        value=0.09133,
        unit="m",
        latex_symbol="D_{t_{sl}}",
        confidence="good, design can change",
        constraints="only positive values",
        margin=(1.0, 1.0),
    )

diameter_throat_vac = Variable(
        name="Diameter of the Throat for Vacuum Optimized Combustion Chambers",
        value=0.08559,
        unit="m",
        latex_symbol="D_{t_{vac}}",
        confidence="good, design can change",
        constraints="only positive values",
        margin=(1.0, 1.0),
    )

contraction_angle_combustion_chamber_vac = Variable(
        name="The contraction angle of the combustion chamber for vacuum optimized combustion chambers",
        value=0.5235988, #basically 30 degrees in radians
        unit="Rad",
        latex_symbol="\theta_{cc_{vac}}",
        confidence="good, design parameter",
        constraints="only positive values",
        margin=(1.0, 1.0),
    )

contraction_angle_combustion_chamber_sl = Variable(
        name="The contraction angle of the combustion chamber for sea level optimized combustion chambers",
        value=0.5235988, #basically 30 degrees in radians
        unit="Rad",
        latex_symbol="\theta_{cc_{sl}}",
        confidence="good, design parameter",
        constraints="only positive values",
        margin=(1.0, 1.0),
    )

length_combustion_chamber_sl = Variable(
        name="Length of the Combustion Chamber for Seal Level Optimized Combustion Chambers",
        value=0.33444,
        unit="m",
        latex_symbol="L_{cc_{sl}}",
        confidence="good, design can change",
        constraints="only positive values",
        margin=(1.0, 1.0),
    )

length_combustion_chamber_vac = Variable(
        name="Length of the Combustion Chamber for Vacuum Optimized Combustion Chambers",
        value=0.33184,
        unit="m",
        latex_symbol="L_{cc_{vac}}",
        confidence="good, design can change",
        constraints="only positive values",
        margin=(1.0, 1.0),
    )

length_nozzle_sl = Variable(
        name="Length of the nozzle for Seal Level Optimized Combustion Chambers",
        value=0.23551,
        unit="m",
        latex_symbol="L_{n_{sl}}",
        confidence="good, design can change",
        constraints="only positive values",
        margin=(1.0, 1.0),
    )

length_nozzle_vac = Variable(
        name="Length of the nozzle for Vacuum Optimized Combustion Chambers",
        value=0.95595,
        unit="m",
        latex_symbol="L_{n_{vac}}",
        confidence="good, design can change",
        constraints="only positive values",
        margin=(1.0, 1.0),
    )

theta_cn_sl = Variable(
        name="Half Angle of the nozzle for Seal Level Optimized Combustion Chambers",
        value=0.40544,
        unit="Rad",
        latex_symbol="\theta_{cn_{sl}}",
        confidence="good, design can change",
        constraints="only positive values",
        margin=(1.0, 1.0),
    )

theta_cn_vac = Variable(
        name="Half Angle of the nozzle for Vacuum Optimized Combustion Chambers",
        value=0.670381,
        unit="Rad",
        latex_symbol="\theta_{cn_{vac}}",
        confidence="good, design can change",
        constraints="only positive values",
        margin=(1.0, 1.0),
    )


diameter_exit_sl = Variable(
        name="Diameter of the Exit for Seal Level Optimized Combustion Chambers",
        value=0.25833,
        unit="m",
        latex_symbol="D_{e_{sl}}",
        confidence="good, design can change",
        constraints="only positive values",
        margin=(1.0, 1.0),
    )

diameter_exit_vac = Variable(
        name="Diameter of the Exit for Vacuum Optimized Combustion Chambers",
        value=0.76550,
        unit="m",
        latex_symbol="D_{e_{vac}}",
        confidence="good, design can change",
        constraints="only positive values",
        margin=(1.0, 1.0),
    )

engine_mixture_ratio = Variable(
        name="The ratio of mass of oxidizer to mass of fuel in the engine",
        value=6.0,
        unit="[-]",
        latex_symbol="O/F",
        confidence="design parameter",
        constraints="only positive values",
        margin=(1.0, 1.0),
    )

engine_specific_impulse_sl_opt = Variable(
        name="The specific impulse of the sea level nozzles in sea level (optimal) conditions",
        value=361.6451,
        unit="s",
        latex_symbol="I_{sp_{sl_{OPT}}}",
        confidence="Result of the Desgin, can change",
        constraints="only positive values",
        margin=(1.0, 1.0),
    )

engine_specific_impulse_sl_non_opt = Variable(
        name="The specific impulse of the ses level nozzles in vacuum (non-optimal) conditions",
        value=393.3471,
        unit="s",
        latex_symbol="I_{sp_{sl_{NON-OPT}}}",
        confidence="Result of the Desgin, can change",
        constraints="only positive values",
        margin=(1.0, 1.0),
    )

engine_specific_impulse_vac_opt = Variable(
        name="The specific impulse of the vacuum nozzles in vacuum (optimal) conditions",
        value=447.9481,
        unit="s",
        latex_symbol="I_{sp_{sl_{OPT}}}",
        confidence="Result of the Desgin, can change",
        constraints="only positive values",
        margin=(1.0, 1.0),
    )

engine_specific_impulse_vac_non_opt = Variable(
        name="The specific impulse of the vacuum nozzles in sea level (non-optimal) conditions",
        value=130.9458,
        unit="s",
        latex_symbol="I_{sp_{sl_{NON-OPT}}}",
        confidence="Result of the Desgin, can change",
        constraints="only positive values",
        margin=(1.0, 1.0),
    )

engine_thrust_sl_opt = Variable(
        name="The thrust of each of the sea level nozzles in sea level (optimal) conditions",
        value=61364.9,
        unit="N",
        latex_symbol="T_{sl_{OPT}}",
        confidence="Design Parameter, can change",
        constraints="only positive values",
        margin=(1.0, 1.0),
    )

engine_thrust_sl_non_opt = Variable(
        name="The thrust of each of the sea level nozzles in vacuum (non-optimal) conditions",
        value=66744.2,
        unit="N",
        latex_symbol="T_{sl_{NON-OPT}}",
        confidence="Result of the Desgin, can change",
        constraints="only positive values",
        margin=(1.0, 1.0),
    )

engine_thrust_vac_opt = Variable(
        name="The thrust of each of the vacuum nozzles in vacuum (optimal) conditions",
        value=66744.2,
        unit="N",
        latex_symbol="T_{vac_{OPT}}",
        confidence="Design Parameter, can change",
        constraints="only positive values",
        margin=(1.0, 1.0),
    )

engine_thrust_vac_non_opt = Variable(
        name="The the thrust of each of the vacuum nozzles in sea level (non-optimal) conditions",
        value=19510.9,
        unit="N",
        latex_symbol="T_{vac_{NON-OPT}}",
        confidence="Result of the Desgin, can change",
        constraints="only positive values",
        margin=(1.0, 1.0),
    )

expansion_ratio_vac = Variable(
        name="The expansion ratio of each of the vacuum nozzles (ratio of the area of the exit and the area of the throat)",
        value=80,
        unit="[-]",
        latex_symbol="A_e/A_t_{vac}",
        confidence="Design Parameter, can change",
        constraints="only positive values",
        margin=(1.0, 1.0),
    )

expansion_ratio_sl = Variable(
        name="The expansion ratio of each of the sea level nozzles (ratio of the area of the exit and the area of the throat)",
        value=8,
        unit="[-]",
        latex_symbol="A_e/A_t_{sl}",
        confidence="Design Parameter, can change",
        constraints="only positive values",
        margin=(1.0, 1.0),
    )

n_chambers_sl = Variable(
        name="The number of sea level thrust chambers used in the design",
        value=11,
        unit="[-]",
        latex_symbol="N_{sl}",
        confidence="Design Parameter, can change",
        constraints="only positive values",
        margin=(1.0, 1.0),
    )

n_chambers_vac = Variable(
        name="The number of vacuum optimized thrust chambers used in the design",
        value=13,
        unit="[-]",
        latex_symbol="N_{vac}",
        confidence="Design Parameter, can change",
        constraints="only positive values",
        margin=(1.0, 1.0),
    )

total_dry_mass = Variable(
        name="The total dry mass of the vehicle",
        value=20642.21346,
        unit="kg",
        latex_symbol="m_s",
        confidence="low",
        constraints="only positive values",
        margin=(0.8, 1.2),
    )

thrust_chamber_sl_mass = Variable(
        name="The mass of a sea level optimized Thrust Chamber",
        value=4.063444786849481, 
        unit="kg",
        latex_symbol="m_{tc_{sl}}",
        confidence="medium, conservative value",
        constraints="only positive values",
        margin=(0.8, 1.2),
    )

thrust_chamber_vac_mass = Variable(
        name="The mass of a vacuum optimized Thrust Chamber",
        value=19.666088892712068, 
        unit="kg",
        latex_symbol="m_{tc_{vac}}",
        confidence="medium, conservative value",
        constraints="only positive values",
        margin=(0.8, 1.2),
    )

l_star = Variable(
        name="Thrust Chamber Characteristic Length",
        value=1.02,
        unit="m",
        latex_symbol="L*",
        confidence="design parameter, most conservative for liquid hydrogen",
        constraints="only positive values",
        margin=(0.8, 1.2),
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
