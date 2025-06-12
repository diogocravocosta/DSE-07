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
    value=3000.0,
    unit="kg",
    latex_symbol="M_c",
    confidence="poor",
    constraints="must be positive",
    margin=(1.0, 1.25),
)

heat_shield_thickness = Variable(
    name="Heat Shield Thickness",
    value=10e-3,  # includes the coolant channels
    unit="m",
    latex_symbol="t_{s}",
    confidence="good",
    constraints="must be greater than 0",
    margin=(0.8, 1.2),
)

heat_shield_mass = Variable(
    name="Heat Shield Mass",
    value=750.0,  # only the mass of the metallic structure
    unit="kg",
    latex_symbol="M_{hs}",
    confidence="poor",
    constraints="must be greater than zero",
    margin=(1.0, 1.5),
)

orbit_insertion_delta_v = Variable(
    name="Orbit Insertion Delta V",
    value=5800.,
    unit="m/s",
    latex_symbol="\\Delta V_{oi}",
    confidence="fair",
    constraints="must be positive",
    margin=(0.95, 1.05),
)

orbit_raising_delta_v = Variable(
    name="Orbit Raising Delta V",
    value=114.,
    unit="m/s",
    latex_symbol="\\Delta V_{or}",
    confidence="good",
    constraints="must be positive",
    margin=(1.0, 1.05),
)

target_orbit_circularization_delta_v = Variable(
    name="Target Orbit Circularization Delta V",
    value=113.,
    unit="m/s",
    latex_symbol="\\Delta V_{toc}",
    confidence="good",
    constraints="must be positive",
    margin=(1.0, 1.05),
)

deorbit_delta_v = Variable(
    name="Deorbit Delta V",
    value=160.,
    unit="m/s",
    latex_symbol="\\Delta V_{do}",
    confidence="good",
    constraints="must be positive",
    margin=(1.0, 1.05),
)

landing_delta_v = Variable(
    name="Landing Delta V",
    value=500.,
    unit="m/s",
    latex_symbol="\\Delta V_{l}",
    confidence="poor",
    constraints="must be positive",
    margin=(0.5, 1.05),
)

drag_compensation_delta_v = Variable(
    name="Drag Compensation Delta V",
    value=30.,
    unit="m/s",
    latex_symbol="\\Delta V_{dc}",
    confidence="poor",
    constraints="must be positive",
    margin=(0.5, 1.5),
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




diameter_exit_sl = Variable(
        name="Diameter of the Exit for Sea Level Optimized Combustion Chambers",
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
        value=8,
        unit="[-]",
        latex_symbol="N_{sl}",
        confidence="Design Parameter, can change",
        constraints="only positive values",
        margin=(1.0, 1.0),
    )

n_chambers_vac = Variable(
        name="The number of vacuum optimized thrust chambers used in the design",
        value=16,
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

gross_mass = Variable(
        name="The total mass of the vehicle including propellant",
        value=246_684.0317,
        unit="kg",
        latex_symbol="m_g",
        confidence="low",
        constraints="only positive values",
        margin=(0.8, 1.2),
    )

payload_mass = Variable(
        name="The payload mass delivered to depot",
        value=10_000,
        unit="kg",
        latex_symbol="m_p",
        confidence="set",
        constraints="only positive values",
        margin=(1.0, 1.0),
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

mass_flow_main = Variable(
        name="Mass Flow of the Main Engine, including all of the Thrust Chambers",
        value=514.4880, 
        unit="kg/s",
        latex_symbol="\dot{m}_{main}",
        confidence="result of other parameters, can change",
        constraints="only positive values",
        margin=(0.8, 1.2),
    )

t_w_sl = Variable(
        name="Thrust to Weight Ratio of the Vehicle during the start of the Landing Burn (only SL engines firing)",
        value=1.6, 
        unit="[-]",
        latex_symbol="T/W_{sl}",
        confidence="Result of Analysis, set value, will not change",
        constraints="only positive values",
        margin=(1, 1),
    )

t_w_vac = Variable(
        name="Thrust to Weight Ratio of the Vehicle during the start of the Vacuum Burn (SL and Vac engines firing)",
        value=1.2, 
        unit="[-]",
        latex_symbol="T/W_{vac}",
        confidence="Result of Analysis, set value, will not change",
        constraints="only positive values",
        margin=(1, 1),
    )

acs_peak_power = Variable(
    name="ACS Peak Power",
    value=2100,
    unit="W",
    latex_symbol="P_{acs,peak}",
    confidence="low",
    constraints="must be positive",
    margin=(0.5, 2.0),
)

acs_average_power = Variable(
    name="ACS Average Power",
    value=1200,
    unit="W",
    latex_symbol="P_{acs,avg}",
    confidence="low",
    constraints="must be positive",
    margin=(0.5, 2.0),
)

acs_dry_mass = Variable(
    name="ACS Dry Mass",
    value=150,
    unit="kg",
    latex_symbol="m_{acs}",
    confidence="low",
    constraints="must be positive",
    margin=(0.5, 2.0),
)

acs_propellant_mass = Variable(
    name="ACS Propellant Mass",
    value=800,
    unit="kg",
    latex_symbol="m_{acs,prop}",
    confidence="low",
    constraints="must be positive",
    margin=(0.5, 2.0),
)

docking_system_mass = Variable(
    name='Docking System Mass',
    value=370,
    unit='kg',
    latex_symbol='m_{docking system}',
    confidence='good',
    constraints='>0',
    margin=(1.0, 1.1),
)

minimum_tank_pressure = Variable(
    name="Minimum Tank Pressure",
    value=2e5,  # in Pa
    unit="Pa",
    latex_symbol="p_{min}",
    confidence="good",
    constraints="must be positive",
    margin=(1.0, 1.05),
)

nose_cone_mass = Variable(
    name="Nose Cone Mass",
    value=650,
    unit="kg",
    latex_symbol="m_{nose cone}",
    confidence="very low",
    constraints="must be positive",
    margin=(1.0, 1.5),
)

hydrogen_power_mass = Variable(
    name="Hydrogen Power Mass",
    value=225,
    unit="kg",
    latex_symbol="m_{H2, power}",
    confidence="low",
    constraints="must be positive",
    margin=(1.0, 1.5),
)

oxygen_power_mass = Variable(
    name="Oxygen Power Mass",
    value=225,
    unit="kg",
    latex_symbol="m_{O2, power}",
    confidence="low",
    constraints="must be positive",
    margin=(1.0, 1.5),
)

power_dry_mass = Variable(
    name="Power subsystem dry mass",
    value=719,
    unit="kg",
    latex_symbol="m_{power,dry}",
    confidence="midterm number",
    constraints="",
    margin=(1.0, 1.0)
)

data_handling_mass = Variable(
    name="Data handling subsystem mass",
    value=100,
    unit="kg",
    latex_symbol="m_{power,dry}",
    confidence="midterm number",
    constraints="",
    margin=(1.0, 1.0)
)

avionics_harness_mass = Variable(
    name="Avionics harness mass",
    value=403.32,
    unit="kg",
    latex_symbol="m_{avionics}",
    confidence="midterm number",
    constraints="",
    margin=(1.0, 1.0)
)

interstage_mass = Variable(
    name="Interstage mass",
    value=461.97,
    unit="kg",
    latex_symbol="m_{interstage}",
    confidence="midterm number",
    constraints="",
    margin=(1.0, 1.0)
)

gnc_mass = Variable(
    name="GNC mass",
    value=100,
    unit="kg",
    latex_symbol="m_{GNC}",
    confidence="midterm number",
    constraints="",
    margin=(1.0, 1.0)
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
