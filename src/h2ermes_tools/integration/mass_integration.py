import numpy as np

import data.constants as cn
import h2ermes_tools.integration.dummy_dry_mass as dds
from h2ermes_tools.propulsion.cycle_sizing import size_turbopump
from h2ermes_tools.landinglegs import size_landing_legs
from h2ermes_tools.structures.tank_sizing import size_tanks
from h2ermes_tools.fatigue_calculation_tool import thickness_optimization_fatigue

tank_material = 1
class MassIntegrator:
    """
    Class to handle mass integration

    dry_mass: float in kg

    gross_mass: float in kg

    propellant_mass: float, mass of propellant used for primary propulsion in kg
        landing_propellant_mass: float, mass of propellant used for landing in kg
        deorbit_propellant_mass: float, mass of propellant used for deorbit in kg
        transfer_propellant_mass: float, mass of propellant used for transfer to depo in kg
        orbit_insertion_propellant_mass: float, mass of propellant used for orbit insertion in kg

    payload_mass: float, payload mass to be transferred to the depo in kg
    h2_boil_off_mass: float, mass of hydrogen lost to boil_off in kg
    o2_boil_off_mass: float, mass of oxygen lost to boil_off in kg
    h2_power_mass: float, mass of hydrogen used for power generation in kg
    o2_power_mass: float, mass of oxygen used for power generation in kg
    coolant_mass: float, mass of hydrogen used for reentry cooling in kg
    acs_propellant_mass: float, mass of propellant used for attitude control system in kg

    subsystem_dry_masses: dict[str, float]
        Dictionary of subsystem dry masses in kg, where keys are subsystem names

    sea_level_isp: float in seconds
    vacuum_isp: float in seconds

    landing_delta_v: float in m/s
    deorbit_delta_v: float in m/s
    circularization_delta_v: float in m/s
    orbit_raising_delta_v: float in m/s
    orbit_insertion_delta_v: float in m/s

    of_ratio: float, oxidizer to fuel mass ratio for the main propulsion system
    """

    def calculate_propellant_mass(self) -> None:
        """
        Calculate the total propellant mass required for the mission phases:
        - Landing
        - Deorbit
        - Transfer to depo (circularization and orbit raising)
        - Orbit insertion
        The total mass is calculated by summing the dry mass, propellant masses, boil_off masses,
        payload mass, coolant mass, and power generation masses.
        """
        # start with dry mass
        total_mass = self.dry_mass
        # calculate and add propellant mass for landing
        self.landing_propellant_mass = (np.exp(self.landing_delta_v / (self.sea_level_isp * cn.g_0)) - 1) * total_mass
        total_mass += self.landing_propellant_mass
        # add coolant mass and acs propellant mass
        total_mass += self.coolant_mass + self.acs_propellant_mass + self.h2_power_mass + self.o2_power_mass
        # calculate and add propellant mass for deorbit
        self.deorbit_propellant_mass = (np.exp(self.deorbit_delta_v / (self.vacuum_isp * cn.g_0)) - 1) * total_mass
        total_mass += self.deorbit_propellant_mass

        # add boil_off mass and payload mass
        total_mass += self.h2_boil_off_mass + self.o2_boil_off_mass + self.payload_mass
        # calculate and add propellant mass for circularization and orbit raising
        self.transfer_propellant_mass = (np.exp((self.circularization_delta_v + self.orbit_raising_delta_v)
                                           / (self.vacuum_isp * cn.g_0)) - 1) * total_mass
        total_mass += self.transfer_propellant_mass

        # calculate and add propellant mass for orbit insertion
        self.orbit_insertion_propellant_mass = (np.exp(self.orbit_insertion_delta_v / (self.vacuum_isp * cn.g_0)) - 1) * total_mass
        self.gross_mass = total_mass + self.orbit_insertion_propellant_mass

        self.propellant_mass = (self.landing_propellant_mass
                                + self.deorbit_propellant_mass
                                + self.transfer_propellant_mass
                                + self.orbit_insertion_propellant_mass)


    def calculate_hydrogen_oxygen_mass(self) -> None:
        """
        Calculate the mass of hydrogen and oxygen in the main and header tanks based on OF ratio.
        """
        main_tank_propellant_mass = self.transfer_propellant_mass + self.orbit_insertion_propellant_mass

        self.main_hydrogen_mass = main_tank_propellant_mass * (1 / (1 + self.of_ratio)) + self.h2_boil_off_mass + self.coolant_mass + self.h2_power_mass
        self.main_oxygen_mass = main_tank_propellant_mass * (self.of_ratio / (1 + self.of_ratio)) + self.o2_boil_off_mass + self.o2_power_mass

        self.total_hydrogen_mass = self.main_hydrogen_mass + self.header_hydrogen_mass

        header_tank_propellant_mass = self.landing_propellant_mass + self.deorbit_propellant_mass
        self.header_hydrogen_mass = header_tank_propellant_mass * (1 / (1 + self.of_ratio))  # + self.coolant_mass, potentially add coolant mass here
        self.header_oxygen_mass = header_tank_propellant_mass * (self.of_ratio / (1 + self.of_ratio))

        self.total_oxygen_mass = self.main_oxygen_mass + self.header_oxygen_mass

    def calculate_dummy_dry_masses(self, oi: 'MassIntegrator') -> None:
        """
        Calculate the dry masses of the subsystems based on the masses from previous iteration.
        """
        self.acs_propellant_mass = dds.acs_propellant_mass(oi.gross_mass)

        self.subsystem_dry_masses = {
            "main_tank": dds.main_tank_mass(oi.main_hydrogen_mass + oi.main_oxygen_mass),
            "header_tank": dds.header_tank_mass(oi.header_hydrogen_mass + oi.header_oxygen_mass),
            "landing_legs": dds.landing_leg_mass(oi.dry_mass),
            "acs": dds.acs_dry_mass(self.acs_propellant_mass),
        }

        self.dry_mass = sum(self.subsystem_dry_masses.values())

    def calculate_dry_mass(self, oi: 'MassIntegrator') -> None:
        """
        Calculate the total dry mass of the vehicle by summing the subsystem dry masses.
        """
        self.subsystem_dry_masses = {
            "turbopump": size_turbopump(oi.min_tank_pressure, oi.total_vacuum_thrust),
            'landing_legs': size_landing_legs(
                n_legs=4,
                mass_land=oi.dry_mass,
                phi=oi.phi,
                r_bottom=oi.bottom_radius,
                material=oi.landingleg_material,
                clearance_height=oi.clearance_height
            ),
            "turbopump": size_turbopump(tank_pressure=oi.min_tank_pressure,
                                        thrust=oi.total_vacuum_thrust),
            "main_tank": size_tanks(material=oi.tank_material,
                                    wet_mass=oi.gross_mass,
                                    LH2_mass=oi.total_hydrogen_mass,
                                    LOX_mass=oi.total_oxygen_mass,
                                    LH2_design_pressure=oi.hydrogen_design_pressure,
                                    LOX_design_pressure=oi.oxygen_design_pressure,
                                    thrust_engines=oi.total_vacuum_thrust),

        }

    def calculate_other_masses(self, oi: 'MassIntegrator') -> None:
        """
        Calculate other masses that are not part of the dry mass.
        """
        self.h2_boil_off_mass = dds.h2_boil_off_mass(oi.main_hydrogen_mass)
        self.o2_boil_off_mass = dds.o2_boil_off_mass(oi.main_oxygen_mass)

    def choose_tank_thickness(self, oi:'MassIntegrator') -> None:
        """
        Calculate thickness of the tank from fatigue and tank sizing, and select the higher one.

        """
        # fatigue thickness
        tank_thickness = thickness_optimization_fatigue(oi.cone_angle,
                                                        oi.tank_radius,
                                                        oi.thickness_tank,
                                                        oi.tank_material,
                                                        oi.fuel_reentry_LH2,
                                                        oi.dry_mass,
                                                        oi.launch_mass,
                                                        oi.payload_mass,
                                                        oi.g_reentry_force_ratio,
                                                        oi.g_launch_force_ratio,
                                                        oi.max_thrust2weight)


