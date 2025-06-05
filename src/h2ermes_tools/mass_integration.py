import numpy as np

import data.constants as cn


class MassIntegrator:
    """
    Class to handle mass integration

    dry_mass: float in kg

    propellant_mass: float, mass of propellant used for primary propulsion in kg
        landing_propellant_mass: float, mass of propellant used for landing in kg
        deorbit_propellant_mass: float, mass of propellant used for deorbit in kg
        transfer_propellant_mass: float, mass of propellant used for transfer to depo in kg
        orbit_insertion_propellant_mass: float, mass of propellant used for orbit insertion in kg

    payload_mass: float, payload mass to be transferred to the depo in kg
    h2_boiloff_mass: float, mass of hydrogen lost to boiloff in kg
    o2_boiloff_mass: float, mass of oxygen lost to boiloff in kg
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
        Calculates the total propellant mass required for the mission.
        And separately the propellant mass left after undocking.
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

        # add boiloff mass and payload mass
        total_mass += self.h2_boiloff_mass + self.o2_boiloff_mass + self.payload_mass
        # calculate and add propellant mass for circularization and orbit raising
        self.transfer_propellant_mass = (np.exp((self.circularization_delta_v + self.orbit_raising_delta_v)
                                           / (self.vacuum_isp * cn.g_0)) - 1) * total_mass
        total_mass += self.transfer_propellant_mass

        # calculate and add propellant mass for orbit insertion
        self.orbit_insertion_propellant_mass = (np.exp(self.orbit_insertion_delta_v / (self.vacuum_isp * cn.g_0)) - 1) * total_mass
        self.total_mass = total_mass + self.orbit_insertion_propellant_mass

        self.propellant_mass = (self.landing_propellant_mass
                                + self.deorbit_propellant_mass
                                + self.transfer_propellant_mass
                                + self.orbit_insertion_propellant_mass)


    def calculate_hydrogen_oxygen_mass(self) -> None:
        self.hydrogen_mass = self.propellant_mass * (1 / (1 + self.of_ratio)) + self.h2_boiloff_mass + self.coolant_mass + self.h2_power_mass
        self.oxygen_mass = self.propellant_mass * (self.of_ratio / (1 + self.of_ratio)) + self.o2_boiloff_mass + self.o2_power_mass

        header_tank_propellant_mass = self.landing_propellant_mass + self.deorbit_propellant_mass
        self.header_hydrogen_mass = header_tank_propellant_mass * (1 / (1 + self.of_ratio))  # + self.coolant_mass, potentially add coolant mass here
        self.header_oxygen_mass = header_tank_propellant_mass * (self.of_ratio / (1 + self.of_ratio))
