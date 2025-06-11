from dataclasses import dataclass

from h2ermes_tools.variables import landing_delta_v



payload_mass = 10_000.0
sea_level_isp = 360
vacuum_isp = 450

landing_delta_v = 525
deorbit_delta_v = 180
circularization_delta_v = 150
orbit_raising_delta_v = 150
orbit_insertion_delta_v = 6000

initial_dry_mass = 30_000.0

@dataclass
class MissionParameters:
    """
    Class to hold mission parameters for the H2ERMES mission.

    Attributes:
        sea_level_isp (float): Sea level specific impulse in seconds.
        vacuum_isp (float): Vacuum specific impulse in seconds.
        of_ratio (float): Oxidizer to fuel mass ratio for the main propulsion system.

        landing_delta_v (float): Delta V required for landing in m/s.
        deorbit_delta_v (float): Delta V required for deorbit in m/s.
        circularization_delta_v (float): Delta V required for circularization in m/s.
        orbit_raising_delta_v (float): Delta V required for orbit raising in m/s.
        orbit_insertion_delta_v (float): Delta V required for orbit insertion in m/s.
    """
    sea_level_isp: float = sea_level_isp
    vacuum_isp: float = vacuum_isp
    of_ratio: float = 6.0

    landing_delta_v: float = landing_delta_v
    deorbit_delta_v: float = deorbit_delta_v
    circularization_delta_v: float = circularization_delta_v
    orbit_raising_delta_v: float = orbit_raising_delta_v
    orbit_insertion_delta_v: float = orbit_insertion_delta_v

@dataclass
class MassParameters:
    """
    Class to hold mass parameters for the H2ERMES mission.

    Attributes:
        dry_mass (float): Initial dry mass of the vehicle in kg.
        payload_mass (float): Mass of the payload in kg.
    """
    dry_mass: float = initial_dry_mass
    payload_mass: float = payload_mass

