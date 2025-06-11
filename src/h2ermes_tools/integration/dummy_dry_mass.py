"""This file contains a dummy dry mass function for testing purposes."""

def main_tank_mass(main_tank_propellant_mass: float) -> float:
    """
    Calculate the dry mass of the main tank based on the propellant mass.

    Args:
        main_tank_propellant_mass (float): The mass of the propellant in the main tank in kg.

    Returns:
        float: The dry mass of the main tank in kg.
    """
    # Dummy calculation for testing purposes
    return 0.1 * main_tank_propellant_mass  # 10% of the propellant mass as dry mass

def header_tank_mass(header_tank_propellant_mass: float) -> float:
    """
    Calculate the dry mass of the header tank based on the propellant mass.

    Args:
        header_tank_propellant_mass (float): The mass of the propellant in the header tank in kg.

    Returns:
        float: The dry mass of the header tank in kg.
    """
    # Dummy calculation for testing purposes
    return 0.2 * header_tank_propellant_mass  # 20% of the propellant mass as dry mass

def landing_leg_mass(dry_mass: float) -> float:
    """
    Calculate the dry mass of the landing legs based on the total dry mass.

    Args:
        dry_mass (float): The total dry mass of the vehicle in kg.

    Returns:
        float: The dry mass of the landing legs in kg.
    """
    # Dummy calculation for testing purposes
    return 0.05 * dry_mass  # 5% of the total dry mass as landing leg mass

def acs_dry_mass(acs_propellant_mass: float) -> float:
    """
    Calculate the dry mass of the ACS based on the propellant mass.

    Args:
        acs_propellant_mass (float): The mass of the propellant in the ACS in kg.

    Returns:
        float: The dry mass of the ACS in kg.
    """
    # Dummy calculation for testing purposes
    return 1 * acs_propellant_mass  # 15% of the propellant mass as dry mass

def acs_propellant_mass(gross_mass: float) -> float:
    """
    Calculate the propellant mass of the ACS based on the dry mass.

    Args:
        gross_mass (float):

    Returns:
        float: The propellant mass of the ACS in kg.
    """
    # Dummy calculation for testing purposes
    return 0.005 * gross_mass  # Assuming 15% dry mass to propellant mass ratio