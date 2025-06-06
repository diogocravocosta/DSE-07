# 3 tanks of oxygen and 3 tanks of hydrogen are good for 8 days, scaling with power output

def lb_to_kg(lb: float) -> float:
    """
    Convert pounds to kilograms.

    Args:
            lb: mass in pounds
    Returns:
            value in kilograms
    """
    return lb * 0.45359237


def power_plant_mass(average_power: float, mission_length_days: float, safety_factor: float) -> float:
    """
    From the mission length and average power usage, calculates total reactant consumption of the fuel cells.
    Adds this to the mass of the 2 fuel cells for the combined mass power-producing equipment.
    Does not include power distribution and regulating equipment mass.
    Based on consumption by the Space Shuttle.

    Args:
            average_power: Power in W consumed on average during the mission
            mission_length_days: Total nominal length of the mission in days
            It is assumed that the average power level applies from the beginning until the end.
            safety_factor: Reactant mass multiplier to extend capabilities beyond nominal mission length

    Returns:
            Combined mass of the fuel cells and reactants
    """

    fuel_cell_mass = lb_to_kg(255)
    number_of_fuel_cells = 2

    oxygen_per_tank   = lb_to_kg(781)
    hydrogen_per_tank = lb_to_kg(92)

    # Shuttle design power output is 14 kW
    # Shuttle needs 3 tanks of oxygen and hydrogen each for an 8-day mission
    shuttle_power = 14000  # W
    shuttle_daily_oxygen   = 3 / 8 * oxygen_per_tank
    shuttle_daily_hydrogen = 3 / 8 * hydrogen_per_tank

    daily_oxygen   = average_power / shuttle_power * shuttle_daily_oxygen
    daily_hydrogen = average_power / shuttle_power * shuttle_daily_hydrogen

    total_oxygen   = daily_oxygen   * mission_length_days * safety_factor
    total_hydrogen = daily_hydrogen * mission_length_days * safety_factor
    #print(total_hydrogen + total_oxygen)
    #print(fuel_cell_mass * number_of_fuel_cells)
    return total_oxygen + total_hydrogen + fuel_cell_mass * number_of_fuel_cells


if __name__ == "__main__":
    our_power = 3573 # W

    print(power_plant_mass(average_power=our_power, mission_length_days=2, safety_factor=2))