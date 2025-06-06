import pytest
import numpy as np

from h2ermes_tools.power_mass import power_plant_mass, lb_to_kg


def test_power():
    shuttle_power = 14000 # kW
    # Nominal Shuttle consumption is 8 tanks of oxygen and hydrogen over 8 days
    calculated_mass = power_plant_mass(average_power=shuttle_power, mission_length_days=8, safety_factor=1)
    # Mass of 2 power cells and contents of 3 oxygen and hydrogen tanks
    ground_truth_mass = 2 * lb_to_kg(255) + 3 * lb_to_kg(781 + 92)

    assert np.isclose(calculated_mass, ground_truth_mass, rtol=1e-3)