import pytest
import numpy as np

import h2ermes_tools.delta_v.helpers as helpers

def test_energy():
    altitude = 2e5
    expected_energy = -30297373.76 # manually calculated

    assert helpers.calculate_circular_orbit_energy(altitude) == pytest.approx(expected_energy, rel=1e-3)

def test_delta_V():
    specific_impulse = 450
    mass_fraction = 0.5

    expected_delta_v = 3059.898229

    assert helpers.delta_v_from_final_mass(mass_fraction, specific_impulse) == pytest.approx(expected_delta_v, rel=1e-3)
