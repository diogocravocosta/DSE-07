import numpy as np
import pytest

from h2ermes_tools.delta_v.landing_burn import landing_burn_no_drag

def test_landing_burn_no_drag():
    initial_velocity = 50.0  # m/s
    thrust_to_weight_ratio = 1.2  # dimensionless

    landing_burn_distance, delta_v = landing_burn_no_drag(initial_velocity, thrust_to_weight_ratio)

    # Expected values based on the formula
    expected_delta_v = 300 # m/s, manually calculated
    expected_landing_burn_distance = 637.105  # m, manually calculated

    assert np.isclose(delta_v, expected_delta_v, rtol=1e-3)
    assert np.isclose(landing_burn_distance, expected_landing_burn_distance, rtol=1e-3)

def test_landing_burn_no_drag_invalid_ratio():
    initial_velocity = 50.0  # m/s
    thrust_to_weight_ratio = 1.0  # dimensionless, invalid for landing burn

    with pytest.raises(ValueError) as exc_info:
        landing_burn_no_drag(initial_velocity, thrust_to_weight_ratio)

    assert str(exc_info.value) == "Thrust-to-weight ratio must be greater than 1 for a landing burn."