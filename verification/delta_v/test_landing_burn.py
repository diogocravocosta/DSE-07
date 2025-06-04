import numpy as np
import pytest
from h2ermes_tools.delta_v.helpers import delta_v_from_final_mass

from h2ermes_tools.delta_v.landing_burn import landing_burn_no_drag, landing_burn_with_drag
import data.constants as cn


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

def test_landing_burn_with_drag_with_falcon_9_values():
    initial_velocity = 950/3.6 # Convert from km/h to m/s
    initial_altitude = 3400 # m
    thrust = 775.65 / 9 * 0.6 # tonnes, one engine at 70% thrust
    dry_mass = 27.2 # tonnes, Falcon 9 first stage dry mass
    approximate_propellant_mass = 5 # tonnes, approximate propellant mass for landing burn
    mass = dry_mass + approximate_propellant_mass # tonnes
    specific_impulse = 283 # seconds, Falcon 9 first stage sea level specific impulse
    thrust_to_weight_ratio = thrust / mass
    ballistic_coefficient = mass*1000 / ((3.7/2)**2 * np.pi * 0.75) # kg/m^2, assuming a diameter of 3.7 m and a drag coefficient of 0.75

    burn_time = 30 # seconds, typical for Falcon 9 landing burn

    df = landing_burn_with_drag(initial_velocity, thrust_to_weight_ratio, ballistic_coefficient, specific_impulse)
    delta_v = delta_v_from_final_mass(df['mass_fraction'].iloc[-1], specific_impulse)
    propellant_mass = dry_mass * (np.exp(delta_v / (cn.g_0 * specific_impulse)) - 1)

    # Check if the time is reasonable
    assert abs(df['time'].iloc[-1] - burn_time) < 5
    # Check if the propellant mass is reasonable
    assert np.isclose(approximate_propellant_mass, propellant_mass , rtol=0.1)
    # Check if the distance is reasonable
    assert np.isclose(-df['distance'].iloc[-1], initial_altitude, rtol=0.2)

def test_cross_verification():
    initial_velocity = 50.0  # m/s
    thrust_to_weight_ratio = 1.2  # dimensionless
    ballistic_coefficient = float('inf')  # kg/m^2, infinite for no drag
    specific_impulse = float('inf')  # seconds, preserves the constant mass used in landing_burn_no_drag

    # Test without drag
    distance_no_drag, delta_v_no_drag = landing_burn_no_drag(initial_velocity, thrust_to_weight_ratio)

    # Test with drag
    df_with_drag = landing_burn_with_drag(initial_velocity, thrust_to_weight_ratio, ballistic_coefficient, specific_impulse)
    delta_v_with_drag = thrust_to_weight_ratio * cn.g_0 * df_with_drag['time'].iloc[-1]

    # Check that the delta_v with drag is less than without drag
    assert np.isclose(delta_v_with_drag, delta_v_no_drag, rtol=1e-3)
    # Check that the distance with drag is less than without drag
    assert np.isclose(-df_with_drag['distance'].iloc[-1], distance_no_drag, rtol=1e-3)