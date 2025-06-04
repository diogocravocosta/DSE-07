import numpy as np
import matplotlib.pyplot as plt

import data.constants as cn

def landing_burn_no_drag(initial_velocity: float,
                         thrust_to_weight_ratio: float,
                         ):
    """
    Calculate the landing burn distance and delta v without considering drag assuming purely vertical motion.

    Args:
        initial_velocity (float): Initial velocity in m/s.
        thrust_to_weight_ratio (float): Thrust-to-weight ratio (dimensionless).

    Returns:
        tuple: (landing_burn_distance, delta_v) where:
            landing_burn_distance (float): Distance required for landing burn in meters.
            delta_v (float): Change in velocity during the landing burn in m/s.
    """
    if thrust_to_weight_ratio <= 1:
        raise ValueError("Thrust-to-weight ratio must be greater than 1 for a landing burn.")

    time = initial_velocity / ((thrust_to_weight_ratio - 1) * cn.g_0)
    landing_burn_distance = initial_velocity * time - 0.5 * (thrust_to_weight_ratio - 1) * cn.g_0 * time**2
    delta_v = thrust_to_weight_ratio * cn.g_0 * time

    return landing_burn_distance, delta_v