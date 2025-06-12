import numpy as np

import data.material as mat

def calculate_tank_mass(bottom_radius: float,
                        top_radius: float,
                        tank_length: float,
                        thickness: float,
                        material: mat.Material,
                        cap_height_radius_ratio: float = 0.5) -> float:
    volume_shell_wall = (
            (np.pi * tank_length / 3)
            * (bottom_radius ** 2 + bottom_radius * top_radius + top_radius ** 2
               - (bottom_radius - thickness) ** 2 - (bottom_radius - thickness) * (top_radius - thickness) - (
                       top_radius - thickness) ** 2)
            + 2 * np.pi / 3 * ((bottom_radius + thickness) ** 3 - bottom_radius ** 3) * cap_height_radius_ratio
            + 2 * np.pi / 3 * ((top_radius + thickness) ** 3 - top_radius ** 3) * cap_height_radius_ratio
    )
    mass = volume_shell_wall * material.rho
    return mass