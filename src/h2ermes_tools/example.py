"""This file contains a small code example for everyone to draw from"""

import numpy as np

import data.constants as cn

def free_fall_velocity(time: float) -> float:
    """
    Calculates free fall velocity after given time
    Args:
        time:

    Returns:
    velocity
    """
    return cn.g_0*time

def dot_product(vector1: list[float], vector2: list[float]) -> float:
    """
    Returns sum of element-wise multiplication of two vectors of the same length.

    Args:
        vector1: first vector
        vector2: second vector

    Returns:
    sum of element-wise multiplication of vector1 and vector2
    """
    if len(vector1) != len(vector2):
        raise ValueError("vector1 and vector2 must be the same length")

    sum = 0.
    for i in range(len(vector1)):
        sum += vector1[i] * vector2[i]
    return sum

class Cylinder:
    def __init__(self, radius: float, height: float) -> None:
        self.radius = radius # m
        self.height = height # m
        self.area = None # m2
        self.volume = None # m3

    def calculate_base_area(self):
        """
        Calculates the area of the cylinder base
        Returns:
        area of base
        """
        self.area = np.pi * self.radius**2

        # For people who know properties and caching and want to use it: Go for it.
        # I just didn't want to confuse people who don't more than necessary

    def calculate_volume(self):
        """
        Calculates the volume of the cylinder
        Returns:
        volume of cylinder
        """
        if self.area is None:
            self.calculate_base_area()
        self.volume =  self.height * self.area

    def print_area_volume(self):
        """Prints the area and volume of the cylinder"""
        print("area:" + str(self.area))
        print("volume:" + str(self.volume))