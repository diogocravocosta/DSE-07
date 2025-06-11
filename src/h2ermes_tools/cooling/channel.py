from math import pi

class CircularChannel:
    """
    A class representing a circular channel for fluid flow.

    Attributes:
        diameter (float): Diameter of the channel in meters.
        length (float): Length of the channel in meters.
        roughness (float): Roughness height of the channel in meters.
    """

    def __init__(self, diameter, length, roughness):
        self.diameter = diameter
        self.hydraulic_diameter = diameter
        self.cross_sectional_area = pi * (diameter / 2) ** 2
        self.length = length
        self.roughness = roughness

    def get_hydraulic_diameter(self):
        """
        Get the hydraulic diameter of the circular channel.

        Returns:
            float: Hydraulic diameter in meters.
        """
        return self.hydraulic_diameter


class RectangularChannel:
    """
    A class representing a rectangular channel for fluid flow.

    Attributes:
        width (float): Width of the channel in meters.
        height (float): Height of the channel in meters.
        length (float): Length of the channel in meters.
        roughness (float): Roughness height of the channel in meters.
    """

    def __init__(self, width, height, length, roughness):
        self.width = width
        self.height = height
        self.hydraulic_diameter = self.get_hydraulic_diameter()
        self.cross_sectional_area = width * height
        self.length = length
        self.roughness = roughness

    def get_hydraulic_diameter(self):
        """
        Calculate the hydraulic diameter of the channel.

        Returns:
            float: Hydraulic diameter in meters.
        """
        return 2 * (self.width * self.height) / (self.width + self.height)
    
    def get_contact_area(self, segment_length=1.0):
        """
        Calculate the contact area of the channel with the cold side wall.

        Returns:
            float: Contact area in square meters.
        """
        return self.width * segment_length
