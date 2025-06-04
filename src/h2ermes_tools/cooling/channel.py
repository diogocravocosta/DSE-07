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
        self.length = length
        self.roughness = roughness

    def hydraulic_diameter(self):
        """
        Calculate the hydraulic diameter of the channel.

        Returns:
            float: Hydraulic diameter in meters.
        """
        return self.diameter
