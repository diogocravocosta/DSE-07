from h2ermes_tools.cooling.channel import RectangularChannel, CircularChannel
import math

def test_rectangular_channel_hydraulic_diameter():
    channel = RectangularChannel(width=0.1, height=0.05, length=1.0, roughness=1e-5)
    expected = 2 * (0.1 * 0.05) / (0.1 + 0.05)  # known formula for hydraulic diameter
    assert math.isclose(channel.get_hydraulic_diameter(), expected, rel_tol=1e-9)

def test_rectangular_channel_cross_sectional_area():
    channel = RectangularChannel(width=0.1, height=0.05, length=1.0, roughness=1e-5)
    assert math.isclose(channel.cross_sectional_area, 0.1 * 0.05, rel_tol=1e-9)

def test_rectangular_channel_contact_area():
    channel = RectangularChannel(width=0.1, height=0.05, length=1.0, roughness=1e-5)
    segment_length = 2.0
    expected = 0.1 * segment_length
    assert math.isclose(channel.get_contact_area(segment_length), expected, rel_tol=1e-9)

def test_circular_channel_hydraulic_diameter():
    channel = CircularChannel(diameter=0.1, length=1.0, roughness=1e-5)
    assert math.isclose(channel.get_hydraulic_diameter(), 0.1, rel_tol=1e-9)

def test_circular_channel_cross_sectional_area():
    channel = CircularChannel(diameter=0.1, length=1.0, roughness=1e-5)
    expected = math.pi * (0.1 / 2) ** 2
    assert math.isclose(channel.cross_sectional_area, expected, rel_tol=1e-9)
