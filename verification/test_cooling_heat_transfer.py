import math
from h2ermes_tools.cooling.heat_transfer import HeatShield
from h2ermes_tools.cooling.material import SS310


def test_estimate_spherical_mass_full_sphere():
    hs = HeatShield(wall_thickness=0.01, material=SS310, heat_shield_diameter=2.0)
    # For a full sphere, h = 2*radius
    radius = 1.0
    hs.diameter = 2 * radius
    # Use a dummy density for test
    hs.heat_shield_density = 1000.0
    hs.wall_thickness = 0.01
    mass = hs.estimate_heat_shield_spherical_mass()
    math.isclose(mass, 125.6637061436)


def test_estimate_heat_shield_spherical_mass():
    hs = HeatShield(wall_thickness=0.01, material=SS310, heat_shield_diameter=2.0)
    hs.heat_shield_density = 1000.0
    hs.wall_thickness = 0.01
    mass = hs.estimate_heat_shield_spherical_mass()
    assert mass > 0
