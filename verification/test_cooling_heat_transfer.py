import math
from h2ermes_tools.cooling.heat_transfer import HeatShield
from h2ermes_tools.cooling.material import SS310


def test_estimate_spherical_mass_full_sphere():
    hs = HeatShield(wall_thickness=0.01, material=SS310, heat_shield_diameter=2.0)
    # For a full sphere, h = 2*radius
    radius = 1.0
    h = 2 * radius
    hs.diameter = 2 * radius
    # Use a dummy density for test
    hs.heat_shield_density = 1000.0
    hs.wall_thickness = 0.01
    mass = hs.estimate_spherical_mass(h=h, radius=radius)
    math.isclose(mass, 125.6637061436)


def test_estimate_spherical_mass_partial():
    hs = HeatShield(wall_thickness=0.01, material=SS310, heat_shield_diameter=2.0)
    radius = 1.0
    h = 0.5 * radius
    hs.diameter = 2 * radius
    hs.heat_shield_density = 1000.0
    hs.wall_thickness = 0.01
    mass = hs.estimate_spherical_mass(h=h, radius=radius)
    assert mass > 0
    # Should be less than full sphere
    mass_full = hs.estimate_spherical_mass(h=2 * radius, radius=radius)
    assert mass < mass_full


def test_estimate_spherical_mass_uses_diameter():
    hs = HeatShield(wall_thickness=0.01, material=SS310, heat_shield_diameter=2.0)
    hs.heat_shield_density = 1000.0
    hs.wall_thickness = 0.01
    h = 1.0
    # Should not raise if radius is None and diameter is set
    mass = hs.estimate_spherical_mass(h=h)
    assert mass > 0
