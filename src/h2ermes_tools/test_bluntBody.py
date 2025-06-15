import pytest
import numpy as np
from aerodynamic_coeffs import BluntBody

def test_intialization():
    cone_length = 10
    cone_max_radius = 5
    cone_min_radius = 2
    base_arc_height = 1
    mass = 1000

    body = BluntBody(cone_length, cone_max_radius, cone_min_radius, base_arc_height, mass)

    # Check attributes
    assert body.cone_length == cone_length
    assert body.mass == mass

    # sphere radius should be >= cone_max_radius
    assert body.sphere_radius >= cone_max_radius

def test_compute_aerodynamics():
    body = BluntBody(10, 5, 2, 1, 1000)
    body.compute_aerodynamics()
     #add all the tests lmao yay
    
def test_stability():
    body = Bluntbody(10,5,2,1,1000)
    body.stability(mode=2)


