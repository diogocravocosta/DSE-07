"""WIP"""

import numpy.testing as npt

from h2ermes_tools.cooling.channel import CircularChannel

def test_diameter():
    channel = CircularChannel(diameter=0.1, length=1.0, roughness=1e-5)
    expected_diameter = 0.1

    npt.assert_almost_equal(channel.diameter, expected_diameter)
