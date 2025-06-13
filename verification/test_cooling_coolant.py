from pyfluids import Fluid, FluidsList, Input
from h2ermes_tools.cooling.coolant import Coolant
from h2ermes_tools.cooling.channel import RectangularChannel, CircularChannel

def make_hydrogen_fluid(T=20.0, P=1e6):
    return Fluid(FluidsList.Hydrogen).with_state(Input.temperature(T), Input.pressure(P))

def test_get_reynolds_number():
    fluid = make_hydrogen_fluid()
    channel = RectangularChannel(width=0.01, height=0.01, length=1.0, roughness=1e-5)
    coolant = Coolant(fluid=fluid, channel=channel, mass_flow=0.01)
    Re = coolant.get_reynolds_number()
    assert Re > 0

def test_get_fluid_speed():
    fluid = make_hydrogen_fluid()
    channel = RectangularChannel(width=0.01, height=0.01, length=1.0, roughness=1e-5)
    coolant = Coolant(fluid=fluid, channel=channel, mass_flow=0.01)
    speed = coolant.get_fluid_speed()
    assert speed > 0

def test_get_heat_transfer_coefficient():
    fluid = make_hydrogen_fluid()
    channel = RectangularChannel(width=0.01, height=0.01, length=1.0, roughness=1e-5)
    coolant = Coolant(fluid=fluid, channel=channel, mass_flow=0.01)
    htc = coolant.get_heat_transfer_coefficient()
    assert htc > 0

def test_add_energy_increases_temperature():
    fluid = make_hydrogen_fluid(T=20.0)
    channel = RectangularChannel(width=0.01, height=0.01, length=1.0, roughness=1e-5)
    coolant = Coolant(fluid=fluid, channel=channel, mass_flow=0.01)
    T0 = coolant.fluid.temperature
    coolant.add_energy(1000.0)
    T1 = coolant.fluid.temperature
    assert T1 > T0

def test_calculate_pressure_drop():
    fluid = make_hydrogen_fluid()
    channel = RectangularChannel(width=0.01, height=0.01, length=1.0, roughness=1e-5)
    coolant = Coolant(fluid=fluid, channel=channel, mass_flow=0.01)
    dP = coolant.calculate_pressure_drop(segment_length=0.5)
    assert dP >= 0

def test_circular_channel_works():
    fluid = make_hydrogen_fluid()
    channel = CircularChannel(diameter=0.01, length=1.0, roughness=1e-5)
    coolant = Coolant(fluid=fluid, channel=channel, mass_flow=0.01)
    assert coolant.get_reynolds_number() > 0
    assert coolant.get_fluid_speed() > 0
    assert coolant.get_heat_transfer_coefficient() > 0

def test_rectangular_channel_works():
    fluid = make_hydrogen_fluid()
    channel = RectangularChannel(width=0.01, height=0.01, length=1.0, roughness=1e-5)
    coolant = Coolant(fluid=fluid, channel=channel, mass_flow=0.01)
    assert coolant.get_reynolds_number() > 0
    assert coolant.get_fluid_speed() > 0
    assert coolant.get_heat_transfer_coefficient() > 0
