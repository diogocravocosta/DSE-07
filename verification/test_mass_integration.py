import numpy as np

from h2ermes_tools.integration.mass_integration import MassIntegrator


def make_mass_integrator():
    """Make a test MassIntegrator object."""
    integrator = MassIntegrator.__new__(MassIntegrator)
    integrator.dry_mass = 30_000.0

    integrator.payload_mass = 10_000.0
    integrator.h2_boil_off_mass = 1000
    integrator.o2_boil_off_mass = 1000
    integrator.h2_power_mass = 225
    integrator.o2_power_mass = 225
    integrator.coolant_mass = 3000.
    integrator.acs_propellant_mass = 1000.0

    integrator.sea_level_isp = 360
    integrator.vacuum_isp = 450

    integrator.landing_delta_v = 525
    integrator.deorbit_delta_v = 180
    integrator.circularization_delta_v = 150
    integrator.orbit_raising_delta_v = 150
    integrator.orbit_insertion_delta_v = 6000

    integrator.of_ratio = 6.0

    return integrator

def test_calculate_propellant_mass():
    """Test the propellant mass calculation."""
    integrator = make_mass_integrator()
    integrator.calculate_propellant_mass()

    expected_landing_mass = 4810.0451
    expected_deorbit_mass = 1634.4715
    expected_transfer_mass = 3720.8679
    expected_insertion_mass = 163_885.6476

    expected_propellant_mass = 174051.0321

    expected_gross_mass = 220501.032

    assert np.isclose(integrator.landing_propellant_mass, expected_landing_mass, rtol=1e-4)
    assert np.isclose(integrator.deorbit_propellant_mass, expected_deorbit_mass, rtol=1e-4)
    assert np.isclose(integrator.transfer_propellant_mass, expected_transfer_mass, rtol=1e-4)
    assert np.isclose(integrator.orbit_insertion_propellant_mass, expected_insertion_mass, rtol=1e-4)
    assert np.isclose(integrator.propellant_mass, expected_propellant_mass, rtol=1e-4)
    assert np.isclose(integrator.gross_mass, expected_gross_mass, rtol=1e-4)

def test_calculate_hydrogen_oxygen_mass():
    integrator = make_mass_integrator()

    integrator.transfer_propellant_mass = 3720.8679
    integrator.orbit_insertion_propellant_mass = 163885.6476
    integrator.deorbit_propellant_mass = 1634.4715
    integrator.landing_propellant_mass = 4810.0451

    expected_main_hydrogen_mass = 28168.7878
    expected_main_oxygen_mass = 144887.7275

    expected_header_hydrogen_mass = 920.6452
    expected_header_oxygen_mass = 5523.8714

    integrator.calculate_hydrogen_oxygen_mass()

    assert np.isclose(integrator.main_hydrogen_mass, expected_main_hydrogen_mass, rtol=1e-4)
    assert np.isclose(integrator.main_oxygen_mass, expected_main_oxygen_mass, rtol=1e-4)
    assert np.isclose(integrator.header_hydrogen_mass, expected_header_hydrogen_mass, rtol=1e-4)
    assert np.isclose(integrator.header_oxygen_mass, expected_header_oxygen_mass, rtol=1e-4)
