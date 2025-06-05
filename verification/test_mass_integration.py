import numpy as np

from h2ermes_tools.mass_integration import MassIntegrator


def make_mass_integrator():
    """Make a test MassIntegrator object."""
    integrator = MassIntegrator.__new__(MassIntegrator)
    integrator.dry_mass = 30_000.0

    integrator.payload_mass = 10_000.0
    integrator.h2_boiloff_mass = 1000
    integrator.o2_boiloff_mass = 1000
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

    expected_total_mass = 220501.032

    assert np.isclose(integrator.landing_propellant_mass, expected_landing_mass, rtol=1e-4)
    assert np.isclose(integrator.deorbit_propellant_mass, expected_deorbit_mass, rtol=1e-4)
    assert np.isclose(integrator.transfer_propellant_mass, expected_transfer_mass, rtol=1e-4)
    assert np.isclose(integrator.orbit_insertion_propellant_mass, expected_insertion_mass, rtol=1e-4)
    assert np.isclose(integrator.total_mass, expected_total_mass, rtol=1e-4)
