import numpy as np

import data.material as dm
import data.constants as cn

import h2ermes_tools.variables as vr
import h2ermes_tools.integration.mass_integration as mi
from h2ermes_tools.integration.main import run_integration

def add_sensitivity_analysis_initial_values(integrator: mi.MassIntegrator) -> None:
    """Add initial values to the MassIntegrator object."""
    multiplier = 4

    integrator.dry_mass = vr.total_dry_mass.value * multiplier
    integrator.h2_boil_off_mass = 1000 * multiplier
    integrator.gross_mass = vr.gross_mass.value * multiplier

    integrator.total_vacuum_thrust = integrator.gross_mass * cn.g_0 * integrator.vacuum_twr * multiplier

    integrator.header_hydrogen_mass = 1000 * multiplier
    integrator.header_oxygen_mass = 5000 * multiplier

    integrator.phi_top = np.deg2rad(10)  # rad
    integrator.middle_radius = 4.8  # m
    integrator.hydrogen_tank_height = 12.65  # m
    integrator.top_radius = 2.5  # m

if __name__ == '__main__':
    run_integration(add_sensitivity_analysis_initial_values)