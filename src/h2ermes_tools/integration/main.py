import numpy as np

import h2ermes_tools.variables as vr
import h2ermes_tools.integration.mass_integration as mi

# Add unchanging variables to the MassIntegrator
def add_unchanging_variables(integrator: mi.MassIntegrator) -> None:
    """Add unchanging variables to the MassIntegrator object."""
    integrator.payload_mass = vr.payload_mass.value  # Unchanging
    integrator.h2_power_mass = vr.h2_power_mass.value  # Unchanging
    integrator.o2_power_mass = vr.o2_power_mass.value  # Unchanging
    integrator.acs_propellant_mass = vr.acs_propellant_mass.value  # Unchanging

    integrator.sea_level_isp = vr.sea_level_isp.value  # Unchanging
    integrator.vacuum_isp = 429.7477666666666  # Unchanging, from Diogo

    integrator.landing_delta_v = vr.landing_delta_v.value  # Unchanging
    integrator.deorbit_delta_v = vr.deorbit_delta_v.value  # Unchanging
    integrator.circularization_delta_v = vr.circularization_delta_v.value  # Unchanging
    integrator.orbit_raising_delta_v = vr.orbit_raising_delta_v.value  # Unchanging
    integrator.orbit_insertion_delta_v = vr.orbit_insertion_delta_v.value  # Unchanging

    integrator.of_ratio = vr.of_ratio.value  # Unchanging


# Define initial values for the MassIntegrator
def add_initial_values(integrator: mi.MassIntegrator) -> None:
    """Add initial values to the MassIntegrator object."""
    integrator.dry_mass = vr.total_dry_mass.value # Changing
    integrator.h2_boil_off_mass = 1000 # Changing

    integrator.landing_delta_v = 525 # Changing
