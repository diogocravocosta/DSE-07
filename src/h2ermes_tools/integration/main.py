import numpy as np

import h2ermes_tools.variables as vr
import h2ermes_tools.integration.mass_integration as mi
import data.material as dm

# Add unchanging variables to the MassIntegrator
def add_unchanging_variables(integrator: mi.MassIntegrator) -> None:
    """Add unchanging variables to the MassIntegrator object."""
    integrator.payload_mass = vr.payload_mass.value  # Unchanging
    integrator.h2_power_mass = vr.hydrogen_power_mass.value  # Unchanging
    integrator.o2_power_mass = vr.oxygen_power_mass.value  # Unchanging
    integrator.acs_propellant_mass = vr.acs_propellant_mass.value  # Unchanging
    integrator.coolant_mass = vr.coolant_mass.value  # Unchanging

    integrator.sea_level_isp = vr.engine_specific_impulse_sl_opt.value  # Unchanging
    integrator.vacuum_isp = 429.7477666666666  # Unchanging, from Diogo

    integrator.landing_delta_v = vr.landing_delta_v.value  # Unchanging
    integrator.deorbit_delta_v = vr.deorbit_delta_v.value  # Unchanging
    integrator.circularization_delta_v = vr.target_orbit_circularization_delta_v.value  # Unchanging
    integrator.orbit_raising_delta_v = vr.orbit_raising_delta_v.value  # Unchanging
    integrator.orbit_insertion_delta_v = vr.orbit_insertion_delta_v.value  # Unchanging

    integrator.of_ratio = vr.engine_mixture_ratio.value  # Unchanging

    # Set parameters which will probably change
    integrator.phi = np.deg2rad(10) # rad
    integrator.bottom_radius = 5 # m
    integrator.landing_leg_material = dm.SS14845
    integrator.clearance_height = 2.5 # m
    integrator.min_tank_pressure = 2e5 # Pa

# Define initial values for the MassIntegrator
def add_initial_values(integrator: mi.MassIntegrator) -> None:
    """Add initial values to the MassIntegrator object."""
    integrator.dry_mass = vr.total_dry_mass.value
    integrator.h2_boil_off_mass = 1000

def main() -> None:
    """Main function to run the mass integration."""
    # Create a MassIntegrator instance
    old_integrator = mi.MassIntegrator()

    # Add unchanging variables to the integrator
    add_unchanging_variables(old_integrator)

    # Add initial values to the integrator
    add_initial_values(old_integrator)

    # Perform the mass integration
    old_integrator.calculate_propellant_mass()
    old_integrator.calculate_hydrogen_oxygen_mass()

    for i in range(10):
        new_integrator = mi.MassIntegrator()
        add_unchanging_variables(new_integrator)
        new_integrator.calculate_dry_mass(old_integrator)

    # Print the results
    print(f"Total mass: {old_integrator.gross_mass:.2f} kg")
    print(f"Dry mass: {old_integrator.dry_mass:.2f} kg")
    print(f"Propellant mass: {old_integrator.propellant_mass:.2f} kg")
    print(f"Hydrogen mass: {old_integrator.total_hydrogen_mass:.2f} kg")
    print(f"Oxygen mass: {old_integrator.total_oxygen_mass:.2f} kg")

if __name__ == "__main__":
    main()