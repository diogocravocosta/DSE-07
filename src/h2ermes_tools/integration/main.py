import numpy as np

import data.material as dm
import data.constants as cn

import h2ermes_tools.variables as vr
import h2ermes_tools.integration.mass_integration as mi


# Add unchanging variables to the MassIntegrator
def add_unchanging_variables(integrator: mi.MassIntegrator) -> None:
    """Add unchanging variables to the MassIntegrator object."""
    integrator.payload_mass = 10_000
    integrator.h2_power_mass = vr.hydrogen_power_mass.value  
    integrator.o2_power_mass = vr.oxygen_power_mass.value  
    integrator.acs_propellant_mass = vr.acs_propellant_mass.value  
    integrator.coolant_mass = vr.coolant_mass.value  

    integrator.sea_level_isp = vr.engine_specific_impulse_sl_opt.value  
    integrator.vacuum_isp = 429.7477666666666  # from Diogo

    integrator.landing_delta_v = vr.landing_delta_v.value  
    integrator.deorbit_delta_v = vr.deorbit_delta_v.value + 50
    integrator.circularization_delta_v = vr.target_orbit_circularization_delta_v.value  
    integrator.orbit_raising_delta_v = vr.orbit_raising_delta_v.value  
    integrator.orbit_insertion_delta_v = vr.orbit_insertion_delta_v.value  

    integrator.of_ratio = vr.engine_mixture_ratio.value  

    integrator.vacuum_twr = vr.t_w_vac.value
    
    integrator.clearance_height = 2.5 # m
    integrator.min_tank_pressure = 2e5 # Pa

    integrator.hydrogen_design_pressure = 2e5
    integrator.oxygen_design_pressure = 2.5e5
    integrator.boiloff_design_pressure = 10e5

    integrator.landing_leg_material = dm.Ti6Al4V
    integrator.tank_material = dm.random_steel
    integrator.header_tank_material = dm.random_steel

    integrator.g_reentry_force_ratio = 8
    integrator.g_launch_force_ratio = 6

    # Set parameters which will probably change in actuality
    integrator.bottom_radius = 5 # m


    # Add unchanging subsystem dry masses
    integrator.unchanging_subsystem_dry_masses = {
        "acs": vr.acs_dry_mass.value,
        "heat shield": vr.heat_shield_mass.value,
        "thrust chambers": vr.thrust_chamber_sl_mass.value * vr.n_chambers_sl.value + vr.thrust_chamber_vac_mass.value * vr.n_chambers_vac.value,
        "docking system": vr.docking_system_mass.value,
        "nose cone": vr.nose_cone_mass.value,
        "power": vr.power_dry_mass.value,
        "data handling": vr.data_handling_mass.value,
        "avionics harness": vr.avionics_harness_mass.value,
        "interstage": vr.interstage_mass.value,
        "gnc": vr.interstage_mass.value
    }


# Define initial values for the MassIntegrator
def add_initial_values(integrator: mi.MassIntegrator) -> None:
    """Add initial values to the MassIntegrator object."""
    integrator.dry_mass = vr.total_dry_mass.value
    integrator.h2_boil_off_mass = 1000
    integrator.gross_mass = vr.gross_mass.value

    integrator.total_vacuum_thrust = integrator.gross_mass * cn.g_0 * integrator.vacuum_twr

    integrator.header_hydrogen_mass = 1000
    integrator.header_oxygen_mass = 5000

    integrator.phi_top = np.deg2rad(10) # rad
    integrator.middle_radius = 4.8 # m
    integrator.hydrogen_tank_height = 12.65 # m
    integrator.top_radius = 2.5 # m


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

    for i in range(20):
        new_integrator = mi.MassIntegrator()
        add_unchanging_variables(new_integrator)

        new_integrator.choose_tank_thickness_dimensions(old_integrator)
        new_integrator.calculate_tank_mass()

        new_integrator.calculate_dry_mass(old_integrator)
        new_integrator.calculate_non_propellant_consumables(old_integrator)

        new_integrator.calculate_propellant_mass()
        new_integrator.calculate_hydrogen_oxygen_mass()

        new_integrator.update_thrust()

        # Print the results
        print(f"Hydrogen tank thickness: {new_integrator.lh2_thickness}")
        print(f"Oxygen tank thickness: {new_integrator.lox_thickness}")
        print("--------------------------------------------------------------")
        print(f"Total mass: {new_integrator.gross_mass:.2f} kg")
        print(f"Dry mass: {new_integrator.dry_mass:.2f} kg")
        print(f"Hydrogen mass: {new_integrator.total_hydrogen_mass:.2f} kg")
        print(f"Oxygen mass: {new_integrator.total_oxygen_mass:.2f} kg")
        print(f"ACS propellant mass: {new_integrator.acs_propellant_mass} kg")
        print("--------------------------------------------------------------")

        old_integrator = new_integrator

if __name__ == "__main__":
    main()