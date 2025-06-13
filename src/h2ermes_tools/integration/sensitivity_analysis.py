import numpy as np

import data.material as dm
import data.constants as cn

import h2ermes_tools.variables as vr
import h2ermes_tools.integration.mass_integration as mi
from h2ermes_tools.integration.main import run_integration


# Add unchanging variables to the MassIntegrator
def add_sensitivity_unchanging_variables(integrator: mi.MassIntegrator) -> None:
    """Add unchanging variables to the MassIntegrator object."""
    integrator.payload_mass = 8_000
    integrator.h2_power_mass = vr.hydrogen_power_mass.value
    integrator.o2_power_mass = vr.oxygen_power_mass.value
    integrator.acs_propellant_mass = vr.acs_propellant_mass.value
    integrator.coolant_mass = vr.coolant_mass.value

    integrator.sea_level_isp = vr.engine_specific_impulse_sl_opt.value
    integrator.vacuum_isp = 429.7477666666666  # from Diogo
    integrator.sea_level_thrust_chambers_number = vr.n_chambers_sl.value
    integrator.vacuum_thrust_chambers_number = vr.n_chambers_vac.value

    integrator.landing_delta_v = vr.landing_delta_v.value
    integrator.deorbit_delta_v = vr.deorbit_delta_v.value + 50
    integrator.circularization_delta_v = vr.target_orbit_circularization_delta_v.value
    integrator.orbit_raising_delta_v = vr.orbit_raising_delta_v.value
    integrator.orbit_insertion_delta_v = vr.orbit_insertion_delta_v.value

    integrator.of_ratio = vr.engine_mixture_ratio.value

    integrator.vacuum_twr = vr.t_w_vac.value

    integrator.clearance_height = 2.5  # m

    integrator.hydrogen_design_pressure = 2e5
    integrator.oxygen_design_pressure = 2.5e5
    integrator.boiloff_design_pressure = 10e5

    integrator.landing_leg_material = dm.Ti6Al4V
    integrator.tank_material = dm.random_steel
    integrator.header_tank_material = dm.random_steel

    integrator.g_reentry_force_ratio = 8
    integrator.g_launch_force_ratio = 6

    # Set parameters which will probably change in actuality
    integrator.bottom_radius = 5  # m

    # Add unchanging subsystem dry masses
    integrator.unchanging_subsystem_dry_masses = {
        "acs"            : vr.acs_dry_mass.value,
        "heat shield"    : vr.heat_shield_mass.value,
        "thrust chambers": vr.thrust_chamber_sl_mass.value * integrator.sea_level_thrust_chambers_number
                           + vr.thrust_chamber_vac_mass.value * integrator.vacuum_thrust_chambers_number,
        "docking system" : vr.docking_system_mass.value,
        "nose cone"      : vr.nose_cone_mass.value,
        "power"          : vr.power_dry_mass.value,
        "avionics"       : vr.avionics_mass.value,
        "wiring"         : vr.wiring_mass.value,
        "interstage"     : vr.interstage_mass.value,
    }

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
    run_integration(add_sensitivity_analysis_initial_values,
                    add_sensitivity_unchanging_variables)