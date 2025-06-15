import numpy as np

import data.material as dm
import data.constants as cn

import h2ermes_tools.variables as vr
import h2ermes_tools.integration.mass_integration as mi
from h2ermes_tools.integration.main import run_integration


# Add unchanging variables to the MassIntegrator
def add_sensitivity_unchanging_variables(
    integrator: mi.MassIntegrator,
    payload_mass=8_000,
    h2_power_mass=vr.hydrogen_power_mass.value,
    o2_power_mass=vr.oxygen_power_mass.value,
    acs_propellant_mass=vr.acs_propellant_mass.value,
    coolant_mass=vr.coolant_mass.value,
    sea_level_isp=vr.engine_specific_impulse_sl_opt.value,
    vacuum_isp=429.7477666666666,  # from Diogo
    sea_level_thrust_chambers_number=vr.n_chambers_sl.value,
    vacuum_thrust_chambers_number=vr.n_chambers_vac.value,
    landing_delta_v=vr.landing_delta_v.value,
    deorbit_delta_v=vr.deorbit_delta_v.value + 50,
    circularization_delta_v=vr.target_orbit_circularization_delta_v.value,
    orbit_raising_delta_v=vr.orbit_raising_delta_v.value,
    orbit_insertion_delta_v=vr.orbit_insertion_delta_v.value,
    of_ratio=vr.engine_mixture_ratio.value,
    vacuum_twr=vr.t_w_vac.value,
    clearance_height=2.5,  # m
    hydrogen_design_pressure=2e5,
    oxygen_design_pressure=2.5e5,
    boiloff_design_pressure=10e5,
    landing_leg_material=dm.Ti6Al4V,
    tank_material=dm.random_steel,
    header_tank_material=dm.random_steel,
    g_reentry_force_ratio=8,
    g_launch_force_ratio=6,
    bottom_radius=5,  # m
    unchanging_subsystem_dry_masses=None,
) -> None:
    """Add unchanging variables to the MassIntegrator object."""
    integrator.payload_mass = payload_mass
    integrator.h2_power_mass = h2_power_mass
    integrator.o2_power_mass = o2_power_mass
    integrator.acs_propellant_mass = acs_propellant_mass
    integrator.coolant_mass = coolant_mass

    integrator.sea_level_isp = sea_level_isp
    integrator.vacuum_isp = vacuum_isp
    integrator.sea_level_thrust_chambers_number = sea_level_thrust_chambers_number
    integrator.vacuum_thrust_chambers_number = vacuum_thrust_chambers_number

    integrator.landing_delta_v = landing_delta_v
    integrator.deorbit_delta_v = deorbit_delta_v
    integrator.circularization_delta_v = circularization_delta_v
    integrator.orbit_raising_delta_v = orbit_raising_delta_v
    integrator.orbit_insertion_delta_v = orbit_insertion_delta_v

    integrator.of_ratio = of_ratio

    integrator.vacuum_twr = vacuum_twr

    integrator.clearance_height = clearance_height

    integrator.hydrogen_design_pressure = hydrogen_design_pressure
    integrator.oxygen_design_pressure = oxygen_design_pressure
    integrator.boiloff_design_pressure = boiloff_design_pressure

    integrator.landing_leg_material = landing_leg_material
    integrator.tank_material = tank_material
    integrator.header_tank_material = header_tank_material

    integrator.g_reentry_force_ratio = g_reentry_force_ratio
    integrator.g_launch_force_ratio = g_launch_force_ratio

    integrator.bottom_radius = bottom_radius

    if unchanging_subsystem_dry_masses is None:
        integrator.unchanging_subsystem_dry_masses = {
            "acs": vr.acs_dry_mass.value,
            "heat shield": vr.heat_shield_mass.value,
            "thrust chambers": vr.thrust_chamber_sl_mass.value
            * sea_level_thrust_chambers_number
            + vr.thrust_chamber_vac_mass.value * vacuum_thrust_chambers_number,
            "docking system": vr.docking_system_mass.value,
            "nose cone": vr.nose_cone_mass.value,
            "power": vr.power_dry_mass.value,
            "avionics": vr.avionics_mass.value,
            "wiring": vr.wiring_mass.value,
            "interstage": vr.interstage_mass.value,
        }
    else:
        integrator.unchanging_subsystem_dry_masses = unchanging_subsystem_dry_masses


def add_sensitivity_analysis_initial_values(integrator: mi.MassIntegrator) -> None:
    """Add initial values to the MassIntegrator object."""
    multiplier = 1.0

    integrator.dry_mass = vr.total_dry_mass.value * multiplier
    integrator.h2_boil_off_mass = 2000 * multiplier
    integrator.gross_mass = vr.gross_mass.value * multiplier

    integrator.total_vacuum_thrust = (
        integrator.gross_mass * cn.g_0 * integrator.vacuum_twr * multiplier
    )

    integrator.header_hydrogen_mass = 1000 * multiplier
    integrator.header_oxygen_mass = 5000 * multiplier

    integrator.phi_top = np.deg2rad(10)  # rad
    integrator.middle_radius = 4.8  # m
    integrator.hydrogen_tank_height = 12.65  # m
    integrator.top_radius = 2.5  # m


if __name__ == "__main__":
    h2_power_mass_values = (
        np.linspace(
            min(vr.hydrogen_power_mass.margin), max(vr.hydrogen_power_mass.margin), 3
        )
        * vr.hydrogen_power_mass.value
    )

    acs_propellant_mass_values = (
        np.linspace(
            min(vr.acs_propellant_mass.margin), max(vr.acs_propellant_mass.margin), 3
        )
        * vr.acs_propellant_mass.value
    )

    sl_isp_values = (
        np.linspace(
            min(vr.engine_specific_impulse_sl_opt.margin),
            max(vr.engine_specific_impulse_sl_opt.margin),
            3,
        )
        * vr.engine_specific_impulse_sl_opt.value
    )
    vac_isp_values = (
        np.linspace(
            min(vr.engine_specific_impulse_vac_opt.margin),
            max(vr.engine_specific_impulse_vac_opt.margin),
            5,
        )
        * vr.engine_specific_impulse_vac_opt.value
    )

    landing_delta_v_values = (
        np.linspace(min(vr.landing_delta_v.margin), max(vr.landing_delta_v.margin), 5)
        * vr.landing_delta_v.value
    )

    # Example usage of the first value in the range for demonstration
    gross_masses = []
    dry_masses = []
    total_hydrogen_masses = []
    total_oxygen_masses = []
    acs_propellant_masses = []

    for ldv in landing_delta_v_values:
        (
            gross_mass,
            dry_mass,
            total_hydrogen_mass,
            total_oxygen_mass,
            acs_propellant_mass,
        ) = run_integration(
            add_sensitivity_analysis_initial_values,
            lambda integrator: add_sensitivity_unchanging_variables(
                integrator, landing_delta_v=ldv
            ),
        )
        gross_masses.append(gross_mass)
        dry_masses.append(dry_mass)
        total_hydrogen_masses.append(total_hydrogen_mass)
        total_oxygen_masses.append(total_oxygen_mass)
        acs_propellant_masses.append(acs_propellant_mass)

    # write the results to a file
    with open("sensitivity_analysis_results.txt", "w") as f:
        f.write("Landing DeltaV:\n")
        f.write(str(landing_delta_v_values) + "\n\n")
        f.write("Gross Masses:\n")
        f.write(str(gross_masses) + "\n\n")
        f.write("Dry Masses:\n")
        f.write(str(dry_masses) + "\n\n")
        f.write("Total Hydrogen Masses:\n")
        f.write(str(total_hydrogen_masses) + "\n\n")
        f.write("Total Oxygen Masses:\n")
        f.write(str(total_oxygen_masses) + "\n\n")
        f.write("ACS Propellant Masses:\n")
        f.write(str(acs_propellant_masses) + "\n")
