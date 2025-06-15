import numpy as np
import itertools
import csv
import concurrent.futures

import data.material as dm
import data.constants as cn

import h2ermes_tools.variables as vr
import h2ermes_tools.integration.mass_integration as mi
from h2ermes_tools.integration.main import run_integration


# Helper to get a Variable object by name from vr
def get_var(name):
    return getattr(vr, name)


def add_unchanging_variables(integrator: mi.MassIntegrator):
    # ...existing code for unchanging variables...
    integrator.payload_mass = 8_000
    integrator.h2_power_mass = vr.hydrogen_power_mass.value
    integrator.o2_power_mass = vr.oxygen_power_mass.value
    integrator.acs_propellant_mass = vr.acs_propellant_mass.value
    integrator.coolant_mass = vr.coolant_mass.value
    integrator.sea_level_isp = vr.engine_specific_impulse_sl_opt.value
    integrator.vacuum_isp = 429.7477666666666
    integrator.sea_level_thrust_chambers_number = vr.n_chambers_sl.value
    integrator.vacuum_thrust_chambers_number = vr.n_chambers_vac.value
    integrator.landing_delta_v = vr.landing_delta_v.value
    integrator.deorbit_delta_v = vr.deorbit_delta_v.value + 50
    integrator.circularization_delta_v = vr.target_orbit_circularization_delta_v.value
    integrator.orbit_raising_delta_v = vr.orbit_raising_delta_v.value
    integrator.orbit_insertion_delta_v = vr.orbit_insertion_delta_v.value
    integrator.of_ratio = vr.engine_mixture_ratio.value
    integrator.vacuum_twr = vr.t_w_vac.value
    integrator.clearance_height = 2.5
    integrator.hydrogen_design_pressure = 2e5
    integrator.oxygen_design_pressure = 2.5e5
    integrator.boiloff_design_pressure = 10e5
    integrator.landing_leg_material = dm.Ti6Al4V
    integrator.tank_material = dm.random_steel
    integrator.header_tank_material = dm.random_steel
    integrator.g_reentry_force_ratio = 8
    integrator.g_launch_force_ratio = 6
    integrator.bottom_radius = 5
    integrator.unchanging_subsystem_dry_masses = {
        "acs": vr.acs_dry_mass.value,
        "heat shield": vr.heat_shield_mass.value,
        "thrust chambers": vr.thrust_chamber_sl_mass.value
        * integrator.sea_level_thrust_chambers_number
        + vr.thrust_chamber_vac_mass.value * integrator.vacuum_thrust_chambers_number,
        "docking system": vr.docking_system_mass.value,
        "nose cone": vr.nose_cone_mass.value,
        "power": vr.power_dry_mass.value,
        "avionics": vr.avionics_mass.value,
        "wiring": vr.wiring_mass.value,
        "interstage": vr.interstage_mass.value,
    }


def add_initial_values(integrator: mi.MassIntegrator):
    multiplier = 1.0
    integrator.dry_mass = vr.total_dry_mass.value * multiplier
    integrator.h2_boil_off_mass = 1000 * multiplier
    integrator.gross_mass = vr.gross_mass.value * multiplier
    integrator.total_vacuum_thrust = (
        integrator.gross_mass * cn.g_0 * integrator.vacuum_twr * multiplier
    )
    integrator.header_hydrogen_mass = 1000 * multiplier
    integrator.header_oxygen_mass = 5000 * multiplier
    integrator.phi_top = np.deg2rad(10)
    integrator.middle_radius = 4.8
    integrator.hydrogen_tank_height = 12.65
    integrator.top_radius = 2.5


def add_vars(integrator, variable_names, combo):
    add_unchanging_variables(integrator)
    for name, val in zip(variable_names, combo):
        setattr(integrator, name, val)


def single_run(args):
    combo, variable_names = args
    def add_vars_wrapper(integrator):
        add_vars(integrator, variable_names, combo)
    results = run_integration(add_initial_values, add_vars_wrapper)
    row = dict(zip(variable_names, combo))
    if isinstance(results, tuple):
        for i, val in enumerate(results):
            row[f"result_{i}"] = val
    else:
        row["result"] = results
    return row


def run_sensitivity(variable_names, n_steps=3, output_csv="sensitivity_results.csv"):
    # Build value grids for each variable
    value_grids = []
    for name in variable_names:
        var = get_var(name)
        values = np.linspace(min(var.margin), max(var.margin), n_steps) * var.value
        value_grids.append(values)
    combos = list(itertools.product(*value_grids))
    args_list = [(combo, variable_names) for combo in combos]

    # Parallel execution
    with concurrent.futures.ProcessPoolExecutor() as executor:
        results_list = list(executor.map(single_run, args_list))

    # Write to CSV
    if results_list:
        with open(output_csv, "w", newline="") as f:
            writer = csv.DictWriter(f, fieldnames=results_list[0].keys())
            writer.writeheader()
            writer.writerows(results_list)
        print(f"Results written to {output_csv}")
    else:
        print("No results to write.")


if __name__ == "__main__":
    # Example: sweep over hydrogen_power_mass and coolant_mass, 5 steps each
    run_sensitivity(["hydrogen_power_mass", "coolant_mass"], n_steps=3)
