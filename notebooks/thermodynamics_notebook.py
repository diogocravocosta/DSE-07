import marimo

__generated_with = "0.13.6"
app = marimo.App(width="medium")


@app.cell
def _(mo):
    mo.md(
        r"""
    # Thermodynamics Notebook

    This notebook will estimate

    - the kinetic energy of a re-entering second stage
    - the thermal energy the liquid hydrogen can contain

    Then the two numbers will be compared to estimate the effectiveness of an active cooling re-entry system.

    Kinetic energy: $E_k = \frac{1}{2} m V^2$

    ### Fluid properties of liquid hydrogen

    Assuming liquid hydrogen at 20 K and at 1 bar:

    - density $\rho = 70.85\space kg/m^3$
    - specific heat at constant pressure $c_p = 9.56\space kJ/kg \cdot K$
    - latent heat of vaporization $L_v = 223.2\space kJ/kg$

    Reference data is from NIST:
    - https://webbook.nist.gov/cgi/fluid.cgi?T=20&PLow=0.5&PHigh=2&PInc=0.1&Digits=5&ID=C1333740&Action=Load&Type=IsoTherm&TUnit=K&PUnit=bar&DUnit=kg%2Fm3&HUnit=kJ%2Fkg&WUnit=m%2Fs&VisUnit=uPa*s&STUnit=N%2Fm&RefState=DEF
    - https://www.nuclear-power.com/hydrogen-specific-heat-latent-heat-vaporization-fusion/

    Ortho-para configuration?

    ### Fluid properties of liquid methane

    https://webbook.nist.gov/cgi/fluid.cgi?T=91&PLow=1&PHigh=2&PInc=1&Digits=5&ID=C74828&Action=Load&Type=IsoTherm&TUnit=K&PUnit=bar&DUnit=kg%2Fm3&HUnit=kJ%2Fkg&WUnit=m%2Fs&VisUnit=uPa*s&STUnit=N%2Fm&RefState=DEF

    ### Further reading

    - https://www.reddit.com/r/spacex/comments/a9y9r0/an_energy_budget_for_starship_reentry/
    - https://tfaws.nasa.gov/TFAWS12/Proceedings/Aerothermodynamics%20Course.pdf

    ### Videos

    - https://www.youtube.com/watch?v=7BA7iVTRyO4
    - https://www.youtube.com/watch?v=LZX8mlNRx2c
    """
    )
    return


@app.cell
def _():
    def calculate_kinetic_energy(mass: float, velocity: float) -> float:
        return 0.5 * mass * velocity**2

    def calculate_heating_energy(mass: float, specific_heat: float, temperature_change: float) -> float:
        return mass * specific_heat * temperature_change

    def calculate_vaporization_energy(mass: float, latent_heat: float) -> float:
        return mass * latent_heat
    return (
        calculate_heating_energy,
        calculate_kinetic_energy,
        calculate_vaporization_energy,
    )


@app.cell
def _(
    calculate_heating_energy,
    calculate_kinetic_energy,
    calculate_vaporization_energy,
    np,
):
    # --- Initial parameters from your example ---
    initial_system_mass = 60000.0  # kg (this is the mass contributing to KE)
    velocity = 6000.0  # m/s

    # This is the total amount of mass that will be targeted for vaporization
    # from the initial_system_mass.
    target_total_mass_to_vaporize = 40000.0 # kg

    # Material/process properties from your example
    # (used for calculate_heating_energy and calculate_vaporization_energy)
    specific_heat_capacity = 9.56e3  # J/(kg*K) or J/(kg*°C)
    delta_T_for_heating = 5.0  # K or °C (temperature change for the heating phase)
    latent_heat_of_vaporization = 223.2e3  # J/kg

    # --- Simulation control parameters ---
    # Define how much mass is processed (heated and vaporized) in each simulation step.
    # This determines the granularity and number of steps in the simulation.
    # Let's choose a value, e.g., 1000 kg per step.
    mass_chunk_to_process_per_step = 1000.0  # kg

    # --- Simulation variables ---
    current_system_mass = initial_system_mass
    total_mass_vaporized_in_simulation = 0.0
    step_count = 0

    # --- Initial State Output ---
    print(f"--- Initial Conditions ---")
    print(f"Initial System Mass: {current_system_mass:,.2f} kg".replace('.', ','))
    initial_kinetic_energy = calculate_kinetic_energy(current_system_mass, velocity)
    print(f"Initial Kinetic Energy: {initial_kinetic_energy:,.2e} J".replace('.', ','))
    print(f"Target Total Mass to Vaporize: {target_total_mass_to_vaporize:,.2f} kg".replace('.', ','))
    print(f"Mass Chunk per Step: {mass_chunk_to_process_per_step:,.2f} kg".replace('.', ','))
    print(f"--- Starting Time-Stepping Simulation ---")
    print("| Step | Mass Vaporized   | Remaining Mass   | Kinetic Energy   | Total Vaporized  | Heating Energy  | Vaporization Energy |")
    print("|------|------------------|------------------|------------------|------------------|-----------------|---------------------|")
    print("|      | this Step (kg)   | (kg)             | (J)              | (kg)             | this Step (J)   | this Step (J)       |")
    print("|------|------------------|------------------|------------------|------------------|-----------------|---------------------|")

    # --- Simulation Loop ---
    while total_mass_vaporized_in_simulation < target_total_mass_to_vaporize and current_system_mass > 0:
        step_count += 1

        # Determine the actual amount of mass to vaporize in this step
        mass_to_vaporize_this_step = mass_chunk_to_process_per_step

        # Adjust if it exceeds the remaining target amount
        if total_mass_vaporized_in_simulation + mass_to_vaporize_this_step > target_total_mass_to_vaporize:
            mass_to_vaporize_this_step = target_total_mass_to_vaporize - total_mass_vaporized_in_simulation

        # Adjust if it exceeds the current available system mass
        if mass_to_vaporize_this_step > current_system_mass:
            mass_to_vaporize_this_step = current_system_mass

        # If no mass can be vaporized (e.g., target met, or system empty), break the loop
        if mass_to_vaporize_this_step <= 1e-9: # Using a small epsilon to handle potential float issues
            break

        # Calculate energies required for processing this mass chunk
        heating_energy_this_step = calculate_heating_energy(mass_to_vaporize_this_step, specific_heat_capacity, delta_T_for_heating)
        vaporization_energy_this_step = calculate_vaporization_energy(mass_to_vaporize_this_step, latent_heat_of_vaporization)

        # Update system state: reduce mass and track total vaporized
        current_system_mass -= mass_to_vaporize_this_step
        total_mass_vaporized_in_simulation += mass_to_vaporize_this_step

        # Calculate the new kinetic energy with the reduced mass
        current_kinetic_energy = calculate_kinetic_energy(current_system_mass, velocity)

        # Derive the new velocity
        velocity = np.sqrt(2 * current_kinetic_energy / current_system_mass) if current_system_mass > 0 else 0

        # Print step details (using standard Python float formatting with '.', then replacing for display)
        # For J values, scientific notation might be more readable for large numbers.
        print(f"| {step_count:<4} | {mass_to_vaporize_this_step:>16,.2f} | {current_system_mass:>16,.2f} | {current_kinetic_energy:>16,.2e} | {total_mass_vaporized_in_simulation:>16,.2f} | {heating_energy_this_step:>15,.2e} | {vaporization_energy_this_step:>19,.2e} |".replace('.', ','))

        # Safety break: if system mass is depleted but target not fully met (should be covered by logic above)
        if current_system_mass <= 1e-9 and total_mass_vaporized_in_simulation < target_total_mass_to_vaporize - 1e-9:
            print("Warning: System mass depleted before achieving the full target vaporization amount.".replace('.', ','))
            break

    print("|------|------------------|------------------|------------------|------------------|-----------------|---------------------|")
    print(f"--- Simulation Ended ---")
    print(f"Total steps: {step_count}".replace('.', ','))
    print(f"Total mass vaporized: {total_mass_vaporized_in_simulation:,.2f} kg".replace('.', ','))
    print(f"Final system mass: {current_system_mass:,.2f} kg".replace('.', ','))

    final_kinetic_energy = calculate_kinetic_energy(current_system_mass, velocity) # Should be same as last current_kinetic_energy
    print(f"Final Kinetic Energy: {final_kinetic_energy:,.2e} J".replace('.', ','))

    if initial_kinetic_energy > 1e-9: # Avoid division by zero
        percentage_ke_reduction = (initial_kinetic_energy - final_kinetic_energy) / initial_kinetic_energy
        print(f"Total Kinetic Energy Reduction: {percentage_ke_reduction:.2%}".replace('.', ','))
    else:
        print("Initial kinetic energy was zero or negligible; percentage reduction cannot be calculated.".replace('.', ','))

    # --- Comparison with original single calculations (for context) ---
    # Original calculation for total energy to heat 40000 kg
    # total_heating_energy_original = calculate_heating_energy(40000.0, specific_heat_capacity, delta_T_for_heating)
    # Original calculation for total energy to vaporize 40000 kg
    # total_vaporization_energy_original = calculate_vaporization_energy(40000.0, latent_heat_of_vaporization)

    # print(f"\\n--- For Reference: Original Single Calculations (for 40,000 kg mass chunk) ---")
    # print(f"Total energy to heat 40,000 kg: {total_heating_energy_original:,.2e} J".replace('.', ','))
    # print(f"Total energy to vaporize 40,000 kg: {total_vaporization_energy_original:,.2e} J".replace('.', ','))
    # if initial_kinetic_energy > 1e-9:
    #     print(f"Original calculation's heating energy as % of initial KE: {total_heating_energy_original / initial_kinetic_energy:.2%}".replace('.', ','))
    #     print(f"Original calculation's vaporization energy as % of initial KE: {total_vaporization_energy_original / initial_kinetic_energy:.2%}".replace('.', ','))
    return


@app.cell
def _():
    import marimo as mo
    import numpy as np
    return mo, np


if __name__ == "__main__":
    app.run()
