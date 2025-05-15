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
    - the mass-stepping simulation of the energy exchange between the vehicle and the liquid hydrogen
    - amount of heat soaked by the vehicle during re-entry
    """
    )
    return


@app.cell
def _(
    entry_speed,
    heat_soak_fraction,
    hydrogen_mass,
    mo,
    target_temperature,
    vehicle_dry_mass,
):
    mo.md(
        rf"""
    ## Simulation Description

    The simulation is stepping in mass. If the vehicle starts with 60 tonnes of dry mass and 40 tonnes of propellant mass, it first calculates the kinetic energy at entry speed of the two masses combined. Then, at each step, it is assumed that a small amount of hydrogen is heated to a target temperature that corresponds to the steady state heat shield temperature (the heat shield is assumed regeneratively cooled) and that this heating energy is taken from the kinetic energy. In practice, only a small portion of the lost kinetic energy is soaked as heat into the vehicle, so a parameter is defined that account for this.

    Next, the heated hydrogen is vented from the vehicle. Effectively, the vehicle loses even more kinetic energy because it reduced in total mass. Then, another step is started.

    In summary, these are the steps of the simulation

    1. Calculate the initial total kinetic energy of the vehicle and the propellant using its entry speed
    2. Calculate the heating energy of the hydrogen

    At the end of the simulation, the two energies are combined answering the question: how much energy is used to heat up the hydrogen and how much kinetic energy is lost.

    | Variable | Adjustment | Value |
    |:---|:---:|:---|
    | Vehicle Dry Mass | {vehicle_dry_mass} | {vehicle_dry_mass.value} |
    | Hydrogen Mass | {hydrogen_mass} | {hydrogen_mass.value} |
    | Entry Speed | {entry_speed} | {entry_speed.value} |
    | Heat Soak Fraction | {heat_soak_fraction} | {heat_soak_fraction.value} |
    | Maximum Heat Shield Temperature | {target_temperature} | {target_temperature.value} |
    """
    )
    return


@app.cell
def _(
    Fluid,
    FluidsList,
    Input,
    hydrogen_initial_temperature,
    hydrogen_pressure,
    target_temperature,
):
    # Calculate the specific heating energy of the hydrogen
    hydrogen = Fluid(FluidsList.Hydrogen).with_state(Input.pressure(hydrogen_pressure.value), Input.temperature(hydrogen_initial_temperature.value))

    initial_internal_energy = hydrogen.internal_energy
    final_internal_energy = hydrogen.heating_to_temperature(temperature=target_temperature.value).internal_energy
    specific_heating_energy = final_internal_energy - initial_internal_energy
    return (specific_heating_energy,)


@app.cell
def _(
    change_plot_style,
    convert_fig_to_svg,
    entry_speed,
    heat_soak_fraction,
    hydrogen_mass,
    np,
    plot_svg,
    plt,
    specific_heating_energy,
    vehicle_dry_mass,
):
    # Mass-stepping simulation
    dm = 1  # mass step
    masses = np.linspace(vehicle_dry_mass.value + hydrogen_mass.value, vehicle_dry_mass.value, int(hydrogen_mass.value / dm) + 1)
    kinetic_energies = np.empty(masses.shape)
    velocities = np.empty(masses.shape)

    velocities[0] = entry_speed.value

    for _i, _m in enumerate(masses[:-1]):
        # Calculate the kinetic energy of the vehicle at the current mass step
        kinetic_energies[_i] = 0.5 * _m * velocities[_i]**2

        # Calculate the energy lost by the vehicle due to the heating of the hydrogen
        _heating_energy = specific_heating_energy * dm

        # Correct for the heat soak fraction, more kinetic energy is actually lost than the heat soaked by the hydrogen
        _heating_energy = _heating_energy / heat_soak_fraction.value

        # Calculate the change in velocity as a result of the heating energy
        velocities[_i+1]  = np.sqrt(2 * (kinetic_energies[_i] - _heating_energy) / _m)

    _fig, _ax1 = plt.subplots(1, 1, figsize=(10, 5))

    _ax1.plot(masses, velocities)  # Convert J/kg to kJ/kg
    _ax1.set_xlabel("Masses (kg)")
    _ax1.set_ylabel("Velocities (m/s)")
    _ax1.set_title("Mass vs Velocity")
    _ax1.grid()

    plt.close(_fig)  # Close the plot to prevent default PNG rendering

    # Change plotting style
    change_plot_style()

    # Display the SVG in Marimo using mo.Html
    _svg_data = convert_fig_to_svg(_fig)
    plot_svg(_svg_data)
    return


@app.cell
def _(mo):
    # sliders
    vehicle_dry_mass = mo.ui.slider(10e3, 100e3, 1e3, value=60e3)
    hydrogen_mass = mo.ui.slider(1e3, 50e3, 1e3, value=40e3)
    hydrogen_initial_temperature = mo.ui.slider(13.8, 20, 0.1, value=13.8)
    hydrogen_pressure = mo.ui.slider(1e5, 10e5, 0.1e5, value=1.0e5)
    entry_speed = mo.ui.slider(6000, 8000, 10, value=7000)
    heat_soak_fraction = mo.ui.slider(0.0, 1.0, 0.01, value=0.05)
    target_temperature = mo.ui.slider(500, 2000, 100, value=1500)

    return (
        entry_speed,
        heat_soak_fraction,
        hydrogen_initial_temperature,
        hydrogen_mass,
        hydrogen_pressure,
        target_temperature,
        vehicle_dry_mass,
    )


@app.cell
def _(mo):
    mo.md(r"""## Theory""")
    return


@app.cell
def _(mo):
    mo.md(
        r"""
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
    - https://physics.stackexchange.com/questions/689218/how-can-hydrogen-have-negative-enthalpy

    ### Videos

    - https://www.youtube.com/watch?v=7BA7iVTRyO4
    - https://www.youtube.com/watch?v=LZX8mlNRx2c

    ### Tasks

    - [ ] Implement a vaporization-stepping simulation where the energy exchange between liquid hydrogen and vehicle kinetic energy happens in small steps
    - [ ] Implement CoolProp to evaluate energy needed to heat, vaporize and heat again to heat sink temperature
    - [ ] Implement the vehicle mass as a heat sink
    """
    )
    return


@app.cell
def _(mo):
    mo.md(
        r"""
    ### Liquid Hydrogen Properties

    The below plots provide some information derived from CoolProp.
    """
    )
    return


@app.cell
def _(
    Fluid,
    FluidsList,
    Input,
    change_plot_style,
    convert_fig_to_svg,
    np,
    plot_svg,
    plt,
):
    # assume subcooled hydrogen stored at 1 bar
    _hydrogen_temperature = 13.8
    _hydrogen_pressure = 1e5
    subcooled_hydrogen = Fluid(FluidsList.Hydrogen).with_state(Input.pressure(_hydrogen_pressure), Input.temperature(_hydrogen_temperature))

    # target heating to 1000 K
    _target_temperature = 1000

    # energy needed to heat subcooled hydrogen to target temperature
    _initial_internal_energy = subcooled_hydrogen.internal_energy
    _final_internal_energy = subcooled_hydrogen.heating_to_temperature(_target_temperature).internal_energy
    _heating_energy = _final_internal_energy - _initial_internal_energy

    print("The density of subcooled hydrogen is:", round(subcooled_hydrogen.density, 3), "kg/m^3")
    print("The internal energy of subcooled hydrogen is:", round(subcooled_hydrogen.internal_energy, 3), "kJ/kg")
    print("The energy needed to bring subcooled hydrogen to target temperature is:", round(_heating_energy, 3), "kJ/kg")

    # plot hydrogen properties as a function of temperature
    temperatures = np.linspace(13.8, 1000, 1000)
    internal_energies = np.empty(temperatures.shape)  # internal energy [kJ/kg]
    specific_heats = np.empty(temperatures.shape)  # mass specific constant pressure specific heat [J/kg/K].
    conductivities = np.empty(temperatures.shape)  # thermal conductivity [W/m/K]

    for _i, _T in enumerate(temperatures):
        # Create a new fluid object for each temperature
        fluid = Fluid(FluidsList.Hydrogen).with_state(Input.temperature(_T), Input.pressure(_hydrogen_pressure))
        # Calculate internal energy and specific heat
        internal_energies[_i] = fluid.internal_energy
        specific_heats[_i] = fluid.specific_heat
        conductivities[_i] = fluid.conductivity

    _fig, (_ax1, _ax2, _ax3) = plt.subplots(3, 1, figsize=(10, 9))

    # fig.suptitle("Hydrogen Properties vs Temperature", fontsize=16)

    _ax1.plot(temperatures, internal_energies*1e-3)  # Convert J/kg to kJ/kg
    # _ax1.set_xlabel("Temperature (K)")
    _ax1.set_ylabel("Internal Energy (kJ/kg)")
    # _ax1.set_title("Internal Energy of Hydrogen vs Temperature")
    _ax1.grid()

    _ax2.plot(temperatures, specific_heats*1e-3)  # Convert J/kg/K to kJ/kg/K
    # _ax2.set_xlabel("Temperature (K)")
    _ax2.set_ylabel("Specific Heat (kJ/kg*K)")
    # _ax2.set_title("Specific Heat of Hydrogen vs Temperature")
    _ax2.grid()

    _ax3.plot(temperatures, conductivities)  # W/m/K
    # _ax3.set_xlabel("Temperature (K)")
    _ax3.set_ylabel("Thermal Conductivity (W/m*K)")
    # _ax3.set_title("Thermal Conductivity of Hydrogen vs Temperature")
    _ax3.grid()

    plt.close(_fig)  # Close the plot to prevent default PNG rendering

    # Change plotting style
    change_plot_style()

    # Display the SVG in Marimo using mo.Html
    _svg_data = convert_fig_to_svg(_fig)
    plot_svg(_svg_data)
    return


@app.cell
def _(mo):
    mo.md(r"""# Helpers""")
    return


@app.cell
def _():
    import io
    import matplotlib.pyplot as plt
    import marimo as mo
    import numpy as np
    from pyfluids import Fluid, FluidsList, Input
    return Fluid, FluidsList, Input, io, mo, np, plt


@app.cell
def _(io, mo):
    def change_plot_style():
        """ this breaks the web version
        plt.rcParams.update({
            "text.usetex": True,
            "font.family": "Computer Modern Roman"
        })
        """
        pass

    def convert_fig_to_svg(fig):
        # Save the plot to an in-memory SVG buffer
        svg_buffer = io.StringIO()
        fig.savefig(svg_buffer, format="svg")
        svg_buffer.seek(0)
        svg_data = svg_buffer.getvalue()
        return svg_data

    def plot_svg(svg_data):
        # Display the SVG in Marimo using mo.Html
        return mo.Html(f"""
            <div>
                {svg_data}
            </div>
        """)

    return change_plot_style, convert_fig_to_svg, plot_svg


if __name__ == "__main__":
    app.run()
