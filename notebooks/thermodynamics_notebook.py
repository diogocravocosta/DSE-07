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
    - amount of heat soaked by the vehicle during reentry
    """
    )
    return


@app.cell
def _(
    coolant,
    coolant_initial_temperature,
    coolant_mass,
    coolant_pressure,
    entry_speed,
    heat_soak_fraction,
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
    | Vehicle Dry Mass | {vehicle_dry_mass} | {vehicle_dry_mass.value} kg |
    | Coolant Mass | {coolant_mass} | {coolant_mass.value} kg |
    | Coolant Initial Temperature | {coolant_initial_temperature} | {coolant_initial_temperature.value} K |
    | Coolant Pressure | {coolant_pressure} | {round(coolant_pressure.value*1e-5, 1)} bar |
    | Entry Speed | {entry_speed} | {entry_speed.value} m/s |
    | Heat Soak Fraction | {heat_soak_fraction} | {heat_soak_fraction.value} - |
    | Maximum Heat Shield Temperature | {target_temperature} | {target_temperature.value} K |
    | Coolant Type | {coolant} | - |
    """
    )
    return


@app.cell
def _(
    Fluid,
    Input,
    coolant,
    coolant_initial_temperature,
    coolant_pressure,
    target_temperature,
):
    # Calculate the specific heating energy of the coolant
    _coolant = Fluid(coolant.value).with_state(Input.pressure(coolant_pressure.value), Input.temperature(coolant_initial_temperature.value))

    initial_internal_energy = _coolant.internal_energy
    final_internal_energy = _coolant.heating_to_temperature(temperature=target_temperature.value).internal_energy
    specific_heating_energy = final_internal_energy - initial_internal_energy
    specific_heating_energy
    return (specific_heating_energy,)


@app.cell
def _(
    change_plot_style,
    convert_fig_to_svg,
    coolant_mass,
    entry_speed,
    heat_soak_fraction,
    np,
    plot_svg,
    plt,
    specific_heating_energy,
    vehicle_dry_mass,
):
    # Mass-stepping simulation
    dm = 1  # mass step
    masses = np.linspace(vehicle_dry_mass.value + coolant_mass.value, vehicle_dry_mass.value, int(coolant_mass.value / dm) + 1)
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

    _fig, _ax1 = plt.subplots(1, 1, figsize=(10, 4))

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
def _(FluidsList, mo):
    # sliders
    vehicle_dry_mass = mo.ui.slider(10e3, 150e3, 1e3, value=60e3)
    coolant_mass = mo.ui.slider(1e3, 10e3, 1e2, value=5e3)
    coolant_pressure = mo.ui.slider(1e5, 15e5, 0.1e5, value=1.0e5)
    entry_speed = mo.ui.slider(6000, 8000, 10, value=7000)
    heat_soak_fraction = mo.ui.slider(0.0, 1.0, 0.01, value=0.02)
    target_temperature = mo.ui.slider(300, 2000, 100, value=1500)
    coolant = mo.ui.dropdown({
        "Hydrogen": FluidsList.Hydrogen,
        "Methane": FluidsList.Methane,
        "Water": FluidsList.Water,
        "Liquid Oxygen": FluidsList.Oxygen,
    }, value="Hydrogen")
    return (
        coolant,
        coolant_mass,
        coolant_pressure,
        entry_speed,
        heat_soak_fraction,
        target_temperature,
        vehicle_dry_mass,
    )


@app.cell
def _(FluidsList, coolant, mo):
    if coolant.value == FluidsList.Hydrogen:
        coolant_initial_temperature = mo.ui.slider(13.8, 30, 0.1, value=13.8)
    elif coolant.value == FluidsList.Methane:
        coolant_initial_temperature = mo.ui.slider(90.7168, 111.7, 0.1, value=91.7)
    elif coolant.value == FluidsList.Water:
        coolant_initial_temperature = mo.ui.slider(273.153, 293.15, 1, value=273.153)
    elif coolant.value == FluidsList.Oxygen:
        coolant_initial_temperature = mo.ui.slider(54.3705, 111.7, 0.1, value=54.3705)
    else:
        raise ValueError(f"Coolant {coolant.value} not supported")
    return (coolant_initial_temperature,)


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
    ### Coolant Fluid Properties

    The below plots provide some information derived from CoolProp.
    """
    )
    return


@app.cell
def _(
    Fluid,
    Input,
    change_plot_style,
    convert_fig_to_svg,
    coolant,
    coolant_initial_temperature,
    coolant_pressure,
    np,
    plot_svg,
    plt,
    target_temperature,
):
    # assume a nice low
    _coolant_temperature = coolant_initial_temperature.value
    _coolant_pressure = coolant_pressure.value
    _coolant = Fluid(coolant.value).with_state(Input.pressure(coolant_pressure.value), Input.temperature(coolant_initial_temperature.value))

    # target heating to 1400 K
    _target_temperature = target_temperature.value

    # energy needed to heat subcooled hydrogen to target temperature
    _initial_internal_energy = _coolant.internal_energy
    _final_internal_energy = _coolant.heating_to_temperature(_target_temperature).internal_energy
    _heating_energy = _final_internal_energy - _initial_internal_energy

    print("The density of subcooled coolant is:", round(_coolant.density, 3), "kg/m^3")
    print("The internal energy of coolant hydrogen is:", round(_coolant.internal_energy, 3), "kJ/kg")
    print("The energy needed to bring subcooled coolant to target temperature is:", round(_heating_energy, 3), "kJ/kg")

    # plot hydrogen properties as a function of temperature
    _temperatures = np.linspace(_coolant_temperature, _target_temperature, int(_target_temperature - _coolant_temperature + 1))
    _internal_energies = np.empty(_temperatures.shape)  # internal energy [kJ/kg]
    _specific_heats = np.empty(_temperatures.shape)  # mass specific constant pressure specific heat [J/kg/K].
    _conductivities = np.empty(_temperatures.shape)  # thermal conductivity [W/m/K]

    for _i, _T in enumerate(_temperatures):
        # Create a new fluid object for each temperature
        _fluid = _coolant.with_state(Input.temperature(_T), Input.pressure(_coolant_pressure))
        # Calculate internal energy and specific heat
        _internal_energies[_i] = _fluid.internal_energy
        _specific_heats[_i] = _fluid.specific_heat
        _conductivities[_i] = _fluid.conductivity

    _fig, (_ax1, _ax2, _ax3) = plt.subplots(3, 1, figsize=(10, 9))

    # fig.suptitle("Hydrogen Properties vs Temperature", fontsize=16)

    _ax1.plot(_temperatures, _internal_energies*1e-3)  # Convert J/kg to kJ/kg
    # _ax1.set_xlabel("Temperature (K)")
    _ax1.set_ylabel("Internal Energy (kJ/kg)")
    # _ax1.set_title("Internal Energy of Hydrogen vs Temperature")
    _ax1.grid()

    _ax2.plot(_temperatures, _specific_heats*1e-3)  # Convert J/kg/K to kJ/kg/K
    # _ax2.set_xlabel("Temperature (K)")
    _ax2.set_ylabel("Specific Heat (kJ/kg*K)")
    # _ax2.set_title("Specific Heat of Hydrogen vs Temperature")
    _ax2.grid()

    _ax3.plot(_temperatures, _conductivities)  # W/m/K
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
    mo.md(
        rf"""
    Heating energies assuming

    - Target temperature: 1400 K
    - Initial temperature: subcooled to freezing point

    are:

    - Hydrogen: 14,627,289.8 J/kg
    - Methane: 4,648,808.6 J/kg
    - Oxygen: 1,256,278.5 J/kg
    - Water: 4,315,890.9 J/kg

    **Hydrogen needs {round(14627289.8 / 4648808.6, 2)} times more energy to reach the heatshield temperature.**
    """
    )
    return


@app.cell
def _(mo):
    mo.md(
        r"""
    ### Vehicle dry mass as heat sink

    For this analysis, the entire vehicle shall be assumed to be made up of one metallic material with constant specific heat capacity and heated up to the heat shield temperature. This amount of energy shall be compared to the internal energy increase of the coolant.
    """
    )
    return


@app.cell
def _(coolant_initial_temperature, target_temperature, vehicle_dry_mass):
    _material = "Aluminium"
    _material_specific_heat = 600  # J/kg*K
    _material_mass = vehicle_dry_mass.value * 0.8  # correction to remove engine mass etc.
    _temperature_difference = target_temperature.value - coolant_initial_temperature.value

    _energy_needed = round(_material_mass * _material_specific_heat * _temperature_difference, 3)
    _energy_needed
    return


@app.cell
def _(mo):
    mo.md(
        rf"""
    Aluminum Al 2195 (Space Shuttle), assuming constant $c_p = 0.9$ kJ/kg*K, vehicle dry mass 60,000 kg, heat shield temperature 300 K: 12,363,840,000.0 J, this is roughly 12.4 GJ.

    Stainless steel (304L), assuming constant $c_p = 0.6$ kJ/kg*K, vehicle dry mass 60,000 kg, heat shield temperature 1400 K: 39,922,560,000.0 J, this is roughly 39.9 GJ

    For both materials, only 80% of the vehicle dry mass is assumed to be the heated mass.

    Note that hydrogen has 14,627,289.8 J/kg of heating energy, so roughly 14.6 MJ/kg. Therefore for stainless steel, the steel is equivalent to roughly 3000 kg of liquid hydrogen!

    Moreover, for stainless steel, one can consider the heat being radiated out.
    """
    )
    return


@app.cell
def _(mo):
    mo.md(
        r"""
    The Stefan-Boltzmann law is:

    $$ \frac{P}{A} = \sigma T^4 $$

    Based on the assumptions in the code, the energy radiated out is roughly 67.5 GJ.
    """
    )
    return


@app.cell
def _(entry_speed, np, vehicle_dry_mass):
    _sigma = 5.6703e-8  # Stefan-Boltzmann constant
    _radiation_temperature = 1400  # K

    # to calculate surface area, assume a cylinder
    _height = 20  # m
    _diameter = 7  # m
    _area = np.pi * _diameter * _height + 2 * np.pi * _diameter**2 / 4  # m^2

    # reentry time
    _time = 600  # s

    _energy = _sigma * _radiation_temperature**4 * _area * _time

    round(_energy*1e-9, 1), round(1/2*vehicle_dry_mass.value*entry_speed.value**2*1e-9, 1)
    return


@app.cell
def _(mo):
    mo.md(r"""So in conclusion, if the vehicle made out of stainless steel and it is allowed to operate at 1400 K, approximately 100 GJ of heat shall be handled by the structure itself.""")
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


@app.cell
def _(Fluid, FluidsList, Input):
    # from pyfluids import Fluid, FluidsList, Input

    # define the fluid
    _hydrogen = Fluid(FluidsList.Hydrogen).with_state(Input.pressure(1e5), Input.temperature(13.8))  # 1 bar, 13.8 K

    # assume the hydrogen is heated to 20 °C
    _initial_energy = _hydrogen.internal_energy  # internal energy at 13.8 K
    _final_energy = _hydrogen.heating_to_temperature(temperature=293.15)  # internal energy at 20 °C
    _heating_energy = _final_energy.internal_energy - _initial_energy  # J/kg

    # calculate the boiled-off mass
    _total_heat_flux = 3000  # W/m^2
    _surface_area = 70.6  # m^2
    _time = 24*3600  # s
    _total_heat = _total_heat_flux * _surface_area * _time  # J
    _boiled_off_mass = _total_heat / _heating_energy  # kg

    print("Heating energy of hydrogen from 13.8 K to 20 °C:", _heating_energy, "J/kg")
    print("Boiled-off mass of hydrogen:", _boiled_off_mass, "kg")
    return


if __name__ == "__main__":
    app.run()
