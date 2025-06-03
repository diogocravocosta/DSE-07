import marimo

__generated_with = "0.13.15"
app = marimo.App(width="medium")


@app.cell
def _(mo):
    mo.md(
        r"""
    | Variable | Adjustment | Value |
    |:---|:---:|:---|
    | Coolant Initial Temperature | {coolant_inlet_temperature} | {coolant_inlet_temperature.value} K |
    | Coolant Pressure | {coolant_inlet_pressure} | {round(coolant_inlet_pressure.value*1e-5, 1)} bar |


    | Parameter                 | Hess & Kunz    | McCarthy <br> & Wolf | Miller <br> Seader <br> & Trebes | Taylor           |
    |:--------------------------|:---------------|:---------------------|:---------------------------------|:-----------------|
    | \(X/D\) (dimensionless)   | *              | 6-50                 | 5-47                             | 2-252            |
    | \(T_s/T_b\) (dimensionless) | 3-14           | 2-11                 | 2-28                             | 1-23             |
    | \(T_{inlet}\) (K)         | 31,67 - 66,67  | 75,00 - 311,11       | 28,33 - 38,33                    | 25,00 - 100,00   |
    | Pressure (MPa)            | 1,620 - 5,137  | 0,221 - 9,335        | 3,158 - 17,140                   | 3,661 - 17,237   |
    | Heat Flux <br> (MW/m²)   | 1,000 - 16,613 | 0,057 - 24,194       | 2,237 - 39,250                   | 0,057 - 45,131   |
    | Mass Flow Rate <br> (kg/s) | *              | 0,000454 - <br> 0,058967 | 0,0907 - 0,3175                  | *                |
    | \(T_s\) (K)               | *              | 461,11 - 1244,44     | 59,44 - 961,11                   | 7,78 - 3127,78   |
    | ID (mm)                   | 5,08 - 7,62    | 5,08 - 10,16         | 5,08                             | *                |

    \* *not given*
    """
    )
    return


@app.cell
def _(mo):
    mo.md(
        r"""
    # Coolant heat flux

    The coolant heat flux is based on the Nusselt number.

    $$ h_g = 0.023 \frac{(\rho v)^{0.8}}{D^{0.2}} Pr^{0.4} \frac{\kappa}{\mu^{0.8}} $$

    The most basic, Dittus-Boelter equation is:

    $$ \text{Nu} = 0.023\, \text{Re}^{0.8} \text{Pr}^{0.4}$$

    Where $\text{Nu}$ is the Nusselt number, a dimensionless number that relates the convective heat transfer to the conductive heat transfer across a boundary, $\text{Re}$ is the Reynolds number, a dimensionless number that relates the inertial forces to the viscous forces, and $\text{Pr}$ is the Prandtl number, a dimensionless number that relates the momentum diffusivity (viscosity) to the thermal diffusivity.

    For liquid hydrogen, there are several empirical relations available. Based on the summary table, the Taylor relationship is deemed the most accurate. It is provided below:

    $$ \text{Nu}_b = 0.023\, \text{Re}_b^{0.8} \text{Pr}_b^{0.4} (T_s / T_b)^{-\left(0.57 - \frac{1.59}{x/D}\right)} $$

    The $b$ subscript indicates bulk properties.
    """
    )
    return


app._unparsable_cell(
    r"""
    def get_nusselt_number_taylor(reynolds_number, prandtl_number, temperature_ratio, dimensionless_length) -> float:
        \"\"\"
        Calculate the Nusselt number based on the Taylor relationship.
    
        Parameters:
        - reynolds_number: Reynolds number (dimensionless)
        - prandtl_number: Prandtl number (dimensionless)
        - temperature_ratio: Ratio of surface temperature to bulk temperature (dimensionless)
        - dimensionless_length: Dimensionless length (x/D, where x is the length and D is the diameter)
    
        Returns:
        - Nusselt number calculated using the Taylor empirical relationship
        \"\"\"
        return 0.023 * (reynolds_number ** 0.8) * (prandtl_number ** 0.4) * (temperature_ratio ** (-0.57 - (1.59 / dimensionless_length)))


    def get_reynolds_number(fluid_density, fluid_speed, hydraulic_diameter, dynamic_viscosity)

    _coolant = Fluid(FluidsList.Hydrogen).with_state(Input.temperature(coolant_inlet_temperature.value), Input.pressure(coolant_inlet_pressure.value))
    _temperature_ratio = 0.55  # sample value
    _dimensionless_length = 1.0  # sample value
    """,
    name="_"
)


@app.cell
def _(change_plot_style, convert_fig_to_svg, plot_svg, plt):
    _fig, _ax1 = plt.subplots(1, 1, figsize=(10, 4))

    _ax1.plot([1,2,3], [4,5,6])  # Convert J/kg to kJ/kg
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
    coolant_inlet_temperature = mo.ui.slider(13.8, 50, 0.1, value=13.8)
    coolant_inlet_pressure = mo.ui.slider(1e5, 50e5, 0.1e5, value=10e5)
    return


@app.cell
def _(mo):
    mo.md(r"""# Helpers""")
    return


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
def _():
    from pyfluids import Fluid, FluidsList, Input
    import io
    import matplotlib.pyplot as plt
    import marimo as mo
    import numpy as np
    return io, mo, plt


if __name__ == "__main__":
    app.run()
