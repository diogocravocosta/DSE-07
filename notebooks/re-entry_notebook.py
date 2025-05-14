import marimo

__generated_with = "0.13.8"
app = marimo.App(width="medium")


@app.cell
def _(mo):
    mo.md(
        r"""
    # Re-entry Simulations Notebook

    This notebook explores the re-entry of a spacecraft into the atmosphere. It performs trajectory calculations and heating calculations.
    """
    )
    return


@app.cell
def _(mo):
    mo.md(
        r"""
    ## Ballistic Entry

    Ballistic entry calculations are implemented per Section 5.3 of the Re-entry Systems book by Erwin Mooij.

    The analytical solution for the normalised velocity of the ballistic re-entry problem is given by:

    $$ \frac{V}{V_E} = \exp \left(\frac{1}{2} \frac{g \rho}{K \beta \sin \gamma_E}\right)$$

    $K$ is the ballistic coefficient and $\beta$ represents the exponential atmosphere density model. It can also be written as $\beta = \frac{1}{H_s}$ where $H_s$ is the scale height. Therefore the above equation can be rewritten as:

    $$ \frac{V}{V_E} = \exp \left(\frac{1}{2} \frac{g \rho H_s}{K \sin \gamma_E}\right)$$

    Generally, the following three equations of motion are applicable:

    $$ m \frac{dV}{dt} = -D - mg \sin \gamma $$

    $$ mV \frac{d\gamma}{dt} = L - mg \cos \gamma \left(1 - \frac{V^2}{V_c^2}\right) $$

    $$ \frac{dR}{dt} = \frac{dh}{dt} = V \sin \gamma $$

    The third equation is a kinematic relation.

    ### Deceleration

    Maximum deceleration is given by:

    $$ \bar{a}_{\text{max}} = - \left( \frac{dV}{dt} \right)_{\text{max}} = - \frac{\beta \sin \gamma_E}{2e} V_E^2 $$

    ### Ballistic parameter

    It can be defined as 
    """
    )
    return


@app.cell
def _(
    ballistic_coefficient,
    entry_flight_path_angle,
    entry_speed,
    mo,
    scale_height,
):
    mo.md(
        f"""
    ### Ballistic Entry Parameters

    | Parameter | Adjustment | Value |
    |:---|:---:|:---|
    | Ballistic Coefficient | {ballistic_coefficient} | {ballistic_coefficient.value} kg/ ... |
    | Entry Speed | {entry_speed} | {entry_speed.value} m/s |
    | Entry Flight Path Angle | {entry_flight_path_angle} | {entry_flight_path_angle.value} degrees |
    | Scale Height | {scale_height} | {scale_height.value} m |
    """
    )
    return


@app.cell
def _(change_plot_style, convert_fig_to_svg, plot_svg, plt):
    # the tank is plotted as two vertical lines and two semicircles.
    fig, ax = plt.subplots()
    # ax.plot()

    ax.set_xlim(-1, 11)
    ax.set_ylim(-6, 60)
    ax.set_xlabel("x [m]")
    ax.set_ylabel("y [m]")
    ax.set_title("V/V_e ...")
    ax.grid(True)
    plt.close(fig)  # Close the plot to prevent default PNG rendering

    # Change plotting style
    change_plot_style()

    # Display the SVG in Marimo using mo.Html
    svg_data = convert_fig_to_svg(fig)
    plot_svg(svg_data)
    return


@app.cell
def _(mo):
    # sliders
    ballistic_coefficient = mo.ui.slider(0, 10000, 100, value=2000)
    entry_speed = mo.ui.slider(6e3, 10e3, 100, value=8e3)
    entry_flight_path_angle = mo.ui.slider(-90, 0, -0.1, value=-10)
    scale_height = mo.ui.slider(0, 10000, 100, value=7200)
    return (
        ballistic_coefficient,
        entry_flight_path_angle,
        entry_speed,
        scale_height,
    )


@app.cell
def _(np):
    def get_ballistic_normalised_velocity(
        ballistic_coefficient: float,
        entry_speed: float,
        entry_flight_path_angle: float,
        scale_height: float,
    ) -> float:
        """
        Calculate the normalised velocity for a ballistic entry.

        Parameters:
        ----------
        ballistic_coefficient : float
            Ballistic coefficient of the spacecraft.
        entry_speed : float
            Entry speed of the spacecraft.
        entry_flight_path_angle : float
            Entry flight path angle of the spacecraft.
        scale_height : float
            Scale height of the atmosphere.

        Returns:
        -------
        float
            Normalised velocity.
        """
        # Constants
        g = 9.81  # m/s^2
        rho = 1.225  # kg/m^3

        # Calculate normalised velocity
        normalised_velocity = np.exp(
            (0.5 * g * rho * scale_height) / (ballistic_coefficient * np.sin(np.radians(entry_flight_path_angle)))
        )

        return normalised_velocity
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
    import marimo as mo
    import numpy as np
    import matplotlib.pyplot as plt
    import io
    return io, mo, np, plt


if __name__ == "__main__":
    app.run()
