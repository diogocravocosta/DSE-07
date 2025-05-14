import marimo

__generated_with = "0.13.6"
app = marimo.App(width="medium")


@app.cell
def _(mo):
    mo.md(
        r"""
    # Re-entry Gliding Entry Simulations Notebook

    This notebook explores the re-entry of a spacecraft into the atmosphere for the gliding mode. It performs trajectory calculations and heating calculations. It
    """
    )
    return


@app.cell
def _(mo):
    mo.md(
        r"""
    ## Gliding Entry

    Ballistic entry calculations are implemented per Section 5.4 of the Re-entry Systems book by Erwin Mooij.

    The analytical solution for the normalised velocity of the gliding re-entry problem is given by:

    $$ \frac{V}{V_{c,0}} = \sqrt{\frac{\frac{W/S}{C_L}}{\frac{1}{2} \rho_0 e^{-\beta h} V_{c,0}^2 + \frac{W/S}{C_L}}} $$


    $\frac{W/S}{C_L}$ is the lift parameter and $\beta$ represents the exponential atmosphere density model. It can also be written as $\beta = \frac{1}{H_s}$ where $H_s$ is the scale height. Therefore the above equation can be rewritten as:

    Generally, the following three equations of motion are applicable:

    $$ m \frac{dV}{dt} = -D - mg \sin \gamma $$

    $$ mV \frac{d\gamma}{dt} = L - mg \cos \gamma \left(1 - \frac{V^2}{V_c^2}\right) $$

    $$ \frac{dR}{dt} = \frac{dh}{dt} = V \sin \gamma $$

    ## Equilibrium Flight-Path Angle

    For a large part of the flight, altitude decreases at a constant rate with velocity, which corresponds with a constant flight-path angle.

    The equilibrium flight-path angle is given by:

    $$ \bar{\gamma} \approx \sin \bar{\gamma} = -\frac{1}{\beta R_e} \cdot \frac{2}{L/D} \cdot \frac{V_c^2}{V^2} $$

    ## Flight Range and Flight Duration

    The flight range is given by:

    $$ \frac{R_f}{R_e} = -\frac{1}{2} \frac{L}{D} \ln \left(1 - \frac{V_E^2}{V_c^2}\right) $$

    Flight time is given by:

    $$ t_{\text{flight}} = \frac{1}{2} \frac{V_c}{g} \frac{L}{D} \ln \left(\frac{1 + V_E/V_c}{1 - V_E/V_c}\right) $$
    """
    )
    return


@app.cell
def _(altitude, lift_drag_ratio, lift_parameter, mo, scale_height):
    mo.md(
        f"""
    ### Gliding Entry Parameters

    | Parameter | Adjustment | Value |
    |:---|:---:|:---|
    | Lift parameter | {lift_parameter} | {lift_parameter.value} [-]|
    | Lift-Drag Ratio | {lift_drag_ratio} | {lift_drag_ratio.value} [-] |
    | Scale Height | {scale_height} | {scale_height.value} m |
    | Altitude | {altitude} | {altitude.value} km |
    """
    )
    return


@app.cell
def _(
        V_Vc_ratio,
        height,
        flight_range_ratio,
        entry_circular_ratio,
        a_weight_amax_ratio,
        a_no_weight_amax_ratio,
        change_plot_style,
        convert_fig_to_svg,
        plot_svg,
        plt):
    # the tank is plotted as two vertical lines and two semicircles.
    fig, (ax1, ax2, ax3) = plt.subplots(3, 1, figsize=(5, 15))

    ax1.plot(V_Vc_ratio, height / 1000, label=r"$V/V_c$", color="blue")
    ax1.set_xlim(0.001, 1)
    ax1.set_ylim(0, 120)
    ax1.set_xlabel(r"$\frac{V}{V_c}$" + "[-]")
    ax1.set_ylabel(r"$h$ [km]")
    ax1.set_title("Velocity Ratio during gliding re-entry")
    ax1.grid(True)

    ax2.plot(entry_circular_ratio, flight_range_ratio, label=r"$R/R_e$", color="blue")
    ax2.set_xlim(0.001, 1)
    ax2.set_ylim(0, 4)
    ax2.set_xlabel(r"$\frac{V_E}{V_c}$" + "[-]")
    ax2.set_ylabel(r"$\frac{R}{R_E}$" + "[-]")
    ax2.set_title("Range flight ratio during gliding re-entry")
    ax2.grid(True)

    ax3.plot(V_Vc_ratio, a_weight_amax_ratio, label=r"$a/a_max$", color="blue")
    ax3.plot(V_Vc_ratio, a_no_weight_amax_ratio, label=r"$a/a_max$", color="blue", linestyle="dotted")
    ax3.set_xlim(0.001, 1)
    ax3.set_ylim(0.001, 1)
    ax3.set_xlabel(r"$\frac{V}{V_c}$" + "[-]")
    ax3.set_ylabel(r"$\frac{a}{a_max}$" + "[-]")
    ax3.set_title("Acceleration ratio during gliding re-entry")
    ax3.grid(True)

    plt.close(fig)  # Close the plot to prevent default PNG rendering

    # Change plotting style
    change_plot_style()

    # Display the SVG in Marimo using mo.Html
    svg_data = convert_fig_to_svg(fig)
    plot_svg(svg_data)
    return


@app.cell
def _(get_velocity_ratio, get_flight_path_angle, get_flight_range, get_flight_time, get_max_decelaration,
      lift_parameter, np, scale_height, altitude, lift_drag_ratio, entry_circular_ratio):
    beta = 1 / scale_height.value

    # Constants
    g = 9.81  # m/s^2
    rho = 1.225  # kg/m^3
    Re = 6378000

    # Circular velocity
    Vc = np.sqrt(g * Re)
    entry_circular_ratio = np.linspace(0, 1, 100)

    # velocity ratio

    V_Vc_ratio, height = get_velocity_ratio(lift_parameter.value, beta, altitude.value, Vc, rho)

    # equilibrium flight path
    flight_path_eq = get_flight_path_angle(V_Vc_ratio, beta, lift_drag_ratio.value)  # TODO: add to table

    flight_range_ratio = get_flight_range(lift_drag_ratio.value, entry_circular_ratio)

    flight_duration = get_flight_time(Vc, lift_drag_ratio.value, entry_circular_ratio, g)  # TODO: add to table

    decelaration_weight, decelaration_no_weight, max_deceleration_weight, max_deceleration_no_weight = get_max_decelaration(
        g, lift_drag_ratio.value, beta, Re, V_Vc_ratio)
    a_weight_amax_ratio = decelaration_weight / max_deceleration_weight
    a_no_weight_amax_ratio = decelaration_no_weight / max_deceleration_no_weight

    return V_Vc_ratio, height, flight_path_eq, flight_range_ratio, entry_circular_ratio, flight_duration, a_weight_amax_ratio, a_no_weight_amax_ratio


@app.cell
def _(mo):
    # sliders
    lift_parameter = mo.ui.slider(0, 20000, 100, value=2000)
    lift_drag_ratio = mo.ui.slider(0, 10, 1, value=2)
    scale_height = mo.ui.slider(0, 10000, 100, value=7200)
    altitude = mo.ui.slider(0, 120, 10, value=90)
    return altitude, lift_drag_ratio, lift_parameter, scale_height


@app.cell
def _(np):
    def get_velocity_ratio(
            lift_parameter: float,
            beta: float,
            altitude: float,
            Vc: float,
            rho: float
    ) -> float:
        """
        Calculate the velocity ratio for gliding re-entry

        Parameters:
        ----------
        lift_parameter : float
            Lift parameter of the spacecraft.
        altitude : float
            Altitude of spacecraft.
        scale_height : float
            Scale height of the atmosphere.

        Returns:
        -------
        float
            Velocity ratio.
        """

        # Altitude
        h = np.linspace(0, altitude * 1000, 10000)

        # Calculate normalised velocity
        V_Vc_ratio = np.sqrt((lift_parameter) / (0.5 * rho * np.exp(-beta * h) * Vc ** 2 + lift_parameter))

        return V_Vc_ratio, h

    def get_flight_path_angle(
            V_Vc_ratio: float,
            beta: float,
            lift_drag_ratio: float
    ) -> float:
        """
        Calculate the velocity ratio for gliding re-entry

        Parameters:
        ----------
        V_Vc_ratio : float
            Contribution of circular velocity
        beta : float
            1 / height scale
        lift_drag_ratio : float
            The L/D ratio of the spacecraft

        Returns:
        -------
        float
             Equivalent flight path angle.
        """
        # Constants
        Re = 6378000

        flight_path_angle = -(1 / (beta * Re)) * (2 / lift_drag_ratio) * V_Vc_ratio ** (-2)

        return flight_path_angle

    def get_flight_range(lift_drag_ratio,
                         entry_circular_velocity_ratio):
        """
        Calculate the ratio between flight range and earth radius.

        Parameters:
        ----------
        lift_drag_ratio : float
            The L/D ratio of the spacecraft

        entry_circular_velocity_ratio : float
            The ratio between entry velocity and circular

        Returns:
        -------
        float
             Equivalent flight path angle.
        """
        Rf_Re = -0.5 * lift_drag_ratio * np.log(1 - (entry_circular_velocity_ratio ** 2))

        return Rf_Re

    def get_flight_time(Vc,
                        lift_drag_ratio,
                        entry_circular_ratio,
                        g):
        t_flight = ((0.5 * Vc * lift_drag_ratio) / (2 * g)) * np.log(
            (1 + entry_circular_ratio) / (1 - entry_circular_ratio))

        return t_flight

    def get_max_decelaration(g,
                             lift_drag_ratio,
                             beta,
                             Re,
                             v_vc_ratio):
        a_weight = g * lift_drag_ratio ** (-1) * (1 - (v_vc_ratio ** 2) - (2 / (beta * Re * v_vc_ratio ** 2)))
        a_no_weight = g * lift_drag_ratio ** (-1) * (1 - (v_vc_ratio ** 2))
        a_max_weight = g * lift_drag_ratio ** (-1) * (1 - 2 * np.sqrt(2 / (beta * Re)))
        a_max_no_weight = g * lift_drag_ratio ** (-1)
        return a_weight, a_no_weight, a_max_weight, a_max_no_weight

    def get_stagnation_heatflux(nose_radius,
                                lift_parameter,
                                Vc,
                                V_Vc_ratio,
                                n,
                                m,
                                rho):
        cstar = 1.1097e8
        c1 = cstar * (1 / np.sqrt(rho)) * (1 / Vc ** 3)
        c2 = c1 * (1 / nose_radius ** n) * (2 * lift_parameter / Vc ** 2) ** (1 - n) * Vc ** m

        qc = c2 * (((1 / V_Vc_ratio ** 2) - 1) ** (1 - n)) * V_Vc_ratio ** m
        return qc

    return (get_velocity_ratio, get_flight_path_angle, get_flight_range, get_flight_time, get_max_decelaration)


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
