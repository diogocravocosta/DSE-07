import marimo

__generated_with = "0.13.8"
app = marimo.App(width="full")


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

    ### Equations of Motion

    The following three equations of motion are applicable to the ballistic re-entry problem:

    $$ m \frac{dV}{dt} = -D - mg \sin \gamma $$

    $$ mV \frac{d\gamma}{dt} = L - mg \cos \gamma \left(1 - \frac{V^2}{V_c^2}\right) $$

    $$ \frac{dR}{dt} = \frac{dh}{dt} = V \sin \gamma $$

    The third equation is a kinematic relation.

    It is assumed that $\frac{\text{d} \gamma}{\text{d} h} \approx 0$, flight path angle does not change with altitude.

    ### Normalised Velocity

    The analytical solution for the normalised velocity of the ballistic re-entry problem is given by:

    $$ \frac{V}{V_E} = \exp \left(\frac{1}{2} \frac{g \rho}{K \beta \sin \gamma_E}\right)$$

    $K$ is the weight referenced ballistic coefficient and $\beta$ represents the exponential atmosphere density model with constant temperature. It can also be written as $\beta = \frac{1}{H_s}$ where $H_s$ is the scale height. It follows from the exponential atmosphere model that $\rho = \rho_0 e^{-\beta h}$. This can be substituted into the equation to solve for normalised velocity as a function of altitude. The relationship is plotted below.

    ### Deceleration

    Deceleration as a function of velocity is given by:

    $$ \bar{a} = -\frac{\text{d} V}{\text{d} t} = \beta \sin \gamma_E V_E^2 \left( \frac{V}{V_E} \right) ^2 \ln \frac{V}{V_E} $$

    Maximum deceleration is given by:

    $$ \bar{a}_{\text{max}} = - \left( \frac{dV}{dt} \right)_{\text{max}} = - \frac{\beta \sin \gamma_E}{2e} V_E^2 $$

    The altitude at which maximum deceleration occurs is given by:

    $$ h' = \frac{1}{\beta} \ln \left( - \frac{\rho_0 g}{K \beta \sin \gamma_E} \right) $$

    ### Thermal Loads

    The normalized heat flux is given by:

    $$ \frac{q_c}{q_{c,max}} = \left(\frac{3e}{1-n}\right)^{1-n} \left|\ln{\frac{V}{V_E}}\right|^{1-n} \left(\frac{V}{V_E}\right)^3 $$

    where $n$ is a parameter for either:

    - laminar boundary layer, $n=0.2$
    - turbulent boundary layer, $n=0.5$

    The maximum heat flux can be very roughly estimated using the following relationship:

    $$ q_{c,max} = c_1 V_E^3 \sqrt{-\frac{1}{3e} \frac{1}{R_N} \frac{K \beta}{g} \sin \gamma_E} $$

    where $c_1 = c^* \frac{1}{R_N^n} \left( \frac{\rho}{\rho_0} \right)^{1-n} \left( \frac{V}{V_c} \right)^m$ with $c^* = 1.1097\cdot 10^8 \sqrt{m}$ (constant), $R_N$ as the nose radius, $V_c$ as the circular velocity ($V_c = \sqrt{g R} \approx 7920$), and $m = 3$ as a costant specific to Earth.

    The altitude, at which maximum heat flux occurs, is given by:

    $$ h'' - h' = - \frac{\ln{\frac{2}{3} (1-n)}}{\beta} $$

    where $h'$ is the altitude at which maximum deceleration occurs and $h''$ is the altitude at which maximum heat flux occurs at the stagnation point. Other arbitrary points on the vehicle are not considered in this analysis.

    ### Ballistic parameter

    There are two different definitions of the ballistic parameter:

    - $K = \frac{m g}{C_D S}$, _weight referenced_ ballistic parameter
    - $K_m = \frac{m}{C_D S}$, _mass refererenced_ ballistic parameter
    """
    )
    return


@app.cell
def _(drag_coefficient, mo, vehicle_diameter, vehicle_mass):
    mo.md(
        rf"""
    The below table provides a rapid way to estimate the ballistic parameter. The reference area and drag coefficient are assumed for a capsule-like vehicle, so the reference area should correspond to circular area at the bottom of a capsule.

    | Variable | Adjustment | Value |
    |:---|:---:|:---|
    | Vehicle Mass | {vehicle_mass} | {vehicle_mass.value} kg |
    | Drag Coefficient | {drag_coefficient} | {drag_coefficient.value} |
    | Vehicle Diameter | {vehicle_diameter} | {vehicle_diameter.value} m |
    | Reference Area | - | {round(3.14159*vehicle_diameter.value**2/4, 1)} m^2 |
    | Ballistic Parameter | - | {round(vehicle_mass.value * 9.81 / (drag_coefficient.value * 3.14159*vehicle_diameter.value**2/4), 1)} N/m^2 |
    """
    )
    return


@app.cell
def _(
    maximum_ballistic_deceleration,
    maximum_ballistic_deceleration_altitude,
    maximum_ballistic_heat_flux,
    maximum_ballistic_heat_flux_altitude,
    mo,
):
    mo.md(
        f"""
    ## Results

    | Variable | Adjustment | Value |
    |:---|:---:|:---|
    | Maximum deceleration | - | {round(maximum_ballistic_deceleration/9.81, 1)} g |
    | Maximum deceleration altitude | - | {round(maximum_ballistic_deceleration_altitude/1000, 1)} km |
    | Maximum heat flux | - | {round(maximum_ballistic_heat_flux*1e-3, 1)} kW/m^2 |
    | Maximum heat flux altitude | - | {round(maximum_ballistic_heat_flux_altitude/1000, 1)} km |
    """
    )
    return


@app.cell
def _(
    ballistic_parameter,
    boundary_layer,
    entry_flight_path_angle,
    entry_speed,
    mo,
    nose_radius,
    scale_height,
):
    mo.md(
        f"""
    ## Ballistic Entry Variables

    | Variable | Adjustment | Value |
    |:---|:---:|:---|
    | Ballistic Parameter | {ballistic_parameter} | {ballistic_parameter.value} N/m^2 |
    | Entry Speed | {entry_speed} | {entry_speed.value} m/s |
    | Entry Flight Path Angle | {entry_flight_path_angle} | {entry_flight_path_angle.value} ° |
    | Scale Height | {scale_height} | {scale_height.value} m |
    | Boundary Layer | {boundary_layer} | {boundary_layer.value} |
    | Nose radius | {nose_radius} | {nose_radius.value} m |
    """
    )
    return


@app.cell
def _(
    altitude,
    change_plot_style,
    convert_fig_to_svg,
    deceleration,
    normalised_heat_flux,
    normalised_velocity,
    plot_svg,
    plt,
):
    # the tank is plotted as two vertical lines and two semicircles.
    fig, (ax1, ax2, ax3) = plt.subplots(3, 1, figsize=(10, 10))

    ax1.plot(normalised_velocity, altitude/1000, label=r"$V/V_E$", color="blue")
    ax1.set_xlim(0.001, 1)
    ax1.set_ylim(0, 80)
    #ax1.set_xlabel(r"$\frac{V}{V_E}$" + "[-]")
    ax1.set_ylabel(r"$h$ [km]")
    ax1.set_title("Normalized velocity during ballistic re-entry")
    ax1.grid(True)

    ax2.plot(normalised_velocity, deceleration, label=r"$a$", color="red")
    ax2.set_xlim(0.001, 1)
    ax2.set_ylim(0)
    #ax2.set_xlabel(r"$\frac{V}{V_E}$" + "[-]")
    ax2.set_ylabel(r"$a$ [m/s$^2$]")
    ax2.set_title("Deceleration during ballistic re-entry")
    ax2.grid(True)

    ax3.plot(normalised_velocity, normalised_heat_flux, label="$q_c/q_{c,\text{max}}$", color="green")
    ax3.set_xlim(0.001, 1)
    ax3.set_ylim(0)
    ax3.set_xlabel(r"$\frac{V}{V_E}$" + "[-]")
    ax3.set_ylabel(r"$\frac{q_c}{q_{c,\text{max}}}$ [-]")
    ax3.set_title("Normalized heat flux during ballistic re-entry")
    ax3.grid(True)

    plt.close(fig)  # Close the plot to prevent default PNG rendering

    # Change plotting style
    change_plot_style()

    # Display the SVG in Marimo using mo.Html
    svg_data = convert_fig_to_svg(fig)
    plot_svg(svg_data)
    return


@app.cell
def _(
    ballistic_parameter,
    boundary_layer,
    entry_flight_path_angle,
    entry_speed,
    get_ballistic_deceleration,
    get_ballistic_normalised_heat_flux,
    get_ballistic_normalised_velocity,
    get_max_ballistic_deceleration,
    get_max_ballistic_deceleration_altitude,
    get_max_ballistic_heat_flux,
    get_max_ballistic_heat_flux_altitude,
    nose_radius,
    np,
    scale_height,
):
    altitude = np.linspace(0, 80e3, 100)  # altitude from 0 to 80 km

    normalised_velocity = get_ballistic_normalised_velocity(
        altitude,
        ballistic_parameter=ballistic_parameter.value,
        entry_speed=entry_speed.value,
        entry_flight_path_angle=entry_flight_path_angle.value,
        scale_height=scale_height.value,
    )

    deceleration = get_ballistic_deceleration(
        normalised_velocity,
        entry_speed=entry_speed.value,
        entry_flight_path_angle=entry_flight_path_angle.value,
        scale_height=scale_height.value,
    )

    normalised_heat_flux = get_ballistic_normalised_heat_flux(
        normalised_velocity,
        entry_speed=entry_speed.value,
        entry_flight_path_angle=entry_flight_path_angle.value,
        scale_height=scale_height.value,
        n=boundary_layer.value,
    )

    # constants
    cstar = 1.1097e8
    rho_0 = 1.225

    # maximum heat flux calculations
    V_c = 7920  # m/s, usual value taken for Earth, can be calculated as sqrt(g*R)
    c1 = cstar * (1 / np.sqrt(rho_0)) * (1 / V_c**3)
    maximum_ballistic_heat_flux = get_max_ballistic_heat_flux(
        ballistic_parameter=ballistic_parameter.value,
        entry_speed=entry_speed.value,
        entry_flight_path_angle=entry_flight_path_angle.value,
        scale_height=scale_height.value,
        nose_radius=nose_radius.value,
        c1=c1
    )

    maximum_ballistic_heat_flux_altitude = get_max_ballistic_heat_flux_altitude(
        ballistic_parameter=ballistic_parameter.value,
        entry_flight_path_angle=entry_flight_path_angle.value,
        scale_height=scale_height.value,
        n=boundary_layer.value,
    )

    # maximum deceleration calculations
    maximum_ballistic_deceleration = get_max_ballistic_deceleration(
        entry_speed=entry_speed.value,
        entry_flight_path_angle=entry_flight_path_angle.value,
        scale_height=scale_height.value,
    )

    maximum_ballistic_deceleration_altitude = get_max_ballistic_deceleration_altitude(
        ballistic_parameter=ballistic_parameter.value,
        entry_flight_path_angle=entry_flight_path_angle.value,
        scale_height=scale_height.value,
    )
    return (
        altitude,
        deceleration,
        maximum_ballistic_deceleration,
        maximum_ballistic_deceleration_altitude,
        maximum_ballistic_heat_flux,
        maximum_ballistic_heat_flux_altitude,
        normalised_heat_flux,
        normalised_velocity,
    )


@app.cell
def _(mo):
    # sliders
    ballistic_parameter = mo.ui.slider(0, 20000, 100, value=5000)
    entry_speed = mo.ui.slider(6e3, 10e3, 100, value=8e3)
    entry_flight_path_angle = mo.ui.slider(-90, 0, -0.1, value=-10)
    scale_height = mo.ui.slider(6000, 8000, 10, value=7200)
    nose_radius = mo.ui.slider(3, 10, 0.1, value=7)

    vehicle_mass = mo.ui.slider(10e3, 50e3, 1e3, value=20e3)
    drag_coefficient = mo.ui.slider(0.1, 2, 0.05, value=0.5)
    vehicle_diameter = mo.ui.slider(3, 10, 0.1, value=7)

    # dropdowns
    boundary_layer = mo.ui.dropdown({
        "laminar": 0.5,
        "turbulent": 0.2
    }, value="laminar")
    return (
        ballistic_parameter,
        boundary_layer,
        drag_coefficient,
        entry_flight_path_angle,
        entry_speed,
        nose_radius,
        scale_height,
        vehicle_diameter,
        vehicle_mass,
    )


@app.cell
def _(mo):
    mo.md(r"""# Functions""")
    return


@app.cell
def _(np):
    def get_ballistic_normalised_velocity(
        altitude: np.ndarray,
        ballistic_parameter: float,
        entry_speed: float,
        entry_flight_path_angle: float,
        scale_height: float,
    ) -> float:
        """
        Calculate the normalised velocity for a ballistic entry.

        Parameters:
        ----------
        ballistic_parameter : float
            Ballistic parameter of the spacecraft.
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
        rho_0 = 1.225  # kg/m^3

        # Exponential atmosphere model
        beta = 1 / scale_height

        # Calculate normalised velocity
        normalised_velocity = np.exp(
            (0.5 * g * rho_0 * np.exp(-beta*altitude)) / (ballistic_parameter * beta * np.sin(np.radians(entry_flight_path_angle)))
        )

        return normalised_velocity

    def get_ballistic_deceleration(
        normalised_velocity: np.ndarray,
        entry_speed: float,
        entry_flight_path_angle: float,
        scale_height: float,
    ) -> float:
        """
        Calculate the deceleration for a ballistic entry.

        Parameters:
        ----------
        normalised_velocity : np.ndarray
            Normalised velocity of the spacecraft.
        entry_speed : float
            Entry speed of the spacecraft.
        entry_flight_path_angle : float
            Entry flight path angle of the spacecraft.
        scale_height : float
            Scale height of the atmosphere.

        Returns:
        -------
        float
            Deceleration.
        """
        # Constants
        g = 9.81  # m/s^2
        rho_0 = 1.225  # kg/m^3

        # Exponential atmosphere model
        beta = 1 / scale_height

        # Calculate deceleration
        deceleration = beta * np.sin(np.radians(entry_flight_path_angle)) * (entry_speed**2) * (normalised_velocity**2) * np.log(normalised_velocity)

        return deceleration

    def get_ballistic_normalised_heat_flux(
        normalised_velocity: np.ndarray,
        entry_speed: float,
        entry_flight_path_angle: float,
        scale_height: float,
        n: float,
    ) -> float:
        """
        Calculate the normalised heat flux for a ballistic entry.

        Parameters:
        ----------
        normalised_velocity : np.ndarray
            Normalised velocity of the spacecraft.
        entry_speed : float
            Entry speed of the spacecraft.
        entry_flight_path_angle : float
            Entry flight path angle of the spacecraft.
        scale_height : float
            Scale height of the atmosphere.
        n : float
            Parameter for laminar or turbulent boundary layer (0.2 for laminar, 0.5 for turbulent).

        Returns:
        -------
        float
            Normalised heat flux.
        """
        # Exponential atmosphere model
        beta = 1 / scale_height

        # Calculate normalised heat flux
        normalised_heat_flux = ((3 * np.e) / (1 - n))**(1 - n) * np.abs(np.log(normalised_velocity))**(1 - n) * (normalised_velocity**3)

        return normalised_heat_flux

    def get_max_ballistic_heat_flux(
        ballistic_parameter: float,
        entry_speed: float,
        entry_flight_path_angle: float,
        scale_height: float,
        nose_radius: float,
        c1: float,
    ) -> float:
        """
        Calculate the maximum heat flux for a ballistic entry.

        Parameters:
        ----------
        ballistic_parameter : float
            Ballistic parameter of the spacecraft.
        entry_speed : float
            Entry speed of the spacecraft.
        entry_flight_path_angle : float
            Entry flight path angle of the spacecraft.
        scale_height : float
            Scale height of the atmosphere.
        nose_radius : float
            Nose radius of the spacecraft.
        c1 : float
            Constant for maximum heat flux calculation.

        Returns:
        -------
        float
            Maximum heat flux.
        """
        # Constants
        g = 9.81  # m/s^2

        # Exponential atmosphere model
        beta = 1 / scale_height

        # Calculate maximum heat flux
        max_heat_flux = c1 * entry_speed**3 * np.sqrt(-1 / (3 * np.e) * (1 / nose_radius) * (ballistic_parameter * beta / g) * np.sin(np.radians(entry_flight_path_angle)))

        return max_heat_flux

    def get_max_ballistic_deceleration(
        entry_speed: float,
        entry_flight_path_angle: float,
        scale_height: float,
    ) -> float:
        """
        Calculate the maximum deceleration for a ballistic entry.

        Parameters:
        ----------
        entry_speed : float
            Entry speed of the spacecraft.
        entry_flight_path_angle : float
            Entry flight path angle of the spacecraft.
        scale_height : float
            Scale height of the atmosphere.

        Returns:
        -------
        float
            Maximum deceleration.
        """
        # Exponential atmosphere model
        beta = 1 / scale_height

        # Calculate maximum deceleration
        max_deceleration = - (beta * np.sin(np.radians(entry_flight_path_angle)) / (2 * np.e)) * (entry_speed**2)

        return max_deceleration

    def get_max_ballistic_deceleration_altitude(
        ballistic_parameter: float,
        entry_flight_path_angle: float,
        scale_height: float,
    ) -> float:
        """
        Calculate the altitude at which maximum deceleration occurs for a ballistic entry.

        Parameters:
        ----------
        ballistic_parameter : float
            Ballistic parameter of the spacecraft.
        entry_flight_path_angle : float
            Entry flight path angle of the spacecraft.
        scale_height : float
            Scale height of the atmosphere.

        Returns:
        -------
        float
            Altitude at which maximum deceleration occurs.
        """
        # Constants
        g = 9.81  # m/s^2
        rho_0 = 1.225  # kg/m^3

        # Exponential atmosphere model
        beta = 1 / scale_height

        # Calculate altitude at maximum deceleration
        max_deceleration_altitude = (1 / beta) * np.log(- (rho_0 * g) / (ballistic_parameter * beta * np.sin(np.radians(entry_flight_path_angle)) * g))

        return max_deceleration_altitude

    def get_max_ballistic_heat_flux_altitude(
        ballistic_parameter: float,
        entry_flight_path_angle: float,
        scale_height: float,
        n: float,
    ) -> float:
        """
        Calculate the altitude at which maximum heat flux occurs for a ballistic entry.

        Parameters:
        ----------
        ballistic_parameter : float
            Ballistic parameter of the spacecraft.
        entry_flight_path_angle : float
            Entry flight path angle of the spacecraft.
        scale_height : float
            Scale height of the atmosphere.
        n : float
            Parameter for laminar or turbulent boundary layer (0.2 for laminar, 0.5 for turbulent).

        Returns:
        -------
        float
            Altitude at which maximum heat flux occurs.
        """
        # Exponential atmosphere model
        beta = 1 / scale_height

        # Calculate the maximum ballistic deceleration altitude
        max_ballistic_deceleration_altitude = get_max_ballistic_deceleration_altitude(
            ballistic_parameter=ballistic_parameter,
            entry_flight_path_angle=entry_flight_path_angle,
            scale_height=scale_height,
        )

        # Calculate altitude at maximum heat flux
        max_heat_flux_altitude = max_ballistic_deceleration_altitude - (np.log(2 / 3 * (1 - n))) / beta

        return max_heat_flux_altitude
    return (
        get_ballistic_deceleration,
        get_ballistic_normalised_heat_flux,
        get_ballistic_normalised_velocity,
        get_max_ballistic_deceleration,
        get_max_ballistic_deceleration_altitude,
        get_max_ballistic_heat_flux,
        get_max_ballistic_heat_flux_altitude,
    )


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
