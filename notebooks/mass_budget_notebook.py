import marimo

__generated_with = "0.13.6"
app = marimo.App(width="medium")


@app.cell
def _(
    cryogenic,
    glom1,
    glom2,
    latitude,
    launch_site,
    mo,
    payload_mass,
    propellant,
):
    mo.md(
        rf"""
    # Mass Budget Tool

    Adjust the values below to perform the sizing.

    | Parameter | Selector | Value | Note |
    |:---|:---:|:---:|:---|
    | Launch Site | {launch_site} | {launch_site.value} | Select the launch site |
    | Propellant | {propellant} | {propellant.value} | Used by the entire vehicle |
    | Cryogenic | {cryogenic} | {cryogenic.value} | Is the propellant cryogenic? |
    | Latitude | - | {latitude(launch_site.value)} ° | Latitude of the launch site |
    | Payload Mass | {payload_mass} | {payload_mass.value} | LH2 mass |
    | Gross Lift-off Mass | - | {round(glom1(payload_mass.value) if payload_mass.value < 7500 else glom2(payload_mass.value))} | Calculated from payload mass |
    """
    )
    return


@app.cell
def _(mo):
    mass_fraction = mo.ui.slider(1, 20, 0.1)
    payload_mass = mo.ui.slider(0, 20000, 10)
    cryogenic = mo.ui.checkbox()
    propellant = mo.ui.dropdown(["LH2/LOX", "LCH4/LOX"], value="LH2/LOX")
    launch_site = mo.ui.dropdown(["Kourou", "Vandenberg", "Cape Canaveral"], value="Kourou")
    return cryogenic, launch_site, payload_mass, propellant


@app.cell(hide_code=True)
def _():
        # The below are functions according to the TU Delft model by B.T.C. Zandbergen

    def glom1(payload_mass_kg: float) -> float:
        """
        Calculates Gross Lift Off Mass (GLOM) in kg for a launch vehicle.
        This function is intended for payload mass into LEO of up to 7500 kg.

        Args:
            payload_mass_kg: The payload mass in kilograms.

        Returns:
            The calculated GLOM in kg.

        Formula:
            GLOM [ton] = 42.827 * M_payload_mass [ton] + 12.996
        """
        if payload_mass_kg < 0:
            raise ValueError("Payload mass cannot be negative.")

        payload_mass_ton = payload_mass_kg / 1000.0
        glom_ton = 42.827 * payload_mass_ton + 12.996
        return glom_ton*1000

    def glom2(payload_mass_kg: float) -> float:
        """
        Calculates Gross Lift Off Mass (GLOM) in kg for a launch vehicle.
        This function is intended for payload mass into LEO from 7500 kg up to 120,000 kg.

        Args:
            payload_mass_kg: The payload mass in kilograms.

        Returns:
            The calculated GLOM in kg.

        Formula:
            GLOM [ton] = 25.138 * M_payload_mass [ton] + 169.358
        """
        if payload_mass_kg < 0:
            raise ValueError("Payload mass cannot be negative.")
    
        payload_mass_ton = payload_mass_kg / 1000.0
        glom_ton = 25.138 * payload_mass_ton + 169.358
        return glom_ton*1000
    return glom1, glom2


@app.cell
def _(payload_mass_kg):
    # The below are functions used to calculate ΔV

    def DeltaV(orbit_height: float, orbit_inclination: float) -> float:
        """
        Calculates the ΔV for a launch vehicle with a payload mass
    
        Args:
            payload_mass_kg: The payload mass in kilograms.

        Returns:
            The calculated delta-v in m/s.
        """
        if payload_mass_kg < 0:
            raise ValueError("Payload mass cannot be negative.")

        payload_mass_ton = payload_mass_kg / 1000.0
        delta_v_m_s = 4.5 * payload_mass_ton + 2.5
        return delta_v_m_s


    def latitude(launch_site: str) -> float:
        """
        Returns the latitude of the launch site in degrees.
    
        Args:
            launch_site: The name of the launch site.

        Returns:
            The latitude of the launch site in degrees.
        """
        if launch_site == "Kourou":
            return 5.236
        elif launch_site == "Vandenberg":
            return 34.736
        elif launch_site == "Cape Canaveral":
            return 28.392
        else:
            raise ValueError("Invalid launch site.")
    return (latitude,)


@app.cell(hide_code=True)
def _(mo):
    mo.md(
        r"""
    # ΔV Tool

    Drag loss equation:

    $$\Delta V_d = \int_{t_0}^{t_f} \frac{D(t)}{m(t)} dt = \int_{t_0}^{t_f} \frac{\frac{1}{2} \rho(h(t)) v(t)^2 C_D(M(t), \alpha(t)) A}{m(t)} dt$$

    Gravity loss equation:

    $$\Delta V_g = \int_{t_0}^{t_f} g(h(t)) \sin(\gamma(t)) dt$$

    Steering loss equation (from Design of Rockets and Space Launch Vehicles, Chapter 6):

    $$\Delta V_s = TBD$$

    Tsiolkovski rocket equation:

    $$\Delta V = I_{sp} g_0 \ln \left( \frac{m_0}{m_f} \right) - \Delta V_d - \Delta V_g - \Delta V_s$$

    **This is too complicated for the DSE scope.**
    """
    )
    return


@app.cell(hide_code=True)
def _(mo):
    mo.md(
        r"""
    ## Gravity Turn

    In the simplied gravity turn, the gravity essentially pulls the rocket down, and the rocket follows a parabolic trajectory. The angle of the rocket with respect to the vertical is given by:

    $$\frac{d\gamma}{dt} = -\frac{g \cos \gamma}{v}$$

    In the cell below, a simple simulation of the gravity turn is performed. The simulation is performed in a 2D plane, and the rocket is assumed to be a point mass.

    The initial assumptions are:

    - The rocket is launched at an angle of 0 degrees with respect to the vertical, $\gamma_0 = 0$
    - The rocket is launched at a velocity of 0 m/s, $V_0 = 0$
    - The rocket is launched at a time of 0 seconds, $t_0 = 0$
    - The rocket is launched at an altitude of 0 m, $h_0 = 0$
    - The rocket is launched at a time step of 0.1 seconds, $\Delta t = 0.1$

    Velocity is a function of the thrust to weight ratio.
    """
    )
    return


@app.cell
def _(change_plot_style, convert_fig_to_svg, np, plot_svg, plt):
    # --- Initial values and Parameters ---
    # initial_pitchover_from_vertical: Small angle rocket pitches away from pure vertical
    # to initiate the gravity turn.
    initial_pitchover_from_vertical = 1e-3  # radians (approx 0.057 degrees)

    gamma_0 = np.pi/2 - initial_pitchover_from_vertical  # Initial flight path angle (from horizontal)
    V_0 = 10.0  # Initial velocity in m/s (non-zero, after initial vertical ascent)
    t_0 = 0.0  # Initial time in seconds
    h_0 = 50.0 # Initial altitude in meters (after initial vertical ascent)

    delta_t = 0.1  # Time step in seconds
    g = 9.81  # Gravitational acceleration in m/s^2
    TWR = 1.5  # Thrust to initial Weight Ratio (assumed constant T/m = TWR*g)

    # --- New Parameters for User Requests ---
    # K_resistance: Factor to slow down the turn.
    # 0 < K_resistance <= 1. (1.0 means no artificial resistance)
    K_resistance = 0.8

    # target_gamma_horizontal: The flight path angle (from horizontal) at which the turn stops.
    # 0.0 means perfectly horizontal flight.
    target_gamma_horizontal = 0.0

    # --- Simulation Time ---
    t_end = 400.0  # End time for simulation in seconds

    # --- Arrays to store results ---
    t = np.arange(t_0, t_end + delta_t, delta_t) # Ensure t_end is included
    V = np.zeros(len(t))
    h = np.zeros(len(t))
    gamma = np.zeros(len(t)) # Flight path angle from horizontal, in radians

    # --- Set initial conditions ---
    V[0] = V_0
    gamma[0] = gamma_0
    h[0] = h_0

    # --- Time stepping simulation ---
    for i in range(1, len(t)):
        # Previous step values
        V_prev = V[i-1]
        gamma_prev = gamma[i-1]
        h_prev = h[i-1]

        # 1. Calculate acceleration and update velocity
        # dV/dt = Thrust_acceleration - g * sin(gamma)
        # Assuming Thrust_acceleration = TWR * g (i.e., m is constant or TWR is T/(m*g) always)
        thrust_accel = TWR * g
        dV_dt = thrust_accel - g * np.sin(gamma_prev)
        V[i] = V_prev + dV_dt * delta_t

        # Ensure velocity doesn't become negative (shouldn't happen if TWR > 1)
        if V[i] < 0:
            V[i] = 0.0

        # 2. Update altitude
        # dh/dt = V * sin(gamma)
        # Using V_prev for consistency in Euler integration (state at start of interval)
        dh_dt = V_prev * np.sin(gamma_prev)
        h[i] = h_prev + dh_dt * delta_t

        # 3. Update flight path angle (gamma from horizontal)
        # Check if the turn should continue
        if gamma_prev <= target_gamma_horizontal:
            # Already at or below the target horizontal angle, so stop turning
            gamma[i] = target_gamma_horizontal
        elif V[i] == 0: # Velocity is zero, cannot calculate turn rate (should generally not happen)
            gamma[i] = gamma_prev # No change in angle
        else:
            # Natural rate of change of gamma for a gravity turn
            # d(gamma)/dt = - (g * cos(gamma)) / V
            # Using V[i] (velocity at the end of the step) as in your original code for the denominator
            d_gamma_dt_natural = - (g * np.cos(gamma_prev) / V[i])

            # Apply resistance factor
            actual_d_gamma_dt = K_resistance * d_gamma_dt_natural

            gamma[i] = gamma_prev + actual_d_gamma_dt * delta_t

            # Ensure gamma does not "overshoot" the target angle due to discrete time step
            # (Since gamma is decreasing towards target_gamma_horizontal)
            if gamma[i] < target_gamma_horizontal:
                gamma[i] = target_gamma_horizontal

    # --- Output a few values as an example (optional) ---
    # print("Time (s), Velocity (m/s), Altitude (m), Gamma (deg)")
    # for j in range(0, len(t), int(len(t)/10)): # Print 10 points
    #     print(f"{t[j]:.1f}, {V[j]:.2f}, {h[j]:.2f}, {np.degrees(gamma[j]):.2f}")



    # plot
    fig, ax = plt.subplots()
    ax.plot(t, gamma*180/np.pi, label="Altitude [m]")
    # ax.plot(t, V, label="Velocity [m/s]")
    # ax.set_xlim(1, 3)
    # ax.set_ylim(1, 10)
    ax.set_xlabel(r"$t$" + r" [s]")
    ax.set_ylabel(r"$\gamma$" + " [°]")
    ax.set_title("Simplified gravity turn simulation")
    ax.grid(True)
    ax.legend()
    plt.close(fig)  # Close the plot to prevent default PNG rendering

    # Change plotting style
    change_plot_style()

    # Display the SVG in Marimo using mo.Html
    svg_data = convert_fig_to_svg(fig)
    plot_svg(svg_data)

    return


@app.cell(hide_code=True)
def _(mo):
    mo.md(r"""# Helpers""")
    return


@app.cell
def _(io, mo, plt):
    def change_plot_style():
        plt.rcParams.update({
            "text.usetex": True,
            "font.family": "Computer Modern Roman"
        })

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
    import matplotlib.pyplot as plt
    import numpy as np
    import io
    return io, mo, np, plt


if __name__ == "__main__":
    app.run()
