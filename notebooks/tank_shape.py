import marimo

__generated_with = "0.13.6"
app = marimo.App(width="medium")


@app.cell(hide_code=True)
def _(
    cylinder_length,
    mo,
    propellant_density,
    propellant_mass,
    propellant_volume,
    tank_diameter,
    tank_length,
):
    mo.md(
        rf"""
    ## Tank Shaping

    | Parameter | Adjustment | Value |
    |:---|:---:|:---|
    | Propellant Mass | {propellant_mass} | {propellant_mass.value} kg |
    | Propellant Density | {propellant_density} | {propellant_density.value} kg/m3 |
    | Propellant Volume | - | {round(propellant_volume, 1)} m3 |
    | Tank Diameter | {tank_diameter} | {tank_diameter.value} m |
    | Tank Length | - | {round(tank_length, 1)} m |
    | Cylinder Length | - | {round(cylinder_length, 1)} m |
    """
    )
    return


@app.cell(hide_code=True)
def _(
    change_plot_style,
    convert_fig_to_svg,
    cylinder_length,
    np,
    plot_svg,
    plt,
    tank_diameter,
):
    # the tank is plotted as two vertical lines and two semicircles.
    fig, ax = plt.subplots()
    ax.plot([0, 0], [0 + tank_diameter.value / 2, cylinder_length + tank_diameter.value / 2], color="blue", lw=2)  # Left side
    ax.plot([tank_diameter.value, tank_diameter.value], [0 + tank_diameter.value / 2, cylinder_length + tank_diameter.value / 2], color="blue", lw=2)  # Right side

    # Semi-circular dome on top
    theta = np.linspace(0, np.pi, 100)
    radius = tank_diameter.value / 2
    center_x = tank_diameter.value / 2

    # Corrected x-coordinates for the dome:
    # Original x: tank_diameter.value / 2 * np.cos(theta)
    # Corrected x: center_x + radius * np.cos(theta)

    # Semi-circular top
    ax.plot(
        center_x + radius * np.cos(theta),  # Centered horizontally
        cylinder_length + radius * np.sin(theta) + tank_diameter.value / 2, # Positioned on top of the rectangle
        color="blue",
        lw=2,
    )

    # Semi-circular bottom
    ax.plot(
        center_x + radius * np.cos(theta),  # Centered horizontally
        -radius * np.sin(theta) + tank_diameter.value / 2,  # Positioned at the bottom
        color="blue",
        lw=2,
    )

    ax.set_xlim(-1, 11)
    ax.set_ylim(-6, 60)
    ax.set_xlabel("x [m]")
    ax.set_ylabel("y [m]")
    ax.set_title("Tank shape")
    ax.set_aspect("equal")
    ax.grid(True)
    plt.close(fig)  # Close the plot to prevent default PNG rendering

    # Change plotting style
    change_plot_style()

    # Display the SVG in Marimo using mo.Html
    svg_data = convert_fig_to_svg(fig)
    plot_svg(svg_data)
    return


@app.cell(hide_code=True)
def _(
    allowable_stress,
    calculate_tank_mass,
    calculate_tank_thickness,
    cylinder_length,
    material_density,
    mo,
    propellant_pressure,
    tank_diameter,
    tank_length,
):
    mo.md(
        rf"""
    ## Tank Mass Estimates

    | Parameter | Adjustment | Value |
    |:---|:---:|:---|
    | Propellant Pressure | {propellant_pressure} | {round(propellant_pressure.value*1e-5, 1)} bar |
    | Allowable Stress | {allowable_stress} | {allowable_stress.value*1e-6} MPa |
    | Tank Thickness | - | {round(calculate_tank_thickness(tank_diameter, cylinder_length, propellant_pressure, allowable_stress.value, )*1e3, 1)} mm |
    | Material Density | {material_density} | {material_density.value} kg/m3 |
    | Tank Mass | - | {round(calculate_tank_mass(tank_diameter.value, tank_length, calculate_tank_thickness(tank_diameter, cylinder_length, propellant_pressure, allowable_stress.value), material_density.value), 1)} kg |

    Drop downs could be implemented here for different materials that have allowable stress and material density pre-defined.
    """
    )
    return


@app.cell(hide_code=True)
def _(mo):
    mo.md(
        r"""
    ### Tank Mass Assumptions

    Assuming a simple thin-walled pressure vessel where thickness is given by $t = \frac{P D_t}{2 \sigma}$, where $P$ is the internal pressure, $D_t$ is the diameter, and $\sigma$ is the allowable stress.
    """
    )
    return


@app.cell(hide_code=True)
def _(mo):
    # sliders
    propellant_mass = mo.ui.slider(0, 50000, 100, value=10000)
    propellant_density = mo.ui.slider(70, 77, 0.1, value=70)
    propellant_pressure = mo.ui.slider(0.1e5, 15e5, 0.1e5, value=5e5)
    tank_diameter = mo.ui.slider(3, 10, 0.1, value=5)
    allowable_stress = mo.ui.slider(270e6, 1400e6, 10e6, value=300e6)
    material_density = mo.ui.slider(2700, 8000, 10, value=7850)
    return (
        allowable_stress,
        material_density,
        propellant_density,
        propellant_mass,
        propellant_pressure,
        tank_diameter,
    )


@app.cell(hide_code=True)
def _(np, propellant_density, propellant_mass, tank_diameter):
    # calculating tank shape
    propellant_volume = propellant_mass.value / propellant_density.value

    # dome volume
    dome_volume = (4 / 3) * np.pi * (tank_diameter.value / 2) ** 3

    # cylinder length
    cylinder_length = (propellant_volume - dome_volume) / (
        np.pi * (tank_diameter.value / 2) ** 2
    )

    # total tank length
    tank_length = cylinder_length + tank_diameter.value

    # feasibility check
    if cylinder_length < 0:
        raise ValueError("Negative tank length. Restart kernel.")
    return cylinder_length, propellant_volume, tank_length


@app.cell
def _(np):
    # calculating tank thickness
    def calculate_tank_thickness(tank_diameter, tank_length, propellant_pressure, allowable_stress):
        # Assuming a thin-walled cylinder for simplicity
        # Using the formula: t = (P * D) / (2 * σ)
        # where P is the internal pressure, D is the diameter, and σ is the allowable stress
        thickness = (propellant_pressure.value * tank_diameter.value) / (2 * allowable_stress)
        return thickness

    def calculate_tank_mass(tank_diameter, tank_length, thickness, material_density):
        # Calculate the volume of the tank material (assuming a cylindrical shape)
        outer_radius = tank_diameter / 2 + thickness
        inner_radius = tank_diameter / 2
        volume = np.pi * (outer_radius**2 - inner_radius**2) * tank_length

        mass = volume * material_density
        return mass
    return calculate_tank_mass, calculate_tank_thickness


@app.cell(hide_code=True)
def _(mo):
    mo.md(
        r"""
    ### Tank shape calculations

    For rapid sizing, the tank is assumed to be cylindrical with two semi-circular domes. All dimensions are derived from the volume the tank needs to contain and setting a diameter. The below equations provide the relationships.

    | Parameter | Symbol |
    |:---|:---|
    | Propellant Volume | $V$ |
    | Tank Diameter | $D_t$ |
    | Tank Length | $L_t$ |
    | Cylinder Length | $L_c$ |
    | Dome Volume | $V_d$ |
    | Cylinder Volume | $V_c$ |

    The volume contained by the domes is:

    $$ V_d = \frac{4}{3} \pi \left(\frac{D_t}{2}\right)^3 $$

    The volume contained by the cylinder is:

    $$ V_c = \pi \left(\frac{D_t}{2}\right)^2 L_c $$

    Given the diameter, the cylinder length is calculated as:

    $$ L_c = \frac{V - V_d}{\pi \left(\frac{D_t}{2}\right)^2} $$
    """
    )
    return


@app.cell
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
    import numpy as np
    import matplotlib.pyplot as plt
    import io
    return io, mo, np, plt


if __name__ == "__main__":
    app.run()
