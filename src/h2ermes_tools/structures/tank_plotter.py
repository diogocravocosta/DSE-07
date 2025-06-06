import matplotlib.pyplot as plt
import numpy as np
from h2ermes_tools.structures import tank_sizing


def draw_hemisphere(ax, center_y, radius, lower=True, color="k", label=None):
    theta = np.linspace(0, np.pi, 100)
    x = radius * np.cos(theta)
    y_offset = (radius / 2) * np.sin(theta)  # vertical radius = R/2 for 2:1 ellipse
    if lower:
        y = center_y - y_offset
    else:
        y = center_y + y_offset
    ax.plot(x, y, color, label=label)
    ax.plot(-x, y, color)


def plot_stacked_tanks(top_tank: tuple[float, float, float],
                       bottom_tank: tuple[float, float, float],
                       lh2_on_top: bool,
                       title: str="Stacked LH2 + LOX Tank Profile"):
    """
    Plot two tanks stacked vertically,
    with the top tank being either LH2 or LOX based on lh2_on_top.
    The tanks are plotted with y as the axis of symmetry (vertical),
    x as radius (horizontal).
    """
    # Unpack tank geometry
    height_top_tank, bottom_radius_top_tank, top_radius_top_tank  = top_tank
    height_bottom_tank, bottom_radius_bottom_tank, top_radius_bottom_tank = bottom_tank

    fig, ax = plt.subplots(figsize=(6, 16))

    # bottom dome (blue)
    draw_hemisphere(ax,
                    center_y=bottom_radius_bottom_tank,
                    radius=bottom_radius_bottom_tank,
                    lower=True,
                    color="blue",
                    label="LH2 bottom dome")

    # lower frustum (green)
    y_start = bottom_radius_bottom_tank
    y_end = bottom_radius_bottom_tank + height_bottom_tank
    ax.plot(
        [-bottom_radius_bottom_tank, -top_radius_bottom_tank],
        [y_start, y_end],
        color="green",
        label="LH2 frustum")
    ax.plot(
        [bottom_radius_bottom_tank, top_radius_bottom_tank],
        [y_start, y_end],
        color="green")

    # Shared dome center at top of bottom frustum
    y_eq = y_end
    # middle dome (lower half, cyan)
    draw_hemisphere(ax,
                        center_y=y_eq,
                        radius=top_radius_bottom_tank,
                        lower=lh2_on_top,
                        color="cyan",
                        label="Shared dome (LOX half)")

    # top frustum (orange)
    y2_start = bottom_radius_bottom_tank + height_bottom_tank
    y2_end = y2_start + height_top_tank
    ax.plot(
        [-bottom_radius_top_tank, -top_radius_top_tank],
        [y2_start, y2_end],
        color="orange", label="LOX frustum")
    ax.plot(
        [bottom_radius_top_tank, top_radius_top_tank],
        [y2_start, y2_end],
        color="orange")

    # LH2 top dome (red, upper half)
    draw_hemisphere(ax,
                    center_y=y2_end,
                    radius=top_radius_top_tank,
                    lower=False,
                    color="red",
                    label="LOX top dome")

    # Interface line (black)
    ax.axhline(y=y_eq,
               color="black",
               linestyle="--",
               linewidth=1,
               label="Frustum interface")

    ax.set_title(title)
    ax.set_xlabel("Radius [m]")
    ax.set_ylabel("Height [m]")
    ax.set_aspect("equal", "box")
    ax.set_xlim(-5, 5)
    ax.grid(True)
    ax.legend(loc="upper right")
    plt.tight_layout()
    plt.show()


if __name__ == "__main__":
    # Geometry
    radius_ratio = 0.5
    top_radius_ratio = 0.8
    bottom_radius_ratio = radius_ratio / top_radius_ratio

    # Tank dimensions
    bottom_radius = 5  # m
    middle_radius = bottom_radius * bottom_radius_ratio
    top_radius = middle_radius * top_radius_ratio
    caps_height_radius_ratio = 0.5  # ratio of the height of the caps to the radius

    lh2_volume = 520 # m^3
    o2_volume = 96 # m^3

    lh2_on_top = False  # Set to True if LH2 tank is on top, False if LOX is on top

    if lh2_on_top:
        top_volume = lh2_volume
        bottom_volume = o2_volume
    else:
        top_volume = o2_volume
        bottom_volume = lh2_volume

    # top tank
    top_height = tank_sizing.calculate_frustum_tank_length(
        top_volume,
        top_radius,
        middle_radius,
        caps_height_radius_ratio,
        subtract_bottom_cap= not lh2_on_top,
        subtract_top_cap= False
    )
    # bottom tank
    bottom_height = tank_sizing.calculate_frustum_tank_length(
        bottom_volume,
        middle_radius,
        bottom_radius,
        caps_height_radius_ratio,
        subtract_bottom_cap=False,
        subtract_top_cap=lh2_on_top
    )

    top_tank = (top_height, middle_radius,top_radius)
    bottom_tank = (bottom_height, bottom_radius, middle_radius)

    plot_stacked_tanks(top_tank, bottom_tank, lh2_on_top=lh2_on_top)
