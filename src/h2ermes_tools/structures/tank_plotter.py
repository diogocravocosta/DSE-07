import matplotlib.pyplot as plt
import numpy as np
from h2ermes_tools.structures import tank_sizing


def draw_hemisphere(ax, center_y, radius, lower=True, color="k", label=None):
    theta = np.linspace(0, np.pi, 100)
    x = radius * np.cos(theta)
    if lower:
        y = center_y - radius * np.sin(theta)
    else:
        y = center_y + radius * np.sin(theta)
    ax.plot(x, y, color, label=label)
    ax.plot(-x, y, color)


def plot_stacked_tanks(lox, lh2, title="Stacked LOX + LH2 Tank Profile"):
    """
    Plot two tanks (LOX at bottom, LH2 on top) as a single cross-section.
    Each tank: (h, R, r, total_length)
    The tanks are plotted with y as the axis of symmetry (vertical), x as radius (horizontal).
    """
    # Unpack tank geometry
    h_LOX, R_LOX, r_LOX  = lox
    h_LH2, R_LH2, r_LH2 = lh2

    fig, ax = plt.subplots(figsize=(6, 16))

    # LOX bottom dome (blue)
    draw_hemisphere(ax, center_y=R_LOX, radius=R_LOX, lower=True, color="blue", label="LOX bottom dome")

    # LOX frustum (green)
    y_start = R_LOX
    y_end = R_LOX + h_LOX
    ax.plot([-R_LOX, -r_LOX], [y_start, y_end], color="green", label="LOX frustum")
    ax.plot([R_LOX, r_LOX], [y_start, y_end], color="green")

    # Shared dome center at top of LOX frustum
    y_eq = y_end
    # LOX top dome (lower half, cyan)
    draw_hemisphere(ax, center_y=y_eq, radius=r_LOX, lower=True, color="cyan", label="Shared dome (LOX half)")
    # LH2 bottom dome (upper half, magenta)
    # draw_hemisphere(ax, center_y=y_eq, radius=r_LOX, lower=False, color="magenta", label="Shared dome (LH2 half)")

    # LH2 frustum (orange)
    y2_start = R_LOX + h_LOX
    y2_end = y2_start + h_LH2
    ax.plot([-r_LOX, -r_LH2], [y2_start, y2_end], color="orange", label="LH2 frustum")
    ax.plot([r_LOX, r_LH2], [y2_start, y2_end], color="orange")

    # LH2 top dome (red, upper half)
    draw_hemisphere(ax, center_y=y2_end, radius=r_LH2, lower=False, color="red", label="LH2 top dome")

    # Interface line (black)
    ax.axhline(y=y_eq, color="black", linestyle="--", linewidth=1, label="Frustum interface")

    ax.set_title(title)
    ax.set_xlabel("Radius [m]")
    ax.set_ylabel("Height [m]")
    ax.set_aspect("equal", "box")
    ax.set_xlim(-5, 5)
    ax.set_ylim(-5, 20)
    ax.grid(True)
    ax.legend(loc="upper right")
    plt.tight_layout()
    plt.show()


if __name__ == "__main__":
    # Get tank geometries
    lox = tank_sizing.calculate_tank_length_LOX(
        tank_model=None,
        volume=tank_sizing.LOX_volume,
        radius_ratio=0.9,
        phi=tank_sizing.phi,
        R=tank_sizing.R,
    )
    lh2 = tank_sizing.calculate_tank_length_LH2(
        tank_model=None,
        volume= tank_sizing.LH2_volume,
        radius_ratio=0.9,
        phi=tank_sizing.phi,
        LOX_radius=lox[1],
    )
    plot_stacked_tanks(lox, lh2)
