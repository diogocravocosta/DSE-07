import matplotlib.pyplot as plt
import numpy as np

import data.constants as cn

def plot_mission_till_rendezvous(initial_orbit:float,
                                 target_orbit:float,
                                 up_phasing_altitude:float,
                                 down_phasing_altitude:float,
                                 reentry_altitude:float,
                                 inner_body_radius:float,
                                 plot_inner_body:bool
                                 ):
    fig, ax = plt.subplots(subplot_kw={'projection': 'polar'})

    # Set radial axis label
    ax.set_xlabel("Radius [km]")

    # plot earth
    if plot_inner_body is True:
        plot_circle(ax=ax,
                    radius=inner_body_radius,
                    color='black',
                    label="Earth radius"
                    )

    # plot initial orbit
    plot_circle(ax,
                initial_orbit + inner_body_radius,
                label="Initial orbit"
                )

    # plot target orbit
    plot_circle(ax,
                target_orbit + inner_body_radius,
                label='Target orbit'
                )

    ax.set_rmax((up_phasing_altitude + inner_body_radius)/1000 * 1.1)
    ax.grid(True)
    ax.legend()
    plt.show()

def plot_circle(ax:plt.Axes,
                radius:float,
                color:str|None = None,
                label:str|None = None
                ):
    steps = 100
    thetas = np.linspace(0, 2*np.pi, steps)
    radia = np.ones(steps) * radius

    ax.plot(thetas, radia/1000, color = color, label=label)

def plot_ellipse(ax:plt.Axes,
                 apoapsis:float,
                 periapsis:float,
                 start_theta:float,
                 end_theta:float
                 color:str|None = None,
                 label:str|None = None
                 ):

def plot_ascent(ax:plt.Axes,
                start_theta:float = 0,
                end_theta:float = 1.819420e-01,
                altitude:float = 2e5) -> None:
    """

    Args:
        ax:
        start_theta:
        end_theta:
        altitude:
    """
if __name__ == '__main__':
    plot_mission_till_rendezvous(initial_orbit=2e5,
                                 target_orbit=6e5,
                                 up_phasing_altitude=7e5,
                                 down_phasing_altitude=2e5,
                                 reentry_altitude=5e4,
                                 inner_body_radius=cn.earth_radius,
                                 plot_inner_body=True
                                 )

    plot_mission_till_rendezvous(initial_orbit=2e5,
                                 target_orbit=6e5,
                                 up_phasing_altitude=7e5,
                                 down_phasing_altitude=2e5,
                                 reentry_altitude=5e4,
                                 inner_body_radius=cn.earth_radius,
                                 plot_inner_body=False
                                 )

    plot_mission_till_rendezvous(initial_orbit=2e5,
                                 target_orbit=6e5,
                                 up_phasing_altitude=7e5,
                                 down_phasing_altitude=2e5,
                                 reentry_altitude=5e4,
                                 inner_body_radius=0,
                                 plot_inner_body=False
                                 )