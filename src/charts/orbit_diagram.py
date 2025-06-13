import matplotlib.pyplot as plt
import numpy as np

import data.constants as cn


def plot_mission_till_rendezvous(initial_orbit:float,
                                 target_orbit:float,
                                 up_phasing_altitude:float,
                                 down_phasing_altitude:float,
                                 reentry_altitude:float,
                                 subtract_earth_radius:bool,
                                 plot_inner_body:bool
                                 ):
    fig, ax = plt.subplots(subplot_kw={'projection': 'polar'})

    if subtract_earth_radius:
        inner_body_radius = 0
    else:
        inner_body_radius = cn.earth_radius

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

    # plot transfer orbit
    plot_ellipse(ax=ax,
                 apoapsis=target_orbit + cn.earth_radius,
                 periapsis=initial_orbit + cn.earth_radius,
                 start_theta=np.pi/2,
                 end_theta=np.pi *3/2,
                 periapsis_theta=np.pi/2,
                 label='Transfer orbit',
                 subtract_earth_radius=subtract_earth_radius
                 )


    # plot down phasing orbit
    plot_ellipse(ax=ax,
                 apoapsis=target_orbit + cn.earth_radius,
                 periapsis=down_phasing_altitude + cn.earth_radius,
                 start_theta=0,
                 end_theta=2*np.pi,
                 periapsis_theta=0,
                 label='Down phasing orbit',
                 subtract_earth_radius=subtract_earth_radius
                 )

    # plot reentry orbit
    plot_ellipse(ax=ax,
                 apoapsis=target_orbit + cn.earth_radius,
                 periapsis=reentry_altitude + cn.earth_radius,
                 start_theta=np.pi,
                 end_theta=2*np.pi,
                 periapsis_theta=0,
                 label='Reentry orbit',
                 subtract_earth_radius=subtract_earth_radius
                 )

    ax.set_rmax((up_phasing_altitude + inner_body_radius)/1000 * 1.1)
    ax.grid(True)
    ax.legend()
    plt.show()

def plot_circle(ax:plt.Axes,
                radius:float,
                start_theta:float=0,
                end_theta:float=2*np.pi,
                color:str|None = None,
                label:str|None = None,
                ):
    step = np.pi/100

    thetas = np.arange(start_theta, end_theta + step, step)
    radia = np.ones_like(thetas) * radius

    ax.plot(thetas, radia/1000, color = color, label=label)

def plot_ellipse(ax:plt.Axes,
                 apoapsis:float,
                 periapsis:float,
                 start_theta:float,
                 end_theta:float,
                 periapsis_theta:float,
                 color:str|None = None,
                 label:str|None = None,
                 subtract_earth_radius:bool = False
                 ):
    semi_major_axis = (apoapsis + periapsis)/2
    eccentricity = (apoapsis - periapsis)/(apoapsis+periapsis)

    step = np.pi/100

    thetas = np.arange(start_theta, end_theta + step, step)
    radia = semi_major_axis * (1-eccentricity**2)/(1+eccentricity*np.cos(thetas - periapsis_theta))

    if subtract_earth_radius:
        radia -= cn.earth_radius

    ax.plot(thetas, radia/1000, color = color, label=label)

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
                                 subtract_earth_radius=False,
                                 plot_inner_body=True
                                 )

    plot_mission_till_rendezvous(initial_orbit=2e5,
                                 target_orbit=6e5,
                                 up_phasing_altitude=7e5,
                                 down_phasing_altitude=2e5,
                                 reentry_altitude=5e4,
                                 subtract_earth_radius=False,
                                 plot_inner_body=False
                                 )

    plot_mission_till_rendezvous(initial_orbit=2e5,
                                 target_orbit=6e5,
                                 up_phasing_altitude=7e5,
                                 down_phasing_altitude=2e5,
                                 reentry_altitude=5e4,
                                 subtract_earth_radius=True,
                                 plot_inner_body=False
                                 )