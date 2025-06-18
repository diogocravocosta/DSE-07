import numpy as np

import data.constants as cn

def calculate_launch_windows(relative_phase_change:float, # deg
                             time_between_passes:float, # s
                             number_of_years:int, # -
                             max_catch_up_angle: float # deg
                             ):
    n = number_of_years * cn.tropical_year // time_between_passes

    angles = np.arange(0, n*relative_phase_change, relative_phase_change)
    remainders = np.mod(angles, 360)
    launch_windows = np.sum(remainders<max_catch_up_angle)
    return launch_windows/number_of_years



if __name__ == "__main__":
    for n_year in range(1, 101):
        launch_windows_per_year = calculate_launch_windows(209.066519260498,
                                 23.4650730466474 * 3600,
                                 n_year,
                                 144.7835
                                 )

    print(f"Launch windows per year:{launch_windows_per_year:.2f} with {n_year} years")

    # import numpy as np
    # import matplotlib.pyplot as plt
    #
    # # Example data generation
    # synodic_period = 23.46 / 24  # in days
    # phase_increment = 209.07
    # num_periods = 500
    # times = np.arange(num_periods) * synodic_period
    # phases = (np.arange(num_periods) * phase_increment) % 360
    #
    # # Plot
    # plt.figure(figsize=(10, 5))
    # plt.plot(times, phases, marker='o', linewidth=0, label='Relative Phase Angle')
    # plt.axhline(144, color='r', linestyle='--', label='Launch Window Threshold')
    # plt.fill_between(times, 0, 144, color='red', alpha=0.2)
    # plt.xlabel('Time (days)')
    # plt.ylabel('Relative Phase Angle (°)')
    # plt.title('Relative Phase Angle vs. Time')
    # plt.legend()
    # plt.grid(True)
    # plt.show()
