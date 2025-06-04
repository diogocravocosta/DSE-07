import matplotlib.pyplot as plt
import pandas as pd

from h2ermes_tools.delta_v.landing_burn import landing_burn_no_drag, landing_burn_with_drag

ballistic_coefficients = [float('inf'), 1000, 500, 200, 100]
thrust_to_weight_ratios = [1.1, 1.2, 1.5, 2.0, 2.5, 3.0, 3.5]
initial_velocities = [50, 100, 150, 200, 250, 300]

# for each ballistic coefficient make a heat map of landing burn delta v with thrust to weight ratio and initial velocity on the axes
def plot_landing_burn_heatmaps():

    for bc in ballistic_coefficients:
        delta_vs = []
        for ttw in thrust_to_weight_ratios:
            row = []
            for iv in initial_velocities:
                _, delta_v = landing_burn_no_drag(iv, ttw)
                row.append(delta_v)
            delta_vs.append(row)

        df = pd.DataFrame(delta_vs, index=thrust_to_weight_ratios, columns=initial_velocities)
        plt.figure(figsize=(10, 6))
        plt.title(f'Landing Burn Delta V Heatmap (BC={bc})')
        plt.xlabel('Initial Velocity (m/s)')
        plt.ylabel('Thrust-to-Weight Ratio')
        plt.imshow(df, aspect='auto', cmap='viridis', origin='lower',
                   extent=[initial_velocities[0], initial_velocities[-1],
                           thrust_to_weight_ratios[0], thrust_to_weight_ratios[-1]])
        plt.colorbar(label='Delta V (m/s)')
        plt.xticks(initial_velocities)
        plt.yticks(thrust_to_weight_ratios)
        plt.grid(False)
        plt.show()

if __name__ == "__main__":
    plot_landing_burn_heatmaps()