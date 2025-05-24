import matplotlib.pyplot as plt
import random

# Color palette
PALETTE = [
    '#0C2340', '#00B8C8', '#0076C2', '#6F1D77', '#EF60A3', '#A50034',
    '#E03C31', '#EC6842', '#FFB81C', '#6CC24A', '#009B77', '#5C5C5C'
]

def plot_pie_dicts(dicts_and_titles, palette=PALETTE, legend_x=0.95, hspace=0.3):
    """
    dicts_and_titles : list of (data_dict, legend_title) tuples
    palette          : list of hex colors to choose from
    legend_x         : x‐position for legends (0–1)
    hspace           : vertical spacing between subplots
    """
    # 1) Collect all unique keys across all dicts
    unique_keys = []
    for data, _ in dicts_and_titles:
        for k in data:
            if k not in unique_keys:
                unique_keys.append(k)

    # 2) Assign each key a random color (no repeats until palette exhausted)
    available = palette.copy()
    random.shuffle(available)
    key_colors = {}
    for k in unique_keys:
        if available:
            key_colors[k] = available.pop()
        else:
            key_colors[k] = random.choice(palette)

    # 3) Set up subplots
    n = len(dicts_and_titles)
    fig, axes = plt.subplots(n, 1, figsize=(8, 4*n))
    if n == 1:
        axes = [axes]

    fig.subplots_adjust(hspace=hspace)

    # 4) For each dict, draw the pie & legend
    for ax, (data, title) in zip(axes, dicts_and_titles):
        # filter zeros & sort ascending
        filtered = {k: v for k, v in data.items() if v > 0}
        items = sorted(filtered.items(), key=lambda iv: iv[1])

        labels = [''] * len(items)
        sizes  = [v for _, v in items]
        colors = [key_colors[k] for k, _ in items]

        wedges, _ = ax.pie(
            sizes,
            labels=labels,
            startangle=90,
            colors=colors,
            wedgeprops=dict(edgecolor='w')
        )
        ax.axis('equal')

        # build legend in descending order
        items_desc   = items[::-1]
        wedges_desc  = wedges[::-1]
        legend_labels = [f"{k}: {v} kg" for k, v in items_desc]

        ax.legend(
            wedges_desc,
            legend_labels,
            title=title,
            loc='center left',
            bbox_to_anchor=(legend_x, 0.5)
        )

    plt.tight_layout()
    plt.show()


delta_v = {
    'Initial orbit injection': 6029,
    'Injection to 600 km transfer orbit': 114,
    'Circularization at 600 km': 112,
    'Deorbit burn': 157,
    'Atmospheric drag': 0,
    'Landing burn': 506,
    'Margin' : round((6029+114+112+157+506)*0.05)
}

average_loads = {
    'GNC': 908,
    'Communication': 205,
    'Thermal control': 360,
    'Data processing': 100,
    'Pumps and compressors': 2000,
    'Actuators and control surfaces': 0
}

peak_loads = {
    "GNC": 1143,
    "Communication": 300,
    "Thermal control": 460,
    "Data processing": 629,
    "Pumps and compressors": 6000,
    "Actuators and control surfaces": 475
}

empty_mass = {
    'Fuel tank': 35219,
    'Thermal control': 600,
    'Thrust vector control': 864,
    'Power': 719,
    'GNC': 100,
    'Attitude control': 4,
    'Oxidizer tank': 12051,
    'Engines': 5457,
    'Payload fairing': 4997,
    'Data handling': 100,
    'Avionics harness': 403,
    'Interstage structure': 462
}

total_mass = {
    'Empty mass': 28041,
    'Fuel': 218643
}

plot_pie_dicts([
    (empty_mass, 'Empty Mass'),
    (total_mass, 'Total Mass'),
])
