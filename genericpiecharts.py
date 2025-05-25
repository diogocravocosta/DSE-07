import matplotlib.pyplot as plt
import random

# Your color palette
PALETTE = [
    '#0C2340', '#00B8C8', '#0076C2', '#6F1D77', '#EF60A3', '#A50034',
    '#E03C31', '#EC6842', '#FFB81C', '#6CC24A', '#009B77', '#5C5C5C'
]

def plot_pie_dicts(
    dicts_info,
    palette=PALETTE,
    legend_x=0.95,
    hspace=0.3
):
    """
    dicts_info : list of (data_dict, title, unit, show_title, sort_flag) tuples
                  - data_dict  : mapping labels→values
                  - title      : legend title
                  - unit       : unit string to append to values
                  - show_title : bool, whether to show legend title
                  - sort_flag  : bool, whether to sort slices by value
    palette    : list of hex colors to choose from
    legend_x   : x‐position for legends (0–1)
    hspace     : vertical spacing between subplots
    """
    # 1) Collect all unique keys across all dicts
    unique_keys = []
    for data, *_ in dicts_info:
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
    n = len(dicts_info)
    fig, axes = plt.subplots(n, 1, figsize=(8*0.85, 4*0.65 * n), dpi=250)
    if n == 1:
        axes = [axes]
    fig.subplots_adjust(hspace=hspace)

    # 4) Draw each pie & legend
    for ax, (data, title, unit, show_title, sort_flag) in zip(axes, dicts_info):
        # filter out zero entries
        filtered = {k: v for k, v in data.items() if v > 0}
        # sort or preserve original order
        if sort_flag:
            items = sorted(filtered.items(), key=lambda iv: iv[1])
        else:
            items = list(filtered.items())

        sizes  = [v for _, v in items]
        colors = [key_colors[k] for k, _ in items]
        wedges, _ = ax.pie(
            sizes,
            labels=[''] * len(items),
            startangle=90,
            colors=colors,
            wedgeprops=dict(edgecolor='w')
        )
        ax.axis('equal')

        # build legend in descending order
        items_desc   = items[::-1]
        wedges_desc  = wedges[::-1]
        legend_labels = [f"{k}: {v} {unit}" for k, v in items_desc]

        legend_kwargs = {
            'handles': wedges_desc,
            'labels': legend_labels,
            'loc': 'center left',
            'bbox_to_anchor': (legend_x, 0.5)
        }
        if show_title:
            legend_kwargs.update({
                'title': title,
            })

        ax.legend(**legend_kwargs)

    plt.tight_layout()
    plt.show()


delta_v = {
    '5% margin': round((6029 + 114 + 112 + 157 + 506) * 0.05),
    'Landing burn': 506,
    'Deorbit burn': 157,
    'Circularization at 600 km': 112,
    'Injection to 600 km transfer orbit': 114,
    'Initial orbit injection': 6029,
    'Atmospheric drag': 0
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
    #'GNC': 100,
    #'Attitude control': 4,
    'Oxidizer tank': 12051,
    'Engines': 5457,
    'Payload fairing': 4997,
    #'Data handling': 100,
    #'Avionics harness': 403,
    'Interstage structure': 462,
    'Avionics': 100+4+100+403,
}

total_mass = {
    'Empty mass': 60906,
    'Propellant and payload': 218643
}

cost = {
    'Ground operations': 1.388,
    'Propellant': 0.226,
    'Flight and mission': 39.218,
    'Transportation': 1.323,
    'Fees and insurance': 1.389,
    'Indirect operational': 2.188,
}

# First: show title, second: sort
plot_pie_dicts([
    # (average_loads, "Average Loads", "W", False, True),
    # (peak_loads,    "Peak Loads",    "W", False, True)
    # (delta_v,       "Delta V", "m/s", False, False),
    #(empty_mass,    "Empty Mass", "kg", False, True),
    #(total_mass,    "Total Mass", "kg", False, True)
    (cost,          "Cost", "M€", False, True)
])