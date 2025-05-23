import matplotlib.pyplot as plt

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

# Google‐style palette
colors = {
    'GNC': '#4285F4',
    'Communication': '#EA4335',
    'Thermal control': '#FBBC05',
    'Data processing': '#34A853',
    'Pumps and compressors': '#AB47BC',
    'Actuators and control surfaces': '#FF7043'
}

def prep(data):
    # filter out zero entries & sort ascending by value
    d = {k: v for k, v in data.items() if v > 0}
    return dict(sorted(d.items(), key=lambda item: item[1]))

avg = prep(average_loads)
pk  = prep(peak_loads)

fig, (ax1, ax2) = plt.subplots(2, 1, figsize=(8, 10))

# tighten vertical spacing
fig.subplots_adjust(hspace=0.3)

# --- Average loads pie + legend ---
wedges1, _ = ax1.pie(
    avg.values(),
    labels=[''] * len(avg),
    startangle=90,
    colors=[colors[k] for k in avg],
    wedgeprops=dict(edgecolor='w')
)
ax1.axis('equal')

# build legend in descending order, shifted closer (x=0.95)
avg_items_desc = list(avg.items())[::-1]
wedges1_desc = list(wedges1)[::-1]
legend_labels1 = [f"{k}: {v} W" for k, v in avg_items_desc]
ax1.legend(
    wedges1_desc,
    legend_labels1,
    title="Average Loads",
    loc="center left",
    bbox_to_anchor=(0.95, 0.5)
)

# --- Peak loads pie + legend ---
wedges2, _ = ax2.pie(
    pk.values(),
    labels=[''] * len(pk),
    startangle=90,
    colors=[colors[k] for k in pk],
    wedgeprops=dict(edgecolor='w')
)
ax2.axis('equal')

# build legend in descending order, shifted closer (x=0.95)
pk_items_desc = list(pk.items())[::-1]
wedges2_desc = list(wedges2)[::-1]
legend_labels2 = [f"{k}: {v} W" for k, v in pk_items_desc]
ax2.legend(
    wedges2_desc,
    legend_labels2,
    title="Peak Loads",
    loc="center left",
    bbox_to_anchor=(0.95, 0.5)
)

plt.tight_layout()
plt.show()
# plt.savefig('pie_chart.png', dpi=300, bbox_inches='tight')