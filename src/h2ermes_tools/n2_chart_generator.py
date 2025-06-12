import matplotlib.pyplot as plt
import numpy as np
from matplotlib.patches import Rectangle

def plot_n2(elements_io: dict, title="N2 Chart", save=False):
    """
    elements_io: dict where
      key   = element name (str)
      value = dict with two keys:
        'inputs'  : list of (source_element_name, label) tuples
        'outputs' : list of (target_element_name, label) tuples

    Example:
      {
        "X": { "inputs":  [], 
               "outputs": [("Y","α"),("Z","β")] },
        "Y": { "inputs": [("X","α")], 
               "outputs":[("Z","δ")] },
        "Z": { "inputs":[("X","β"),("Y","δ")], 
               "outputs": [] }
      }
    """
    # 1. Build an ordered list of elements
    elements = list(elements_io.keys())
    N = len(elements)
    idx = {el: i for i, el in enumerate(elements)}

    # 2. Initialize an empty interaction matrix of strings
    inter = [["" for _ in range(N)] for _ in range(N)]

    # 3. Fill in outputs
    for src, io in elements_io.items():
        i = idx[src]
        for tgt, label in io.get("outputs", []):
            j = idx[tgt]
            inter[i][j] = label

    # 4. Plotting
    aspect_ratio = 1.4142  # Golden ratio for aesthetics
    size_multiplier = 1.6  # Base size multiplier
    figsize = (max(9, N * size_multiplier), max(6, N * size_multiplier / aspect_ratio))
    fig, ax = plt.subplots(figsize=figsize)

    # grid lines
    ax.set_xticks(np.arange(N+1)-0.5, minor=True)
    ax.set_yticks(np.arange(N+1)-0.5, minor=True)
    ax.grid(which="minor", color="black", linewidth=1)
    ax.tick_params(which="minor", size=0)

    # Remove axis ticks
    ax.set_xticks([])
    ax.set_yticks([])

    # labels
    # ax.set_xticks(np.arange(N))
    # ax.set_yticks(np.arange(N))
    # ax.set_xticklabels(elements)
    # ax.set_yticklabels(elements)
    # ax.xaxis.tick_top()
    # ax.xaxis.set_label_position('top')
    ax.invert_yaxis()
    fig.tight_layout()

    # fill cells
    for i in range(N):
        for j in range(N):
            if i == j:
                # diagonal: fill entire cell
                ax.add_patch(Rectangle((j-0.5, i-0.5), 1, 1, facecolor="lightgray", edgecolor="black"))
                ax.text(j, i, elements[i], ha="center", va="center", fontsize=10, weight="bold")
            elif inter[i][j]:
                # interaction: fill entire cell
                ax.add_patch(Rectangle((j-0.5, i-0.5), 1, 1, facecolor="skyblue", edgecolor="black"))
                ax.text(j, i, inter[i][j], ha="center", va="center", fontsize=9)

    plt.title(title, pad=20)
    # plt.tight_layout()

    if save:
        fig.savefig("chart.svg")
        print("Saved N2 chart as chart.svg")

    plt.show()


if __name__ == "__main__":
    # Example usage:
    elements_io = {
        "Provide power": {
            "inputs":  [],
            "outputs": [("Provide structural support","α"),("Provide rigidity","β")]
        },
        "Provide structural support": {
            "inputs": [("Provide power","α")],
            "outputs":[("Provide rigidity","δ")]
        },
        "Provide rigidity": {
            "inputs":[("Provide power","β"),("Provide structural support","δ")],
            "outputs":[]
        },
        "Provide cooling": {
            "inputs":  [],
            "outputs": [("Provide power","γ"),("Provide structural support","ε")]
        },
    }

    plot_n2(elements_io, title=r"H$_2$ERMES Active Heat Shield N2 Chart", save=True)
