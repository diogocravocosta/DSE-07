"""
This script produces the risk maps for the project plan (organizational risks).
Can be reused for the risk map later on as well.
"""
# Import statements
import matplotlib.pyplot as plt
from collections import defaultdict


# Class definition
class RiskMap:
    # data is a dictionary of all points, which have probability, severity and an ID
    # mit is only used for file name when saving, changing from original to mitigated
    def __init__(self, data, mit=False):
        # initialize attributes
        self.data = data
        self.mit = mit

        # initialize plot
        self.fig, self.ax = plt.subplots(1, 1, figsize=(8, 8))

        # Set some of the properties for layout
        self.ax.grid(True)

        self.ax.set_xlim([-0.1, 5.1])
        self.ax.set_ylim([-0.1, 5.1])

        self.ax.set_xlabel('Probability')
        self.ax.set_ylabel('Severity')

        # Used for having the ID's not overlap
        coord_count = defaultdict(int)

        # Plot all the points in the data
        for point in self.data.keys():
            P = self.data[point]['P']
            S = self.data[point]['S']
            ID = self.data[point]['ID']

            # Different color and marker for internal/external risks
            color = 'r' if 'INT' in ID else 'b'
            marker = 'o' if 'INT' in ID else 's'

            self.ax.scatter(P, S, color=color, marker=marker)

            # Again used for not having IDs overlap
            count = coord_count[(P, S)]
            offset = 0.15 * count
            self.ax.text(P, S + offset, ID)

            coord_count[(P, S)] += 1

        # Create points outside the plotting region for the legend
        self.ax.scatter(-1, -1, color='r', marker='o', label='Internal Risks')
        self.ax.scatter(-1, -1, color='b', marker='s', label='External Risks')

        # show the legend
        self.ax.legend(loc='lower left')

        # set tight layout
        self.fig.tight_layout()

    def show(self):
        self.fig.show()

    def save(self):
        filename = '7_H2GO_RiskMap_Original.eps' if not self.mit else '7_H2GO_RiskMap_Mitigated.eps'
        self.fig.savefig(filename, format='eps')


def get_original_data():
    # Values from the project plan. Reversed for order of IDs on the risk map
    return dict(reversed([
        (1, {'P': 5, 'S': 5, 'ID': 'RI-INT-1'}),
        (2, {'P': 5, 'S': 5, 'ID': 'RI-INT-2'}),
        (3, {'P': 4, 'S': 5, 'ID': 'RI-INT-3'}),
        (4, {'P': 3, 'S': 5, 'ID': 'RI-INT-4'}),
        (5, {'P': 4, 'S': 4, 'ID': 'RI-INT-5'}),
        (6, {'P': 3, 'S': 4, 'ID': 'RI-INT-6'}),
        (7, {'P': 4, 'S': 4, 'ID': 'RI-INT-7'}),
        (8, {'P': 4, 'S': 4, 'ID': 'RI-INT-8'}),
        (9, {'P': 5, 'S': 5, 'ID': 'RI-EXT-1'}),
        (10, {'P': 5, 'S': 5, 'ID': 'RI-EXT-2'}),
        (11, {'P': 4, 'S': 5, 'ID': 'RI-EXT-3'}),
        (12, {'P': 4, 'S': 4, 'ID': 'RI-EXT-4'}),
        (13, {'P': 4, 'S': 3, 'ID': 'RI-EXT-5'}),
        (14, {'P': 3, 'S': 4, 'ID': 'RI-EXT-6'}),
        (15, {'P': 3, 'S': 3, 'ID': 'RI-EXT-7'}),
        (16, {'P': 3, 'S': 2, 'ID': 'RI-EXT-8'}),
        (17, {'P': 2, 'S': 5, 'ID': 'RI-EXT-9'}),
        (18, {'P': 3, 'S': 4, 'ID': 'RI-EXT-10'})
    ]))


def get_mitigated_data():
    # Values from the project plan. Reversed for order of IDs on the risk map
    return dict(reversed([
        (1, {'P': 3, 'S': 2, 'ID': 'RI-INT-1'}),
        (2, {'P': 3, 'S': 2, 'ID': 'RI-INT-2'}),
        (3, {'P': 3, 'S': 3, 'ID': 'RI-INT-3'}),
        (4, {'P': 2, 'S': 3, 'ID': 'RI-INT-4'}),
        (5, {'P': 2, 'S': 3, 'ID': 'RI-INT-5'}),
        (6, {'P': 2, 'S': 3, 'ID': 'RI-INT-6'}),
        (7, {'P': 2, 'S': 2, 'ID': 'RI-INT-7'}),
        (8, {'P': 3, 'S': 3, 'ID': 'RI-INT-8'}),
        (9, {'P': 3, 'S': 3, 'ID': 'RI-EXT-1'}),
        (10, {'P': 3, 'S': 3, 'ID': 'RI-EXT-2'}),
        (11, {'P': 4, 'S': 2, 'ID': 'RI-EXT-3'}),
        (12, {'P': 3, 'S': 3, 'ID': 'RI-EXT-4'}),
        (13, {'P': 3, 'S': 2, 'ID': 'RI-EXT-5'}),
        (14, {'P': 3, 'S': 2, 'ID': 'RI-EXT-6'}),
        (15, {'P': 2, 'S': 2, 'ID': 'RI-EXT-7'}),
        (16, {'P': 1, 'S': 1, 'ID': 'RI-EXT-8'}),
        (17, {'P': 2, 'S': 5, 'ID': 'RI-EXT-9'}),
        (18, {'P': 3, 'S': 3, 'ID': 'RI-EXT-10'})
    ]))


if __name__ == '__main__':
    # Create risk maps
    RM = RiskMap(get_original_data())
    RM2 = RiskMap(get_mitigated_data(), mit=True)

    # Show them
    RM.show()
    RM2.show()

    # Save them if wanted
    if True:
        RM.save()
        RM2.save()
