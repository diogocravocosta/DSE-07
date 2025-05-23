from data import read_data
import warnings
import copy
import numpy as np
import matplotlib.pyplot as plt

warnings.simplefilter(action='ignore', category=FutureWarning)
warnings.simplefilter(action='ignore', category=UserWarning)

plt.rcParams.update({'font.size': 14})

class TradeOff:
    def __init__(self, filepath: str, n_runs: int = 10, sens_mode: str = 'uniform', std_frac: float = 0.33) -> None:
        """
        Class to perform the trade-off.
        Args:
            filepath (str): Path to the Excel file used as input
            n_runs (int): Number of extra runs to do for the sensitivity analysis (default 10)
        """
        # Read the data from the Excel sheet
        self.concepts_df, self.criteria_df, self.matrix_df = read_data(filepath)
        self.n_runs = n_runs
        self.sens_mode = sens_mode
        self.std_frac = std_frac

        # Create internal dictionaries
        self.__create_concept_dictionary()
        self.__create_criteria_dictionary()
        self.__create_matrix_dictionary()

    def run(self):
        # Get the original winner (based on the chosen weights for the criteria)
        self.compute_scores()
        self.initial_winner = self.return_winner()

        # Run the sensitivity analysis
        self.run_sensitivity()

    def __create_concept_dictionary(self):
        # Create the dictionary
        self.concepts = {}
        for concept in self.concepts_df.iterrows():
            self.concepts[concept[1][0]] = concept[1][1]

    def __create_criteria_dictionary(self):
        # Create the dictionary
        self.criteria = {}
        for criteria in self.criteria_df.iterrows():
            self.criteria[criteria[1][0]] = criteria[1][1]

        self.normalize_criteria_weights()

    def normalize_criteria_weights(self):
        # Normalize the weights
        total = sum(self.criteria.values())
        # print(self.criteria.values())
        if total == 0:
            # in case of a total of zero, all weights will be set to 1 and then normalized (equal weights)
            for criteria, weight in self.criteria.items():
                self.criteria[criteria] = 1
                self.normalize_criteria_weights()
            return

        for key, val in self.criteria.items():
            self.criteria[key] = val / total

    def __create_matrix_dictionary(self):
        self.matrix = {}
        headers = list(self.matrix_df.columns)[1:]

        for concept in self.matrix_df.iterrows():
            self.matrix[concept[1][0]] = {}
            for idx, criteria in enumerate(headers):
                # value: concept[1][idx+1]
                self.matrix[concept[1][0]][criteria] = concept[1][idx + 1]

    def compute_scores(self):
        self.score_matrix = copy.deepcopy(self.matrix)

        # Get the individual criteria scores
        for concept, scores in self.score_matrix.items():
            for criteria, score in scores.items():
                self.score_matrix[concept][criteria] = self.score_matrix[concept][criteria] * self.criteria[criteria]

        # Compute the total score
        for concept, scores in self.score_matrix.items():
            self.score_matrix[concept]['Total'] = sum(scores.values())

    def return_winner(self):
        winner = max(self.score_matrix, key=lambda k: self.score_matrix[k]['Total'])

        return winner

    def run_sensitivity(self):
        self.__initialize_sensitivity_analysis()
        self.__perform_sensitivity_analysis()

    def __initialize_sensitivity_analysis(self):
        self.winner_matrix = {}
        for concept in self.concepts.keys():
            self.winner_matrix[concept] = 0

        self.winner_matrix[self.initial_winner] = 1

        self.original_criteria = copy.deepcopy(self.criteria)
        self.original_matrix = copy.deepcopy(self.matrix)

    def __perform_sensitivity_analysis(self):
        # Initialize a dictionary to save all used weights for a boxplot
        self.sensitivity_weights_lists = {}
        for criteria in self.criteria.keys():
            self.sensitivity_weights_lists[criteria] = []

        for i in range(self.n_runs):
            # Randomly reassign weights
            for criteria, weight in self.criteria.items():
                if self.sens_mode == 'uniform':
                    self.criteria[criteria] = np.random.uniform(0, 1)
                elif self.sens_mode == 'normal':
                    # raise NotImplementedError('Normal weight selection not yet implemented')
                    self.criteria[criteria] = max(np.random.normal(loc=self.original_criteria[criteria], scale=self.original_criteria[criteria]*self.std_frac, size=1)[0], 0)
                else:
                    raise ValueError('sens_mode must be either uniform or normal. Exiting...')

            # Normalize criteria weights
            self.normalize_criteria_weights()

            # Save the weights for the boxplots
            for criteria in self.criteria.keys():
                self.sensitivity_weights_lists[criteria].append(self.criteria[criteria])

            # Compute new scores
            self.compute_scores()

            # Get new winner and add one point to it
            winner = self.return_winner()
            self.winner_matrix[winner] += 1

            # Print update every 10_000 runs
            if i % 50_000 == 1:
                print(f'Running sensitivity analysis - run {i:_}/{self.n_runs:_}')

        # Compute percentage of wins per concept
        total_wins = sum(self.winner_matrix.values())
        self.winner_matrix_perc = copy.deepcopy(self.winner_matrix)
        for key, val in self.winner_matrix_perc.items():
            self.winner_matrix_perc[key] = val / total_wins * 100

        # Get relative best option (most total wins across the board)
        self.real_winner = max(self.winner_matrix_perc, key=lambda k: self.winner_matrix_perc[k])

    def print_results(self):
        print(f'Original Winner (using original criteria weights): {self.concepts[self.initial_winner]}\n'
              f'Sensitivity analysis result: {self.concepts[self.real_winner]} at {self.winner_matrix_perc[self.real_winner]}% winrate')

    def print_all_results(self):
        print(f'Range of criteria weights:')
        for crit in self.criteria.keys():
            print(f'{crit} weight range: [{min(self.sensitivity_weights_lists[crit])}-{max(self.sensitivity_weights_lists[crit])}]')

        print(f'\nNumber of wins per concept:')
        for concept in self.concepts.keys():
            print(f'{self.concepts[concept]}: {self.winner_matrix[concept]} wins ({self.winner_matrix_perc[concept]}%)')



    def create_weight_boxplots(self):
        self.weight_fig, self.weight_ax = plt.subplots(figsize=(6, 6))

        labels = list(self.sensitivity_weights_lists.keys())
        data = [self.sensitivity_weights_lists[key] for key in labels]

        self.weight_ax.boxplot(data, labels=labels)

        self.weight_ax.set_ylabel('Weight')
        self.weight_ax.set_ylim([-0.1, 1.1])

        self.weight_fig.tight_layout()

        self.weight_fig.show()

    def create_winner_barchart(self, type='perc'):
        self.win_fig, self.win_ax = plt.subplots(figsize=(6, 6))

        if type == 'abs':
            data = self.winner_matrix.values()
            self.win_ax.set_ylabel('Number of wins (-)')
            self.win_ax.set_ylim([0, self.n_runs])
        elif type == 'perc':
            data = self.winner_matrix_perc.values()
            self.win_ax.set_ylabel('Winrate (%)')
            self.win_ax.set_ylim([0, 100])

        self.win_ax.bar(self.concepts.values(), data)


        self.win_ax.set_xticklabels(self.concepts.values(), rotation=80, ha='right')
        self.win_fig.tight_layout()
        self.win_fig.show()

    def save_figures(self):
        self.win_fig.savefig('H2GO_SensAnalysis_Winrate.eps', format='eps')
        self.weight_fig.savefig('H2GO_SensAnalysis_Weights.eps', format='eps')


if __name__ == '__main__':
    to = TradeOff('tradeOffInput.xlsx', n_runs=100)
    to.run()
    to.create_winner_barchart()
