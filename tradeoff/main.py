from tradeOff import TradeOff
import time


def timer(func):
    def wrapper(*args, **kwargs):
        start = time.perf_counter()
        result = func(*args, **kwargs)
        end = time.perf_counter()
        print(f'Ran trade-off in {end - start:.4f} seconds')
        return result
    return wrapper


@timer
def main():
    file_path = 'tradeOffInput.xlsx'

    to = TradeOff(file_path, n_runs=10_000, sens_mode='uniform', std_frac=1/2)
    to.run()
    to.print_results()
    to.create_weight_boxplots()
    to.create_winner_barchart(type='perc')


if __name__ == '__main__':
    main()
