import numpy as np


class TankSizer:
    def __init__(self, volume: float, taper: float = 0.5):
        self.volume = volume
        self.taper = taper

        self.print = True

        self.find_configs()

        if self.print:
            self.print_config()

    def find_configs(self):
        phi_range = np.arange(1, 90, 1)

        r_bottom = (self.volume * 3 * np.tan(np.radians(phi_range)) / (
                    np.pi * (1 - self.taper) * (1 + self.taper + self.taper ** 2))) ** (1 / 3)

        r_top = r_bottom * self.taper

        h = (r_bottom - r_top) / np.tan(np.radians(phi_range))

        for i, val in enumerate(r_bottom):
            if 2 * val >= 10:
                continue
            if 2 * val < 7:
                continue
            highest_idx = i

        try:
            self.phi = phi_range[highest_idx]
            self.r_top = r_top[highest_idx]
            self.r_bottom = r_bottom[highest_idx]
            self.h = h[highest_idx]
        except UnboundLocalError:
            print('No valid configuration found. Try changing the taper parameter')
            self.print = False

    def print_config(self):
        print(f'phi: {self.phi} deg\n'
              f'taper: {self.taper} -\n'
              f'r_top: {self.r_top} m\n'
              f'r_bottom: {self.r_bottom} m\n'
              f'h: {self.h} m\n')


if __name__ == '__main__':
    ts = TankSizer(619, taper=0.1)
