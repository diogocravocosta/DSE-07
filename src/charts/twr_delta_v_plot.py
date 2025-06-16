import matplotlib.pyplot as plt

delta_v_twr_relation = {'delta_v'               : [6210, 5937, 5810, 5739, 5673, 5645, 5629, 5618, 5610],
                        'thrust_to_weight_ratio': [0.6, 0.7, 0.8, 0.9, 1.0, 1.1, 1.2, 1.3, 1.4],
                        'offset'                : [26.5, 17, 11, 6.5, 3.7, 1.5, 0, -1.5, -3.5]
                        }

plt.plot(delta_v_twr_relation['thrust_to_weight_ratio'], delta_v_twr_relation['delta_v'], marker='x')


plt.xlabel('Thrust-to-weight Ratio [-]')
plt.ylabel('Delta V [m/s]')
# plt.title('Delta V vs Thrust to Weight Ratio')
plt.grid(True)
plt.show()