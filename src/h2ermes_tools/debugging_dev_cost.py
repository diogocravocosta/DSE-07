import numpy as np

def calc_development_cost(M, flights_per_vehicle):
    L_d = 0.9
    s_BAU = 0.6
    s_COM = 0.2
    q = 0.08
    c_p = (s_COM*q + 1) / (s_BAU*q + 1)
    print(f'c_p: {c_p:.2f}')
    M_PA_perc = 5.25
    learning_factor = 0.9
    delta_TRL = 4

    all_equipment_names = ["Pressurizant Tank", "Fuel Tank", "Oxidizer Tank", "Thrust Cone", "Skirt", "Thermal Control", 
                           "Engine(s)", "Thrust Vector Control", "Pressurizant System", "Pipes", "Valves", "Stage Harness", 
                           "Payload Adapter", "Payload Fairing", "Comms", "Power", "Data Handling", "GNC", 
                           "Avionics Harness", "Attitude Control Module", "Interstage Structure"]

    all_HW = np.array([0, flights_per_vehicle, flights_per_vehicle, 0, 0, flights_per_vehicle, flights_per_vehicle,
                       15, flights_per_vehicle, 10, 15, 0, 0, 0, flights_per_vehicle, flights_per_vehicle, flights_per_vehicle, 
                       flights_per_vehicle, flights_per_vehicle, flights_per_vehicle, flights_per_vehicle])

    a_all = [19.99465, 19.99465, 19.99465, 2.799300, 2.799300, 2.799300, 31.48271, 33.90978, 11.50618, 8.958770, 8.958770,
             27.45211, 26.01794, 23.59239, 51.11253, 42.01174, 141.6820, 69.05491, 27.45211, 257.8420, 6.70369]
    b_all = [0.71253, 0.71253, 0.71253, 0.91199, 0.91199, 0.91199, 0.78811, 0.60977, 1.06948, 0.68815, 0.68815, 0.44623,
             0.44623, 0.70000, 0.80000, 0.80000, 0.80000, 0.82458, 0.44623, 0.75000, 0.68041]

    # === Filter out zero-mass components ===
    M = np.array(M)
    valid_mask = M != 0

    valid_mask = M != 0

    M = M[valid_mask]
    HW = all_HW[valid_mask]
    equipment_names = [name for i, name in enumerate(all_equipment_names) if valid_mask[i]]
    a = [a_all[i] for i in range(len(a_all)) if valid_mask[i]]
    b = [b_all[i] for i in range(len(b_all)) if valid_mask[i]]


    # === Compute number of test units for development ===
    n = np.ones_like(HW, dtype=float)
    with np.errstate(divide='ignore', invalid='ignore'):
        mask = HW > 0
        n[mask] = np.ceil(flights_per_vehicle / HW[mask])

    DEV = []

    for i in range(len(M)):
        component_cost = a[i] * M[i]**b[i]
        T1 = []
        T1_component = component_cost * (n[i] ** (np.log(learning_factor)/np.log(2)))
        T1.append(T1_component)
        print(f'Component: {equipment_names[i]}, T1_component: {T1_component:.2f} k€')
        DD = (3 + delta_TRL) * T1_component
        DM = 0.3 * T1_component
        EM = 1.3 * T1_component
        PFM = 1.5 * T1_component

        FM1 = (1 - M_PA_perc / 100) 

        ENG = DD * FM1
        STH = 3.1
        MAIT = FM1 * STH * L_d * HW[i]

        DEV_component = c_p * ((ENG + (MAIT + ENG) * M_PA_perc / 100) + (FM1 * STH * L_d * HW[i]))
        DEV.append(DEV_component)
    T1_total = np.sum(T1)
    print(f'T1_total: {T1_total/1000:.2f} M€')
    DEV_total = np.sum(DEV)
    return equipment_names, DEV, DEV_total

# === MAIN TEST ===
if __name__ == "__main__":
    M = [0, 24477, 3817, 0, 0, 0, 336, 0, 437.19, 92.7828, 92.7828, 0, 0, 850, 0.5, 240, 0, 60, 100, 923, 0]
    equipment_names, DEV_list, DEV_total = calc_development_cost(M, flights_per_vehicle=25)

    for name, cost in zip(equipment_names, DEV_list):
        print(f'Component: {name:30} | Development cost: {cost:,.2f} k€')
    
    print(f'\nTotal development cost: {DEV_total/1000:.2f} M€')
