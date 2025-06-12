

def calc_development_cost(M):
    learning_factor = 0.9 # Learning factor for development cost
    equipment_names = ["Pressurizant Tank",
    "Fuel Tank",
    "Oxidizer Tank",
    "Thrust Cone",
    "Skirt",
    "Thermal Control",
    "Engine(s)",
    "Thrust Vector Control",
    "Pressurizant System",
    "Pipes",
    "Valves",
    "Stage Harness",
    "Payload Adapter",
    "Payload Fairing",
    "Comms",
    "Power",
    "Data Handling",
    "GNC",
    "Avionics Harness",
    "Attitude Control Module",
    "Interstage Structure"]
    a = [19.99465, 19.99465, 19.99465, 2.799300, 2.799300, 2.799300, 31.48271, 33.90978, 11.50618, 8.958770, 8.958770, 27.45211, 26.01794, 23.59239, 51.11253, 42.01174, 141.6820, 69.05491, 27.45211, 257.8420, 6.70369]
    b = [0.71253, 0.71253, 0.71253, 0.91199, 0.91199, 0.91199, 0.78811, 0.60977, 1.06948, 0.68815, 0.68815, 0.44623, 0.44623, 0.70000, 0.80000, 0.80000, 0.80000, 0.82458, 0.44623, 0.75000, 0.68041]
    for i in range len(equipment_names): 
        
        component_cost = a[i] * M[i]**b[i]
        first_unit_cost_component = component_cost*(n**((np.log(learning_factor))/np.log(2)))
        T = np.sum(component_cost)
    DD = 3 * T
    DM = 0.3 * T
    EM = 1.3 * T
    PFM = 1.5 * T

    DD_prime = DD + delta_TRL
    FM1 = 1 - M_PA_perc / 100

    ENG = DD_prime * FM1
    STH = DM + EM + PFM
    MAIT = FM1 * STH * L_d * HW_i

    DEV = c_p * ((ENG + (MAIT + ENG) * M_PA_perc) + (FM1 * STH * L_d * HW_i))


M = [0, 24477, 3817]