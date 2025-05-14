import numpy as np

#Nomenclature
#DEV = Development cost (engineering cost)
#MAIT = Manufacturing, Assembly, Integration and Test
#ENG = Engineering 
#FM1 = First flight model
#M/PA = Management and Production Assurance cost
#c_p = profit retention cost reduction factor
#STH = total system test hardware
#L_d = development learning factor
#HW = number of times a subsytem is re-used on the launch vehicle

#FIXED VALUES
L_d = 0.9 #as a default in Drenthe’s SOLSTICE cost estimating model
s_BAU = 0.6 #scope of subcontracted work under Business As Usual
s_COM = 0.2 #scope of subcontracted work under commercial development
q = 0.08 #average subcontractor profit
c_p = (s_COM*q+1)/(s_BAU*q+1)
M_PA_perc = 5.25 
concept_number = input("Enter concept number (e.g., 1, 2, 3): ")

if concept_number == "1" or concept_number == "2":
    delta_TRL = 4
elif concept_number == "3" or concept_number == "4":
    delta_TRL = 2
elif concept_number == "5" or concept_number == "6":
    delta_TRL = 0
else: 
    print ("Invalid concept number. Please enter a number between 1 and 6.")

#print(f"Development learning factor = {L_d}")
#print(f"Profit retention cost reduction factor = {c_p}")


#LISTS
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

if concept_number == "1":
    M = []
elif concept_number == "2":
    M = []
elif concept_number == "3":
    M = []
elif concept_number == "4":
    M = []
elif concept_number == "5":
    M = []
elif concept_number == "6":
    M = []
else: 
    print ("Invalid concept number. Please enter a number between 1 and 6.")
HW = [25, 25, 25, 25, 25, 25, 25, 25, 25, 25, 25, 25, 25, 25, 25, 25, 25, 25, 25, 25, 25] 

def calc_T1 (a, b, M):
    T1 = []
    for i in range(len(a)):
        C = a[i] * M**b[i]
        T1.append(C)
    T1_total = np.sum(T1)
    return T1

T1 = calc_T1(a, b, M)

def calc_DEV(T1, M_PA_perc, delta_TRL, HW, L_d, c_p):
    DEV_list = []
    for i in range(len(a)):
        DD = 3*T1[i]
        DM = 0.3 * T1[i]
        EM = 1.3 * T1[i]
        PFM = 1.5 * T1[i]
        DD_prime = DD + delta_TRL
        FM1 = T1 - M_PA_perc
        ENG = DD_prime * FM1
        STH = DM + EM + PFM
        MAIT = FM1*STH*L_d*HW[i]
        DEV = c_p*((ENG+(MAIT+ENG)*M_PA_perc)+(FM1*STH*L_d*HW[i]))
        DEV_list.append(DEV)
        DEV_total = np.sum(DEV_list)
    return DEV_total












print ("Development cost (engineering cost) = ", calc_DEV(T1, M_PA_perc, delta_TRL, HW, L_d, c_p))

